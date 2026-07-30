suppressPackageStartupMessages({
  library(io)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggnewscale)
  library(grid)
  library(scales)
  library(ggrepel)
  library(patchwork)
  library(devtools)
  library(parallel)
  library(purrr)
  library(GenomicRanges)
  library(cngpld)
  library(rjson)
  library(optparse)
})

# Parse command line arguments #################################################
option_list = list(
	## REQUIRED
	make_option(c("--case_regions"), type="character", default=NULL, help="Path to the case cngpld_sig-regions file (REQUIRED).", metavar="FILE"),
	make_option(c("--control_regions"), type="character", default=NULL, help="Path to the control cngpld_sig-regions file (REQUIRED).", metavar="FILE"),
	make_option(c("--outdir"), type="character", default=NULL, help="Output directory for all results and cache (REQUIRED).", metavar="DIR"),

	## OPTIONAL
  make_option(c("--case_tag"), type="character", default="case", help="Tag for the case cohort output file names (default: %default).", metavar="STRING"),
  make_option(c("--control_tag"), type="character", default="control", help="Tag for the control cohort output file names (default: %default).", metavar="STRING"),
  make_option(c("--case_label"), type="character", default="case", help="Label for the case cohort (default: %default).", metavar="STRING"),
  make_option(c("--control_label"), type="character", default="control", help="Label for the control cohort (default: %default).", metavar="STRING"),
	make_option(c("--genome"), type="character", default="hg19", help="Reference genome version (default: %default).", metavar="STRING"),
	# Annotation and Score Thresholds
  make_option(c("--clip_fdr"), type="numeric", default=0,  help="FDR value of volcano plot to clip (default: %default).", metavar="NUM"),
  make_option(c("--clip_fc"), type="numeric", default=NA,  help="FC value of volcano plot to clip (default: %default).", metavar="NUM"),
	make_option(c("--fdr_threshold"), type="numeric", default=0.1,  help="FDR threshold for annotation scoring (default: %default).", metavar="NUM"),
	make_option(c("--fc_threshold"), type="numeric", default=1.15,  help="Fold-change threshold for annotation scoring (default: %default).", metavar="NUM"),
	make_option(c("--frac_patients_threshold"), type="numeric", default=0.01,  help="Minimum fraction of patients required to retain an interval (default: %default).", metavar="NUM"),
	make_option(c("--score_threshold"), type="numeric", default=0.5, help="Annotation confidence threshold (default: %default).", metavar="NUM"),
	make_option(c("--min_seg_size"), type="integer", default=5e4, help="Minimum segment size (base pairs) for significance (default: %default).", metavar="INT"),
	make_option(c("--n_obs_threshold"), type="integer", default=5, help="Minimum number of observations (samples) required for significance (default: %default).", metavar="INT"),
	make_option(c("--font_size_scale"), type="numeric", default=3.3, help="Scale for font size of annotations (default: %default).", metavar="NUM"),
  # Figure geometry. The canvas grows beyond fig_width/fig_height as needed to fit labels
  # that are drawn outside the panels; the panels themselves keep the size set here.
  make_option(c("--fig_width"), type="numeric", default=8, help="Width (inches) of the plotting canvas (default: %default).", metavar="NUM"),
  make_option(c("--fig_height"), type="numeric", default=7, help="Height (inches) of the plotting canvas (default: %default).", metavar="NUM"),
  make_option(c("--fig_margin"), type="numeric", default=0.02, help="Extra padding (inches) around the figure (default: %default).", metavar="NUM"),
  # Additional files for driver gene annotations
  make_option(c("--drivers_amp_file"), type="character", default=NA, help="Path to drivers.amp.txt file for driver gene annotation (default: %default).", metavar="FILE"),
  make_option(c("--drivers_del_file"), type="character", default=NA, help="Path to drivers.del.txt file for driver gene annotation (default: %default).", metavar="FILE"),
  make_option(c("--gencode"), type="character", default=NULL, help="Path to ABSOLUTE gencode genes for driver gene annotation (default: %default).", metavar="FILE"),
  make_option(c("--chr_arms"), type="character", default=NULL, help="Path to ABSOLUTE chromosome arms data for driver gene annotation (default: %default).", metavar="FILE")
)
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

# Explicitly check for required arguments (those defined with default=NULL)
required_args <- c("case_regions", "control_regions", "outdir")
missing_args <- required_args[sapply(opt[required_args], is.null)]
if (length(missing_args) > 0) {
  stop(paste("Missing arguments:", paste(paste0("--", missing_args), collapse=", ")))
}

get_script_dir <- function() {
  # Check for 'Rscript' execution (most reliable for command-line tools)
  cmd_args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmd_args)
  if (length(match) > 0) {
    path <- sub(needle, "", cmd_args[match])
    return(dirname(normalizePath(path)))
  }
  return(".")
}
project_root <- get_script_dir()
if (is.null(opt$gencode)) {
  opt$gencode <- file.path(project_root, "data-raw", opt$genome, paste0("gencode.", opt$genome, ".genes.RData"))
}
if (is.null(opt$chr_arms)) {
  opt$chr_arms <- file.path(project_root, "data-raw", opt$genome, paste0(opt$genome, "_ChrArmsDat.RData"))
}

cat("Input options:\n")
for (name in names(opt)) {
  cat(sprintf("  %s: %s\n", name, toString(opt[[name]])))
}
case_regions <- opt$case_regions
case_tag <- opt$case_tag
case_label <- opt$case_label
control_regions <- opt$control_regions
control_tag <- opt$control_tag
control_label <- opt$control_label
outdir <- opt$outdir
genome <- opt$genome
fdr_threshold <- opt$fdr_threshold
clip_fdr <- opt$clip_fdr
clip_fc <- opt$clip_fc
if (!is.na(clip_fc) && clip_fc <= 1) {
  stop("--clip_fc must be greater than 1 (fold change is clipped symmetrically to [1/clip_fc, clip_fc]).")
}
fc_threshold <- opt$fc_threshold
frac_patients_threshold <- opt$frac_patients_threshold
score_threshold <- opt$score_threshold
min_seg_size <- opt$min_seg_size
n_obs_threshold <- opt$n_obs_threshold
font_size_scale <- opt$font_size_scale
fig_width <- opt$fig_width
fig_height <- opt$fig_height
fig_margin <- opt$fig_margin
drivers_amp_file <- opt$drivers_amp_file
drivers_del_file <- opt$drivers_del_file
gencode_file <- opt$gencode
chr_arms_file <- opt$chr_arms

dir.create(paste0(outdir), showWarnings = FALSE, recursive = TRUE)

# Load Data ###################################################################

# Exit cleanly (no plot file) when there is nothing plottable: an empty cohort is a
# legitimate result, not an error, and callers should not have to handle a crash.
exit_no_plot <- function(...) {
  cat(paste0(..., "\nNo volcano plot produced.\n"))
  quit(save = "no", status = 0)
}

read_regions <- function(path, tag) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    cat(sprintf("Region file is missing or empty: %s\n", path))
    return(NULL)
  }
  regions <- io::qread(path)
  if (is.null(regions) || nrow(regions) == 0) {
    cat(sprintf("No regions in: %s\n", path))
    return(NULL)
  }
  regions$cohort <- tag
  return(mutate(regions, chromosome = as.character(chromosome)))
}

regions.case <- read_regions(case_regions, case_tag)
regions.control <- read_regions(control_regions, control_tag)

regions.all <- bind_rows(regions.case, regions.control)
if (nrow(regions.all) == 0) {
  exit_no_plot("Neither cohort contributed any region.")
}

regions.all <- regions.all %>%
  filter(type %in% c("Amp", "Del")) %>%
  mutate(group = recode(type, Amp = "Amplification", Del = "Deletion")) %>%
  mutate(
    fc = exp(ldiff),
    # fc is drawn on a log axis and centred on 1, so clip it symmetrically in log
    # space: clip_fc = 5 constrains fc to [1/5, 5]. NA (default) means no clipping.
    fc = if (is.na(clip_fc)) fc else pmin(pmax(fc, 1 / clip_fc), clip_fc),
    fdr = pmax(fdr, clip_fdr),
    chr = paste0(chromosome, arm)
  ) %>%
  arrange(score)  # draw lower score points first.

# fc and fdr are drawn on log axes, so non-finite/non-positive values cannot be placed.
n_undrawable <- sum(!(is.finite(regions.all$fc) & regions.all$fc > 0 & is.finite(regions.all$fdr)))
if (n_undrawable > 0) {
  cat(sprintf("Dropping %d region(s) with non-finite fold change or FDR.\n", n_undrawable))
  regions.all <- regions.all %>% filter(is.finite(fc), fc > 0, is.finite(fdr))
}
if (nrow(regions.all) == 0) {
  exit_no_plot("No Amp/Del region with a finite fold change and FDR.")
}
regions.all$group <- factor(regions.all$group, levels = c("Amplification", "Deletion", "non-significant"))


# Annotate results ############################################################

x_abslog <- function(f, t = exp(1)) { return(abs(log(f) / log(t)) )}
invx <- function(y, t = exp(1), base = 10) { return( base^(-x_abslog(f=t, t=base) * y) ) }  # only choosing the second branch!

# evidence / significance / signal score defined in run.cngpld.R
beta = -log(fdr_threshold)
sig <- function(x) {
  return( 1 - exp(-beta * x) )
}  # 1 - exp(beta * log(f) / log(t)) = 1 - f ^ (beta / log(t))
invsig <- function(y) {
  return( -log1p(-y) / beta )
}
def_score <- function(fdr, fc) {
  return(sig(x_abslog(fdr, t=fdr_threshold)) * sig(x_abslog(fc, t=fc_threshold)) )
}
fdr_sig_threshold_from_score <- function(score) {
  return( function(fc) { return( invx(invsig( score / sig(x_abslog(fc, t = fc_threshold)) ), t = fdr_threshold, base = 10) ) } )
}

# Already set at the end of cngpld run:
# regions.all$score <- sig(x_abslog(regions.all$fdr, t=fdr_threshold)) * sig(x_abslog(regions.all$fc, t=fc_threshold))

# The non-significant group is just for annotation purposes:
idx <- with(
  regions.all,
  frac_patients >= frac_patients_threshold
  & n_obs >= n_obs_threshold
  # & score >= score_threshold  # data is colored by score
  & end - start + 1 > min_seg_size  # intervals are padded
)
regions.all$group[!idx] <- "non-significant"


regions.all$gene <- NA

if (file.exists(drivers_amp_file) & file.exists(drivers_del_file) & file.exists(chr_arms_file) & file.exists(gencode_file)) {
  cat("Annotating intervals with pre-selected nominated driver genes...\n")

  # Annotate intervals with pre-selected nominated driver genes
  lines <- readLines(drivers_amp_file)
  gene_lines <- lines[!grepl("^#", lines) & nzchar(lines)]
  amps <- sapply(strsplit(gene_lines, "\\s+"), `[`, 1)

  lines <- readLines(drivers_del_file)
  gene_lines <- lines[!grepl("^#", lines) & nzchar(lines)]
  dels <- sapply(strsplit(gene_lines, "\\s+"), `[`, 1)

  load(gencode_file) # loads gencode
  load(chr_arms_file) # loads chr.arms.dat
  gencode.dt <- as.data.table(gencode)
  chr.arms.dt <- as.data.table(chr.arms.dat)
  setnames(gencode.dt, c("Chr", "Start", "End"), c("chr", "start", "end"))
  setnames(chr.arms.dt, c("Start.bp", "End.bp"), c("start", "end"))
  gencode.dt[, chr := as.character(chr)]
  chr.arms.dt[, chr := as.character(chr)]
  chr.arms.dt[, arm := rownames(chr.arms.dat)]
  setkey(chr.arms.dt, chr, start, end)
  setkey(gencode.dt, chr, start, end)
  genes <- foverlaps(gencode.dt, chr.arms.dt, type = "within", nomatch = NA)
  genes <- genes[, .(HGNC, Chr = chr, Start = i.start, End = i.end, Arm = arm)]

  pad <- 2 * 1e6
  annotate_drivers <- function(regions, symbols, region_type) {
    for (g in symbols) {
      loci <- genes[genes$HGNC == g, ]
      if (nrow(loci) == 0) {
        next  # symbol is not in the gencode reference
      }
      # A symbol can sit at more than one locus - PPP2R3B, for instance, is annotated on
      # both Xp and Yp. Testing the region against all loci at once recycles the shorter
      # vector: it mis-assigns rows, and errors outright when there is a single region.
      near <- Reduce(`|`, lapply(seq_len(nrow(loci)), function(k) {
        regions$chr == loci$Arm[k] &
          regions$start - pad <= loci$End[k] &
          regions$end + pad >= loci$Start[k]
      }))
      hits <- (
        near
        & regions$type == region_type
        & regions$group != "non-significant"
        & regions$score > score_threshold
      )
      hits[is.na(hits)] <- FALSE
      regions$gene[hits] <- ifelse(
        is.na(regions$gene[hits]) | regions$gene[hits] == "" | regions$gene[hits] == regions$chr[hits],
        g,
        paste(regions$gene[hits], g, sep = ",")
      )
    }
    return(regions)
  }

  regions.all <- annotate_drivers(regions.all, amps, "Amp")
  regions.all <- annotate_drivers(regions.all, dels, "Del")

} else {
  cat("Annotating intervals with chromosome arm/region...\n")

  # iterate over all rows in regions.case and annotate the chromosome arm/region:
  for (i in seq_len(nrow(regions.all))) {
    if (
      (regions.all$score[i] >= score_threshold)
        & (regions.all$group[i] != "non-significant")
        & is.na(regions.all$gene[i])
    ) {
      regions.all$gene[i] <- paste0(regions.all$chromosome[i], regions.all$arm[i], ":", floor(regions.all$start[i] / 1e6), "-", ceiling(regions.all$end[i] / 1e6), "M")
      # regions.all$gene[i] <- paste0(regions.all$chromosome[i], regions.all$arm[i])
    }
  }

}


# Plot results ################################################################

revlog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}

log_minor_break = function (base = 10, ...) {
  function(x) {
    minx         = floor(min(log(x, base=base), na.rm=T)) - 1
    maxx         = ceiling(max(log(x, base=base), na.rm=T)) + 1
    n_major      = maxx - minx + 1
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks =
      rep(log(seq(1, base-1, by=1), base=base), times = n_major) +
      rep(major_breaks, each = base-1)
    return(base^(minor_breaks))
  }
}

log_minor_break_dense <- function(base = 2, n = 4) {
  force(base); force(n)

  function(x) {
    x_log <- log(x, base = base)
    minx <- floor(min(x_log, na.rm = TRUE)) - 1
    maxx <- ceiling(max(x_log, na.rm = TRUE)) + 1
    major_breaks <- seq(minx, maxx, by = 1)

    # Create multipliers between 1 and base, evenly spaced
    minor_factors <- seq(1, base, length.out = n + 2)[-c(1, n + 2)]  # drop endpoints

    # Get all minor positions
    minor_breaks <- outer(minor_factors, base^major_breaks, `*`)
    sort(as.vector(minor_breaks))
  }
}

mix_to_grey_vec <- function(col, w, grey = "#BFBFBF", darker = FALSE) {
  # col: vector of base colors (e.g. group colors)
  # w:   numeric in [0,1]
  # returns: vector of mixed colors
  if (length(col) == 0) {
    return(character(0))  # mapply() would return an empty list, not a colour vector
  }
  w <- ifelse(is.finite(w), pmin(pmax(w, 0), 1), 0)
  mapply(function(cc, ww) {
    ramp <- scales::colour_ramp(c(grey, cc))
    rramp <- scales::colour_ramp(c("#000000", ramp(ww)))
    return(ifelse(darker, rramp(0.66), ramp(ww)))
  }, col, w)
}

text_width_in <- function(label, size) {
  # size is a ggplot text size (mm); .pt converts it to a font size in points.
  return(grid::convertWidth(
    grid::grobWidth(grid::textGrob(label, gp = grid::gpar(fontsize = size * .pt))),
    "in", valueOnly = TRUE
  ))
}

# The cohort labels are drawn outside the panel (coord_cartesian(clip = "off")), so the
# figure edge - not the panel - is what cuts them off. Measure how far they overhang and
# return the outer margins that make room for them. The canvas is then grown by exactly
# the extra margin, so the panels keep the size, and hence the aspect ratio, they have on
# the base canvas.
fit_outer_margin <- function(plot, width, height, pad = 0) {
  grDevices::pdf(NULL, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  gt <- suppressWarnings(patchwork::patchworkGrob(plot))  # the real render reports these
  w <- grid::convertWidth(gt$widths, "in", valueOnly = TRUE)
  h <- grid::convertHeight(gt$heights, "in", valueOnly = TRUE)
  base <- c(t = h[1], r = w[length(w)], b = h[length(h)], l = w[1])  # patchwork's own margin

  # The panel column is the only "null" width; null units convert to 0, so whatever the
  # canvas has left over after the fixed columns is the panel.
  panel_col <- which(grid::unitType(gt$widths) == "null")
  cols <- seq_along(w)
  panel_w <- width - sum(w)
  inner_l <- sum(w[cols < min(panel_col) & cols != 1])              # y axis title and text
  inner_r <- sum(w[cols > max(panel_col) & cols != length(w)])      # legend

  # Panel-relative x of each label anchor (scale_x_continuous expands limits by 5%).
  lim <- log2(c(xmin, xmax))
  expand <- 0.05 * diff(lim)
  rel <- (log2(c(x_left, x_right)) - (lim[1] - expand)) / (diff(lim) + 2 * expand)

  # control_label is right-aligned at x_left, case_label left-aligned at x_right.
  over_l <- text_width_in(control_label, cohort_label_size) - rel[1] * panel_w - inner_l
  over_r <- text_width_in(case_label, cohort_label_size) - (1 - rel[2]) * panel_w - inner_r

  margin <- pmax(base, c(t = 0, r = over_r, b = 0, l = over_l)) + pad
  return(list(margin = margin, grow = margin - base))
}

hits_to_plot <- (regions.all$fdr < 5 * fdr_threshold) & (regions.all$group != "non-significant")
if (!any(hits_to_plot)) {
  exit_no_plot("No significant hits to plot.")
}
# The volcano is centred on fc = 1: the reference line, the cohort labels and every
# non-significant point live around it. If all hits sit on one side (e.g. a single
# strong amplification), a range derived from the hits alone excludes 1 entirely.
# Every value in the other panel is then censored, its x scale trains to c(Inf, -Inf),
# and the duplicated x axis dies in ggplot2's mono_test(). Keep 1 inside the range.
fc_hits <- regions.all$fc[hits_to_plot]
xmin <- 1 / 1.1 * min(1, fc_hits)
xmax <- 1.1     * max(1, fc_hits)

# Anchors of the two cohort labels above the top panel: midpoints in log2 space for the
# left/right halves. Needed outside make_plot() to size the figure margins.
cohort_label_size <- 4.5
log2_mid <- min(abs(log2(xmin) / 5), abs(log2(xmax) / 5))
x_left  <- 2^(-log2_mid)
x_right <- 2^log2_mid

plot_vulcano_combined <- function(amp_regions, del_regions, suffix = "both") {
  group_colors <- c(
    "Amplification" = "#DC4A4B",
    "Deletion" = "#2080C2",
    "non-significant" = "#BFBFBF"
  )

  # Panel floor: fall back to a default when a panel has no significant region, and
  # ignore fdr == 0 (possible with clip_fdr = 0), which would push the log axis to -Inf.
  fdr_axis_min <- function(fdr, default = 0.1) {
    fdr <- fdr[is.finite(fdr) & fdr > 0]
    if (length(fdr) == 0) {
      return(default)
    }
    return(min(fdr) / 1.2)
  }
  sig_amp_fdr = amp_regions$fdr[amp_regions$group != "non-significant"]
  sig_del_fdr = del_regions$fdr[del_regions$group != "non-significant"]
  amp_ymin = fdr_axis_min(sig_amp_fdr)
  del_ymin = fdr_axis_min(sig_del_fdr)
  ymax = 1

  height_ratio = (log(amp_ymin) - log(ymax)) / (log(del_ymin) - log(ymax))
  if (!is.finite(height_ratio) || height_ratio <= 0) {
    height_ratio = 1
  }

  make_plot <- function(regions, is_top=FALSE) {
    ymin = if (is_top) amp_ymin else del_ymin

    xlim = c(xmin, xmax)
    xbreaks = 2^(floor(log2(xmin)):ceiling(log2(xmax)))
    label_func = function(x) trimws(formatC(x, format = "fg", drop0trailing = TRUE))

    if (is_top) {
      ytrans = revlog_trans(10)
      ylim = c(ymax, ymin)
      scale_x <- scale_x_continuous(
        trans = log_trans(2),
        limits = xlim,
        breaks = xbreaks,
        labels = label_func,
        minor_breaks = log_minor_break_dense(base=2, n=7)
      )
      title <- annotate(
        "text",
        x     = c(x_left, x_right),
        y     = c(ymin, ymin),
        label = c(control_label, case_label),
        hjust = c(1, 0),
        vjust = -1,
        size  = cohort_label_size
      )
    } else {
      ytrans = log_trans(10)
      ylim = c(ymin, ymax)
      scale_x <- scale_x_continuous(
        trans = log_trans(2),
        limits = xlim,
        breaks = xbreaks,
        labels = label_func,
        minor_breaks = log_minor_break_dense(base=2, n=7),
        sec.axis = if (!is_top) dup_axis(name = NULL) else NULL
      )
      title = NULL
    }

    # Precompute mixed colour for this panel
    regions <- regions %>%
      mutate(
        base_col = group_colors[as.character(group)],
        mixed_col = mix_to_grey_vec(base_col, score, grey = group_colors["non-significant"]),
        mixed_text_col = mix_to_grey_vec(base_col, score, grey = group_colors["non-significant"], darker = TRUE)
      )

    # Display score contour lines:
    scores <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999)
    sig_threshold_lines <- map_dfr(scores, function(w) {
      tibble(
        score = w,
        fc    = seq(xmin, xmax, length.out = 1000),
        fdr   = fdr_sig_threshold_from_score(w)(fc)
      )
    })

    strip_trailing_zeros <- function(x, accuracy = 0.001) {
      raw <- scales::label_number(accuracy = accuracy, trim = TRUE)(x)
      return(sub("\\.?0+$", "", raw))   # drop trailing zeros, then optional "."
    }

    sig_threshold_labels <- sig_threshold_lines %>%
      filter(fc > 1, !is.na(fdr), fdr >= ymin, fdr <= ymax) %>%
      group_by(score) %>%
      slice_max(fc, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(label = strip_trailing_zeros(score, accuracy = 0.001))


    dummy_type <- regions |> dplyr::distinct(group) |> dplyr::mutate(fc = NA_real_, fdr = NA_real_, frac_patients = NA_real_)

    # MAKE PLOT

    ggplot() +
      # MID LINE
      geom_vline(xintercept = 1) +

      # # SCORE CONTOUR LINES
      # geom_line(
      #   data = sig_threshold_lines,
      #   aes(x = fc, y = fdr, group = score, alpha = score),
      #   color = "#BFBFBF",
      #   linetype = 1,
      #   linewidth = 0.3,
      #   show.legend = FALSE
      # ) +
      # geom_text(
      #   data = sig_threshold_labels,
      #   aes(x = fc, y = fdr, label = label, alpha = score),
      #   color = "#BFBFBF",
      #   hjust = -0.1,
      #   vjust = 0.5,
      #   size = 1.5,
      #   show.legend = FALSE
      # ) +
      scale_alpha_identity(name = "Scale", guide = "none") +

      # SCATTER PLOT
      geom_point(
        data = regions,
        aes(x = fc, y = fdr, size = frac_patients, fill = I(mixed_col)),
        shape  = 21,
        stroke = 0.1,
        # position = position_jitter(width = 0.02, height = 0.02),
        show.legend = c(size = !is_top, fill = FALSE)
      ) +
      # # LEGEND-ONLY POINTS (for group legend)
      geom_point(
        data = dummy_type,
        aes(x = fc, y = fdr, fill = group),
        shape        = 21,
        colour       = "black",
        size         = 3.5,
        stroke       = 0.3,
        inherit.aes  = FALSE,
        show.legend  = c(size = FALSE, fill = !is_top)
      ) +
      scale_fill_manual(
        name = "Type",
        values = group_colors,
        drop = FALSE,
        guide = guide_legend(
          override.aes = list(
            alpha  = 1,
            shape  = 21,
            colour = "black"
          )
        )
      ) +
      scale_size_continuous(
        name   = if (!is_top) "% patients\nwith event" else NULL,
        range  = c(0.1, 5),
        breaks = if (!is_top) c(0.25, 0.5, 0.75) else waiver(),
        labels = if (!is_top) scales::percent_format(accuracy = 1) else waiver(),
        guide  = if (is_top) "none" else guide_legend(
          override.aes = list(
            alpha  = 1,
            shape  = 21,
            colour = "black"
          )
        )
      ) +

      # ANNOTATIONS
      ggnewscale::new_scale("size") +  # Scale font size separately
      geom_text_repel(
        data = regions,
        show.legend = FALSE,
        nudge_y = 0.0,
        max.overlaps = 50,
        segment.size = 0.3,
        min.segment.length = 0,
        force = 4,
        max.time = 1,
        aes(x = fc, y = fdr, label = gene, size = font_size_scale * score, colour = I(mixed_text_col))
      ) +
      scale_size_identity(guide = "none") +

      # LAYOUT
      theme_classic() +
      theme(
        plot.margin = margin(t = 20, r = 0, b = 0, l = 0),  # space for labels
        axis.ticks.x.bottom = if (!is_top) element_blank() else element_line(),
        axis.text.x.bottom = element_blank(),
        axis.text.x.top = if (!is_top) element_text(size=7.5, margin = margin(b = 5.3, t = -100)) else element_blank(),  # get x tick labels from bottom panel
        axis.line.x.bottom = if (!is_top) element_blank() else element_line(),
        axis.line.y = element_blank(),
        axis.text.y.left = element_text(size=7.5),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.25)
      ) +
      scale_x +
      scale_y_continuous(
        trans = ytrans,
        limits = ylim,
        breaks = log_breaks(base = 10),
        minor_breaks = log_minor_break(base=10),
        labels = label_func,
        expand = c(0, 0),  # No padding
        oob = scales::oob_censor
      ) +
      title +
      coord_cartesian(clip = "off") +   # allow drawing outside panel
      xlab(if (!is_top) paste0("CNV summary statistic contrast") else NULL) +
      ylab(paste("false discovery rate", if (is_top) "(AMP)" else "(DEL)"))
  }

  amp_case_plot <- make_plot(amp_regions, is_top=TRUE)
  del_case_plot <- make_plot(del_regions)

  panels <- (amp_case_plot / del_case_plot) +
    plot_layout(guides = "collect", heights = c(height_ratio, 1))
  legend_theme <- theme(legend.position = "right", legend.justification = "top")

  # Widen the canvas until the cohort labels fit on it, keeping the panels' size fixed.
  fit <- tryCatch(
    fit_outer_margin(panels & legend_theme, fig_width, fig_height, pad = fig_margin),
    error = function(e) {
      cat(sprintf("Could not measure the figure margins (%s); using defaults.\n", conditionMessage(e)))
      return(NULL)
    }
  )
  if (is.null(fit)) {
    final_plot <- panels & legend_theme
    canvas <- c(width = fig_width, height = fig_height)
  } else {
    final_plot <- panels +
      plot_annotation(theme = theme(plot.margin = margin(
        t = fit$margin[["t"]], r = fit$margin[["r"]],
        b = fit$margin[["b"]], l = fit$margin[["l"]], unit = "in"
      ))) &
      legend_theme
    canvas <- c(
      width  = fig_width  + fit$grow[["l"]] + fit$grow[["r"]],
      height = fig_height + fit$grow[["t"]] + fit$grow[["b"]]
    )
  }

  # Render to a scratch file first so a failure never leaves a truncated PDF behind.
  filename <- paste0(outdir, "/vulcano_", case_tag, "-vs-", control_tag, ".pdf")
  tmpfile <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmpfile), add = TRUE)
  ggsave(filename = tmpfile, plot = final_plot, width = canvas[["width"]], height = canvas[["height"]])
  if (!file.copy(tmpfile, filename, overwrite = TRUE)) {
    stop("Could not write ", filename)
  }
  return(filename)
}

amp_regions <- regions.all[regions.all$type == "Amp", ]
del_regions <- regions.all[regions.all$type == "Del", ]
for (panel in c("Amp", "Del")) {
  if (sum(regions.all$type == panel) == 0) {
    cat(sprintf("No %s region to plot - that panel stays empty.\n", panel))
  }
}

saved <- plot_vulcano_combined(amp_regions, del_regions, suffix="both")
cat(paste0("Saved volcano plot to: ", saved, "\n"))
