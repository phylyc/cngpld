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
gencode_default_path <- file.path(project_root, "data-raw", "hg19", "gencode.hg19.genes.RData")
chr_arms_default_path <- file.path(project_root, "data-raw", "hg19", "hg19_ChrArmsDat.RData")

# Parse command line arguments #################################################
option_list = list(
	## REQUIRED
	make_option(c("--case_regions"), type="character", default=NULL, help="Path to the case cngpld_sig-regions file (REQUIRED).", metavar="FILE"),
	make_option(c("--control_regions"), type="character", default=NULL, help="Path to the control cngpld_sig-regions file (REQUIRED).", metavar="FILE"),
	make_option(c("--outdir"), type="character", default=NULL, help="Output directory for all results and cache (REQUIRED).", metavar="DIR"),

	## OPTIONAL
  make_option(c("--case"), type="character", default="case", help="Label for the case cohort (default: %default).", metavar="STRING"),
  make_option(c("--control"), type="character", default="control", help="Label for the control cohort (default: %default).", metavar="STRING"),
	make_option(c("--genome"), type="character", default="hg19", help="Reference genome version (default: %default).", metavar="STRING"),
	# Annotation and Score Thresholds
	make_option(c("--fdr_threshold"), type="numeric", default=0.1,  help="FDR threshold for annotation scoring (default: %default).", metavar="NUM"),
	make_option(c("--fc_threshold"), type="numeric", default=1.15,  help="Fold-change threshold for annotation scoring (default: %default).", metavar="NUM"),
	make_option(c("--frac_patients_threshold"), type="numeric", default=0.01,  help="Minimum fraction of patients required to retain an interval (default: %default).", metavar="NUM"),
	make_option(c("--score_threshold"), type="numeric", default=0.5, help="Annotation confidence threshold (default: %default).", metavar="NUM"),
	make_option(c("--min_seg_size"), type="integer", default=5e4, help="Minimum segment size (base pairs) for significance (default: %default).", metavar="INT"),
	make_option(c("--n_obs_threshold"), type="integer", default=5, help="Minimum number of observations (samples) required for significance (default: %default).", metavar="INT"),
  # Additional files for driver gene annotations
  make_option(c("--drivers_amp_file"), type="character", default=NA, help="Path to drivers.amp.txt file for driver gene annotation (default: %default).", metavar="FILE"),
  make_option(c("--drivers_del_file"), type="character", default=NA, help="Path to drivers.del.txt file for driver gene annotation (default: %default).", metavar="FILE"),
  make_option(c("--gencode"), type="character", default=gencode_default_path, help="Path to ABSOLUTE gencode genes for driver gene annotation (default: %default).", metavar="FILE"),
  make_option(c("--chr_arms"), type="character", default=chr_arms_default_path, help="Path to ABSOLUTE chromosome arms data for driver gene annotation (default: %default).", metavar="FILE")
)
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

# Explicitly check for required arguments (those defined with default=NULL)
required_args <- c("case_regions", "control_regions", "outdir")
missing_args <- required_args[sapply(opt[required_args], is.null)]
if (length(missing_args) > 0) {
  stop(paste("Missing arguments:", paste(paste0("--", missing_args), collapse=", ")))
}

print(opt) # Print all parsed options for monitoring
case_regions <- opt$case_regions
case <- opt$case
control_regions <- opt$control_regions
control <- opt$control
outdir <- opt$outdir
genome <- opt$genome
fdr_threshold <- opt$fdr_threshold
fc_threshold <- opt$fc_threshold
frac_patients_threshold <- opt$frac_patients_threshold
score_threshold <- opt$score_threshold
min_seg_size <- opt$min_seg_size
n_obs_threshold <- opt$n_obs_threshold
drivers_amp_file <- opt$drivers_amp_file
drivers_del_file <- opt$drivers_del_file
gencode_file <- opt$gencode
chr_arms_file <- opt$chr_arms

dir.create(paste0(outdir), showWarnings = FALSE, recursive = TRUE)

# # Function to extract file label from path
# get_file_label <- function(filepath) {
# 	if (is.na(filepath)) {
# 		return(NA_character_)
# 	}
# 	base_name <- basename(filepath)
# 	# Removes everything from the first dot to the end
# 	label <- sub("\\..*$", "", base_name)
# 	return(label)
# }

# case <- get_file_label(case_regions)
# case <- sub("^cngpld_sig-regions_", "", case)
# control <- get_file_label(control_regions)
# control <- sub("^cngpld_sig-regions_", "", control)
# print(case)
# print(control)
# fits.fn <- paste0(outdir, "/cngpld/", case, "-vs-", control, ".rds")

# Load Data ###################################################################
regions.case <- io::qread(case_regions)
regions.control <- io::qread(control_regions)
regions.case$cohort <- case
regions.control$cohort <- control

min_fdr <- 1e-3
regions.all <- rbind(
  data.frame(filter(regions.case, type == "Amp"), group = "Amplification"),
  data.frame(filter(regions.case, type == "Del"), group = "Deletion"),
  data.frame(filter(regions.control, type == "Amp"), group = "Amplification"),
  data.frame(filter(regions.control, type == "Del"), group = "Deletion")
) %>%
  mutate(
    fc = exp(ldiff),
    fdr = pmax(fdr, min_fdr),
    chr = paste0(chromosome, arm)
  ) %>%
  arrange(score)  # draw lower score points first.
regions.all$group <- factor(regions.all$group, levels = c("Amplification", "Deletion", "non-significant"))


# Annotate results ############################################################

x_abslog <- function(f, t = exp(1)) { return(abs(log(f) / log(t)) )}
invx <- function(y, t = exp(1), base = 10) { return( base^(-x_abslog(f=t, t=base) * y) ) }  # only choosing the second branch!

# evidence / significance / signal score defined in run.cngpld.R
beta = -log(fdr_threshold)
sig <- function(x) { return( 1 - exp(-beta * x) ) }  # 1 - exp(beta * log(f) / log(t)) = 1 - f ^ (beta / log(t))
invsig <- function(y) { return( -log1p(-y) / beta ) }
def_score <- function(fdr, fc) { return(sig(x_abslog(fdr, t=fdr_threshold)) * sig(x_abslog(fc, t=fc_threshold)) ) }
fdr_sig_threshold_from_score <- function(score) { return( function(fc) { return( invx(invsig( score / sig(x_abslog(fc, t = fc_threshold)) ), t = fdr_threshold, base = 10) ) } ) }

# Already set at the end of cngpld run:
# regions.all$score <- sig(x_abslog(regions.all$fdr, t=fdr_threshold)) * sig(x_abslog(regions.all$fc, t=fc_threshold))

idx <- with(
  regions.all,
  frac_patients >= frac_patients_threshold
  & n_obs >= n_obs_threshold
  # & score >= score_threshold  # data is colored by score
  & end - start + 1 > min_seg_size  # intervals are padded
)

if ("is_significant" %in% colnames(regions.all)) {
  # If it is NOT "S" (meaning it is "NS"), exclude it from significance
  idx <- idx & (regions.all$is_significant == "S")
}

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

  suppressWarnings({
    pad <- 1.5 * 1e6
    for (g in amps) {
      gene <- genes[genes$HGNC == g, ]
      hits <- (
        regions.all$chr == gene$Arm
        & regions.all$start - pad <= gene$End
        & regions.all$end + pad >= gene$Start
        & regions.all$type == "Amp"
        & regions.all$group != "non-significant"
        & regions.all$score > score_threshold
      )
      regions.all$gene[hits] <- ifelse(
        is.na(regions.all$gene[hits]) | regions.all$gene[hits] == "" | regions.all$gene[hits] == regions.all$chr[hits],
        g,
        paste(regions.all$gene[hits], g, sep = ",")
      )
    }

    for (g in dels) {
      gene <- genes[genes$HGNC == g, ]
      hits <- (
        regions.all$chr == gene$Arm
        & regions.all$start - pad <= gene$End
        & regions.all$end + pad >= gene$Start
        & regions.all$type == "Del"
        & regions.all$group != "non-significant"
        & regions.all$score > score_threshold
      )
      regions.all$gene[hits] <- ifelse(
        is.na(regions.all$gene[hits]) | regions.all$gene[hits] == "" | regions.all$gene[hits] == regions.all$chr[hits],
        g,
        paste(regions.all$gene[hits], g, sep = ",")
      )
    }
  })

} else {
  cat("Annotating intervals with chromosome arm/region...\n")

  # iterate over all rows in regions.case and annotate the chromosome arm/region:
  for (i in 1:nrow(regions.all)) {
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
  mapply(function(cc, ww) {
    ramp <- scales::colour_ramp(c(grey, cc))
    rramp <- scales::colour_ramp(c("#000000", ramp(ww)))
    return(ifelse(darker, rramp(0.66), ramp(ww)))
  }, col, w)
}

xmin = 1 / 1.1 * min(regions.all$fc[regions.all$fdr < 5 * fdr_threshold])
xmax = 1.1     * max(regions.all$fc[regions.all$fdr < 5 * fdr_threshold])
# xmin = 0.45
# xmax = 2

plot_vulcano_combined <- function(amp_regions, del_regions, suffix = "both") {
  group_colors <- c(
    "Amplification" = "#DC4A4B",
    "Deletion" = "#2080C2",
    "non-significant" = "#BFBFBF"
  )

  amp_ymin = min(amp_regions$fdr[amp_regions$group != "non-significant"]) / 1.2
  del_ymin = min(del_regions$fdr[del_regions$group != "non-significant"]) / 1.2
  ymax = 1

  height_ratio = (log(amp_ymin) - log(ymax)) / (log(del_ymin) - log(ymax))

  make_plot <- function(regions, is_top=FALSE) {
    ymin = if (is_top) amp_ymin else del_ymin

    xlim = c(xmin, xmax)
    xbreaks = 2^(floor(log2(xmin)):ceiling(log2(xmax)))
    label_func = function(x) trimws(formatC(x, format = "fg", drop0trailing = TRUE))

    # midpoints in log2 space for left/right halves
    log2_xmin <- log2(xmin)
    log2_xmax <- log2(xmax)
    log2_mid_left  <- (log2_xmin + log2(1)) / 5
    log2_mid_right <- (log2(1) + log2_xmax) / 5
    log2_mid <- min(abs(log2_mid_left), abs(log2_mid_right))
    x_left  <- 2^(-log2_mid)
    x_right <- 2^log2_mid

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
        label = c(control, case),
        hjust = c(1, 0),
        vjust = -1,
        size  = 4
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

      # SCORE CONTOUR LINES
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
        name   = if (!is_top) "Fraction of patients\nwith event" else NULL,
        range  = c(0.1, 4),
        breaks = if (!is_top) c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) else waiver(),
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
        aes(x = fc, y = fdr, label = gene, size = 3.3 * score, colour = I(mixed_text_col))
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
      xlab(if (!is_top) paste0("CNV summary statistic fold-change") else NULL) +
      ylab(paste("false discovery rate", if (is_top) "(AMP)" else "(DEL)"))
  }

  amp_case_plot <- make_plot(amp_regions, is_top=TRUE)
  del_case_plot <- make_plot(del_regions)

  final_plot <- (amp_case_plot / del_case_plot) +
    plot_layout(guides = "collect", heights = c(height_ratio, 1)) &
    theme(legend.position = "right")

  ggsave(
    filename = paste0(outdir, "/vulcano_", case, "-vs-", control, ".pdf"),
    plot = final_plot,
    width = 8, height = 7
  )
}

plot_vulcano_combined(regions.all[regions.all$type == "Amp", ], regions.all[regions.all$type == "Del", ], suffix="both")
cat(paste0("Saved volcano plot to: ", outdir, "/vulcano_", case, "-vs-", control, ".pdf\n"))