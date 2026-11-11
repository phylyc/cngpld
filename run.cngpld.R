suppressPackageStartupMessages({
  library(io)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(patchwork)
  library(parallel)
  library(purrr)
  library(GenomicRanges)
  library(cngpld)
  library(gpldiff)
  library(rjson)
  library(optparse)
})

# Parse command line arguments #################################################
option_list = list(
	## REQUIRED
	make_option(c("--case_file"), type="character", default=NULL, help="Path to the case segmentation file (REQUIRED).", metavar="FILE"),
	make_option(c("--control_file"), type="character", default=NULL, help="Path to the control segmentation file (REQUIRED).", metavar="FILE"),
	make_option(c("--outdir"), type="character", default=NULL, help="Output directory for all results and cache (REQUIRED).", metavar="DIR"),

	## OPTIONAL
  make_option(c("--case_tag"), type="character", default="case", help="Tag for the case cohort output file names (default: %default).", metavar="STRING"),
  make_option(c("--control_tag"), type="character", default="control", help="Tag for the control cohort output file names (default: %default).", metavar="STRING"),
  make_option(c("--case_label"), type="character", default="case", help="Label for the case cohort (default: %default).", metavar="STRING"),
  make_option(c("--control_label"), type="character", default="control", help="Label for the control cohort (default: %default).", metavar="STRING"),
	make_option(c("--genome"), type="character", default="hg19", help="Reference genome version (default: %default).", metavar="STRING"),
	make_option(c("--use_cache"), action="store_true", default=FALSE, help="Flag to enable caching of the model (default: FALSE)."),
	make_option(c("--is_paired"), action="store_true", default=FALSE, help="Whether case and control cohorts are paired samples from the same patients (default: FALSE)."),
	# Statistical & CNV Thresholds
	make_option(c("--lodds_cut"), type="numeric", default=NULL, help="Legacy single log-odds threshold; if supplied, disables peak decomposition.", metavar="NUM"),
	make_option(c("--lodds_support"), type="numeric", default=3, help="Log posterior odds threshold defining candidate support regions (default: %default).", metavar="NUM"),
	make_option(c("--lodds_core"), type="numeric", default=5, help="Higher log posterior odds threshold defining peak cores (default: %default).", metavar="NUM"),
	make_option(c("--peak_prominence"), type="numeric", default=1.5, help="Minimum log-odds drop between adjacent cores required to split them (default: %default).", metavar="NUM"),
	make_option(c("--cn_res"), type="numeric", default=100, help="Aggregated log copy-ratio signal resolution 1/cn.res (default: %default).", metavar="NUM"),
	make_option(c("--min_tCR"), type="numeric", default=0.3, help="Absolute copy-ratio threshold for CNV summary stats (default: %default).", metavar="NUM"),
	make_option(c("--max_tCR"), type="numeric", default=3.5, help="Cap maximum total copy-ratio to avoid differential signal being dominated by differences in amount of cDNA amplifications (default: %default).", metavar="NUM"),
	make_option(c("--min_nprobes"), type="integer", default=4,  help="Minimum number of probes supporting a segment (default: %default).", metavar="INT"),
	make_option(c("--sigma2_concentration"), type="numeric", default=0.1, help="Prior concentration of sigma^2 (default: %default).", metavar="NUM"),
	make_option(c("--weight_by_cohort_size"), action="store_true", default=FALSE, help="Weight uncertainty by cohort size (default: FALSE)."),
	# Annotation and Score Thresholds
	make_option(c("--fdr_threshold"), type="numeric", default=0.1,  help="FDR threshold for annotation scoring (default: %default).", metavar="NUM"),
	make_option(c("--fc_threshold"), type="numeric", default=1.15,  help="Fold-change threshold for annotation scoring (default: %default).", metavar="NUM"),
	make_option(c("--frac_patients_threshold"), type="numeric", default=0.1, help="Minimum fraction of patients required to retain an interval (default: %default).", metavar="NUM"),
	make_option(c("--delta_frac_patients_threshold"), type="numeric", default=0.00, help="Minimum increase in the fraction of patients carrying the signal in the enriched cohort relative to the other cohort. Vetoes latent differences that are boundary artifacts of the GP (default: %default).", metavar="NUM"),
	make_option(c("--score_threshold"), type="numeric", default=0.68, help="Annotation confidence threshold (default: %default).", metavar="NUM"),
	make_option(c("--min_seg_size"), type="integer", default=5e4, help="Minimum segment size (base pairs) for significance (default: %default).", metavar="INT"),
	make_option(c("--n_obs_threshold"), type="integer", default=5, help="Minimum number of observations (probes) required for significance (default: %default).", metavar="INT"),
  make_option(c("--n_cores"), type="integer", default=1, help="Number of CPU cores to use (default: %default).", metavar="INT")
)
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

# Explicitly check for required arguments (those defined with default=NULL)
required_args <- c("case_file", "control_file", "outdir")
missing_args <- required_args[sapply(opt[required_args], is.null)]
if (length(missing_args) > 0) {
  stop(paste("Missing arguments:", paste(paste0("--", missing_args), collapse=", ")))
}

cat("Input options:\n")
for (name in names(opt)) {
  cat(sprintf("  %s: %s\n", name, toString(opt[[name]])))
}
case_file <- opt$case_file
case_tag <- opt$case_tag
control_file <- opt$control_file
control_tag <- opt$control_tag
outdir <- opt$outdir
genome <- opt$genome
use_cache <- opt$use_cache
is_paired <- opt$is_paired
lodds_cut <- opt$lodds_cut
lodds_support <- opt$lodds_support
lodds_core <- opt$lodds_core
peak_prominence <- opt$peak_prominence
if (!is.null(lodds_cut)) {
  # Backward-compatible single-threshold behavior.
  lodds_support <- lodds_cut
  lodds_core <- NULL
}
cn_res <- opt$cn_res
min_tCR <- opt$min_tCR
max_tCR <- opt$max_tCR
min_nprobes <- opt$min_nprobes
sigma2_concentration <- opt$sigma2_concentration
weight_by_cohort_size <- opt$weight_by_cohort_size
fdr_threshold <- opt$fdr_threshold
fc_threshold <- opt$fc_threshold
frac_patients_threshold <- opt$frac_patients_threshold
delta_frac_patients_threshold <- opt$delta_frac_patients_threshold
score_threshold <- opt$score_threshold
min_seg_size <- opt$min_seg_size
n_obs_threshold <- opt$n_obs_threshold
n_cores <- opt$n_cores

dir.create(paste0(outdir), showWarnings = FALSE, recursive = TRUE)

# Run analysis ################################################################
seg.case <- cngpld::read_seg(
  case_file
) %>% mutate( logr = pmin(pmax(-max_tCR, log(2) * logr), max_tCR) ) %>% filter( nprobes >= min_nprobes )

seg.control <- cngpld::read_seg(
  control_file
) %>% mutate( logr = pmin(pmax(-max_tCR, log(2) * logr), max_tCR) ) %>% filter( nprobes >= min_nprobes )

fits.fn <- paste0(outdir, "/", case_tag, "-vs-", control_tag, ".rds")
if (file.exists(fits.fn) & use_cache) {
  cat("Using cached model.")
  fits <- io::qread(fits.fn)

} else {

  hparams <- cngpld::default_hparams()
  sigma2_mode <- 0.1 / (1 + 0.1)  # default mode from settings alpha=0.1 and beta=0.1
  hparams$alpha <- sigma2_concentration
  hparams$beta <- sigma2_mode * (1 + sigma2_concentration)

  options(mc.cores = n_cores)
  fits <- cngpld::compare_segs(
    seg.case,
    seg.control,
    genome = genome,
    cn.cut = min_tCR,  # absolute threshold for copy-number log ratio to be considered in summary statistics
    smooth = TRUE,  # whether to median smooth the copy-number data
    cn.res = cn_res,  # aggregated log copy-ratio signal resolution 1/cn.res
    pair = is_paired,  # whether case and control cohorts are paired samples from the same patients
    hparams = hparams,
    adapt = "none",
    weight.N.ref = 100,
    weight_by_cohort_size = weight_by_cohort_size,
    verbose = 1,
  )

  io::qwrite(fits, fits.fn)
}


# Examine results #############################################################

# remove problematic deletion profiles
fits.orig <- fits
idx <- unlist(lapply(fits$del, function(x) is(x$model, "gpldiff")))
fits$del <- fits$del[idx]

regions.case <- summary(
  fits,
  genome = genome,
  lodds.cut = lodds_support,
  lodds.core = lodds_core,
  peak.prominence = peak_prominence
)
regions.control <- summary(
  fits,
  direction = -1,
  genome = genome,
  lodds.cut = lodds_support,
  lodds.core = lodds_core,
  peak.prominence = peak_prominence
)

# summary() returns NULL for a direction when no significant regions are found.
# Guard both directions so downstream steps still produce (possibly empty) tables.
empty_regions <- data.frame(
  type = character(),
  chromosome = character(),
  start = numeric(),
  end = numeric(),
  start_idx = integer(),
  end_idx = integer(),
  n_obs = integer(),
  diff = numeric(),
  ldiff = numeric(),
  posterior = numeric(),
  fdr = numeric(),
  stringsAsFactors = FALSE
)
# Prefer a zero-row copy of a non-NULL result (preserves the exact schema),
# otherwise fall back to the canonical empty template.
regions_template <- if (!is.null(regions.case)) {
  regions.case[FALSE, ]
} else if (!is.null(regions.control)) {
  regions.control[FALSE, ]
} else {
  empty_regions
}
if (is.null(regions.case)) {
  regions.case <- regions_template
}
if (is.null(regions.control)) {
  regions.control <- regions_template
}

# Create plots for each chr arm
dir.create(paste0(outdir, "/chr"), showWarnings = FALSE, recursive = TRUE)

plot_region <- function(chr_arm, profile, profile_str = "amp") {
  qdraw({
    with(profile[[chr_arm]], plot(model, data, which = c("response", "latent", "odds"), xlab = "position (Mbp)"))
  }, width = 5, height = 10, file = paste0(outdir, "/chr/", chr_arm, "_", profile_str, ".pdf"))
}

for (arm in (filter(regions.case, type == "Amp")$chromosome)) {
  plot_region(chr_arm = arm, profile = fits$amp, profile_str = "amp")
}
for (arm in (filter(regions.case, type == "Del")$chromosome)) {
  plot_region(chr_arm = arm, profile = fits$del, profile_str = "del")
}


# Annotate and save results ###################################################

# Split chromosome and arm
split_chr_arm <- function(dt, chr_col = "chromosome") {
  dt <- as.data.table(dt)
  # create chromosome (no arm) and arm (p/q/NA)
  dt[, `:=`(
    chromosome = sub("[pq]$", "", get(chr_col)),  # strip trailing p/q
    arm = fifelse(
      grepl("[pq]$", get(chr_col)),
      substr(get(chr_col), nchar(get(chr_col)), nchar(get(chr_col))),
      NA_character_
    )
  )]
  # move "arm" column next to chromosome
  cols <- names(dt)
  new_order <- c("type", "chromosome", "arm", setdiff(cols, c("type", "chromosome", "arm")))
  setcolorder(dt, new_order)
  dt
}

regions.case <- split_chr_arm(regions.case)
regions.control <- split_chr_arm(regions.control)

# Calculate fraction of patients carrying the type of CNV per interval.
# give regions an ID to merge back later
seg.case <- as.data.table(seg.case)
seg.control <- as.data.table(seg.control)

annotate_frac_patients <- function(regions, col, seg, min_tCR) {
  regions <- as.data.table(regions)
  seg <- as.data.table(seg)

  # total number of samples in this cohort
  n_samples <- uniqueN(seg$sample)

  # key by intervals for fast overlap
  setkey(seg,     chromosome, start, end)
  setkey(regions, chromosome, start, end)

  # initialize result; (col) forces evaluation of the variable as column name
  regions[, region_id := .I]
  regions[, (col) := 0]

  if (nrow(regions) == 0 || n_samples == 0) {
    regions[, region_id := NULL]
    return(regions[])
  }

  count_frac <- function(seg_sub, region_type) {
    regions_sub <- regions[type == region_type]
    if (nrow(seg_sub) == 0 || nrow(regions_sub) == 0) return(invisible(NULL))
    setkey(regions_sub, chromosome, start, end)
    hits <- foverlaps(seg_sub, regions_sub, nomatch = 0L)
    if (nrow(hits) == 0) return(invisible(NULL))
    # hits has both seg + region columns; count unique samples per region
    counts <- hits[, .(n_samples_hit = uniqueN(sample)), by = region_id]
    regions[counts, (col) := i.n_samples_hit / n_samples, on = "region_id"]
    invisible(NULL)
  }

  count_frac(seg[logr >= min_tCR], "Amp")   # amplifications
  count_frac(seg[logr <= -min_tCR], "Del")  # deletions

  regions[, region_id := NULL] # drop temp ID
  regions[]
}

annotate_both_cohorts <- function(regions) {
  regions <- annotate_frac_patients(regions, "frac_patients_case", seg.case, min_tCR)
  regions <- annotate_frac_patients(regions, "frac_patients_control", seg.control, min_tCR)
  regions[, delta_frac_patients := frac_patients_case - frac_patients_control]
  regions[]
}

regions.case <- annotate_both_cohorts(regions.case)
regions.control <- annotate_both_cohorts(regions.control)


# The "score" is a measure of confidence in this interval as a characteristic of the difference between those two cohorts.
# Scales with 1 - fdr and with 1 - fc^k for fc < 1 and some power k.
# This score could be used e.g. for downstream pathway or over-representation analysis.
beta = -log(fdr_threshold)
x_abslog <- function(f, t = exp(1)) { 
  return(abs(log(f) / log(t)) )
}
sig <- function(x) { 
  return( 1 - exp(-beta * x) ) 
}  # = 1 - exp(beta * log(f) / log(t)) = 1 - f ^ (beta / log(t))
def_score <- function(fdr, fc) { 
  return(sig(x_abslog(fdr, t=fdr_threshold)) * sig(x_abslog(fc, t=fc_threshold)) ) 
}

# Label regions as significant ("S"), not significant ("NS"), not relevant ("NR"),
# or vetoed ("VETO"). The veto catches latent differences that are artifacts of
# GP boundary effects: a genuine difference must also be carried by a larger fraction
# of the cohort the region is enriched in than by the other cohort.
# The Bayesian FDR is recomputed from the posterior column over the set of retained
# ("S" and "NS") entries only, since "NR" and "VETO" entries are not relevant or are
# vetoed and should not contribute to the multiple-testing burden. The score, which
# separates "S" from "NS", is then derived from this recomputed FDR.
annotate_regions <- function(regions, direction = c("case", "control")) {
  direction <- match.arg(direction)
  regions <- as.data.frame(regions)
  if (direction == "case") {
    frac_own <- regions$frac_patients_case
    delta <- regions$delta_frac_patients
  } else {
    frac_own <- regions$frac_patients_control
    delta <- -regions$delta_frac_patients
  }

  # Retained ("S" or "NS") entries: neither vetoed nor irrelevant.
  # This classification does not depend on the score, so the FDR can be recomputed
  # on this subset without circularity.
  retained <- rep(FALSE, nrow(regions))
  if (nrow(regions) > 0) {
    is_veto <- delta < delta_frac_patients_threshold
    is_relevant <- frac_own >= frac_patients_threshold &
      regions$n_obs >= n_obs_threshold &
      regions$end - regions$start + 1 > min_seg_size
    retained <- !is_veto & is_relevant
    regions$is_significant <- ifelse(is_veto, "VETO", ifelse(is_relevant, "S", "NR"))
  } else {
    regions$is_significant <- character()
  }

  # Recompute the Bayesian FDR from the posterior (= 1 - local fdp) among retained
  # entries, ordered by decreasing posterior.
  if (any(retained)) {
    idx <- which(retained)
    ord <- idx[order(regions$posterior[idx], decreasing = TRUE)]
    regions$fdr[ord] <- gpldiff:::bayesian_fdr(regions$posterior[ord])
  }

  regions %>%
    mutate(
      score = def_score(fdr = fdr, fc = exp(ldiff)),
      IGV = paste0("chr", chromosome, ":", start, "-", end),
      is_significant = ifelse(
        is_significant == "S" & score < score_threshold,
        "NS",
        is_significant
      )
    ) %>%
    arrange(desc(score))
}

regions.case <- annotate_regions(regions.case, direction = "case")
regions.control <- annotate_regions(regions.control, direction = "control")

io::qwrite(regions.case, file.path(outdir, paste0("cngpld_sig-regions.", case_tag, ".tsv")))
io::qwrite(regions.control, file.path(outdir, paste0("cngpld_sig-regions.", control_tag, ".tsv")))

cat(paste0("Saved outputs to: ", outdir, "\n"))
