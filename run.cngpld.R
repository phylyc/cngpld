suppressPackageStartupMessages({
  library(io)
  library(data.table)
  library(dplyr)
  library(ggplot2)
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
	make_option(c("--case_file"), type="character", default=NULL, help="Path to the case segmentation file (REQUIRED).", metavar="FILE"),
	make_option(c("--control_file"), type="character", default=NULL, help="Path to the control segmentation file (REQUIRED).", metavar="FILE"),
	make_option(c("--outdir"), type="character", default=NULL, help="Output directory for all results and cache (REQUIRED).", metavar="DIR"),

	## OPTIONAL
	make_option(c("--genome"), type="character", default="hg19", help="Reference genome version (default: %default).", metavar="STRING"),
	# Caching (Action is important here for a boolean flag)
	make_option(c("--use-cache"), action="store_true", default=FALSE, help="Flag to enable caching of the model (default: disabled)."),
	# Statistical & CNV Thresholds
	make_option(c("--lodds_cut"), type="numeric", default=3, help="Probability of discovery cut-off (default: %default).", metavar="NUM"),
	make_option(c("--min_tCR"), type="numeric", default=0.4, help="Absolute copy-ratio threshold for CNV summary stats (default: %default).", metavar="NUM"),
	make_option(c("--min_nprobes"), type="integer", default=4,  help="Minimum number of probes supporting a segment (default: %default).", metavar="INT"),
	# Annotation and Score Thresholds
	make_option(c("--fdr_threshold"), type="numeric", default=0.1,  help="FDR threshold for annotation scoring (default: %default).", metavar="NUM"),
	make_option(c("--fc_threshold"), type="numeric", default=1.1,  help="Fold-change threshold for annotation scoring (default: %default).", metavar="NUM"),
	make_option(c("--frac_patients_threshold"), type="numeric", default=0.01,  help="Minimum fraction of patients required to retain an interval (default: %default).", metavar="NUM"),
	make_option(c("--score_threshold"), type="numeric", default=0.5, help="Annotation confidence threshold (default: %default).", metavar="NUM"),
	make_option(c("--min_seg_size"), type="integer", default=1e5, help="Minimum segment size (base pairs) for significance (default: %default).", metavar="INT"),
	make_option(c("--n_obs_threshold"), type="integer", default=5, help="Minimum number of observations (samples) required for significance (default: %default).", metavar="INT")
)
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

# Explicitly check for required arguments (those defined with default=NULL)
required_args <- c("case_file", "control_file", "outdir")
missing_args <- required_args[sapply(opt[required_args], is.null)]
if (length(missing_args) > 0) {
  stop(paste("Missing arguments:", paste(paste0("--", missing_args), collapse=", ")))
}

print(opt) # Print all parsed options for monitoring
case_file <- opt$case_file
control_file <- opt$control_file
outdir <- opt$outdir
genome <- opt$genome
use_cache <- opt$`use-cache` # Use backticks if you keep the hyphen in the variable name
lodds_cut <- opt$lodds_cut
min_tCR <- opt$min_tCR
min_nprobes <- opt$min_nprobes
fdr_threshold <- opt$fdr_threshold
fc_threshold <- opt$fc_threshold
frac_patients_threshold <- opt$frac_patients_threshold
score_threshold <- opt$score_threshold
min_seg_size <- opt$min_seg_size
n_obs_threshold <- opt$n_obs_threshold

dir.create(paste0(outdir), showWarnings = FALSE, recursive = TRUE)

# Function to extract file label from path
get_file_label <- function(filepath) {
	if (is.na(filepath)) {
		return(NA_character_)
	}
	base_name <- basename(filepath)
	# Removes everything from the first dot to the end
	label <- sub("\\..*$", "", base_name)
	return(label)
}

case <- get_file_label(case_file)
control <- get_file_label(control_file)

# Run analysis ################################################################
seg.case <- cngpld::read_seg(
  case_file
) %>% mutate( logr = pmax(log(2) * logr, -3) ) %>% filter( nprobes >= min_nprobes )

seg.control <- cngpld::read_seg(
  control_file
) %>% mutate( logr = pmax(log(2) * logr, -3) ) %>% filter( nprobes >= min_nprobes )

fits.fn <- paste0(outdir, "/", case, "-vs-", control, ".rds")
if (file.exists(fits.fn) & use_cache) {
  # cat("Using cached model.")
  fits <- io::qread(fits.fn)

} else {

  options(mc.cores = 1)
  fits <- cngpld::compare_segs(
    seg.case,
    seg.control,
    genome = genome,
    cn.cut = min_tCR,  # absolute threshold for copy-number log ratio to be considered in summary statistics
    smooth = TRUE,  # whether to median smooth the copy-number data
    cn.res = 1e4,  # aggregated log copy-ratio signal resolution 1/cn.res
    pair = FALSE,  # whether case and control cohorts are paired samples from the same patients
    verbose = 1,
  )

  io::qwrite(fits, fits.fn)
}


# Examine results #############################################################

# remove problematic deletion profiles
fits.orig <- fits
idx <- unlist(lapply(fits$del, function(x) is(x$model, "gpldiff")))
fits$del <- fits$del[idx]

regions.case <- summary(fits, genome = genome, lodds.cut = lodds_cut)
regions.control <- summary(fits, direction = -1, genome = genome, lodds.cut = lodds_cut)
if (is.null(regions.control)) {
  regions.control <- regions.case[FALSE,]
}

# Create plots for each chr arm
plot_region <- function(chr_arm, profile, profile_str = "amp") {
  qdraw({
    with(profile[[chr_arm]], plot(model, data, which = c("response", "latent", "odds"), xlab = "position (Mbp)"))
  }, width = 5, height = 10, file = paste0(outdir, "/", chr_arm, "_", profile_str, ".pdf"))
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

annotate_frac_patients <- function(regions, seg, min_tCR) {
  # total number of samples in this cohort
  n_samples <- uniqueN(seg$sample)

  # key by intervals for fast overlap
  setkey(seg,     chromosome, start, end)
  setkey(regions, chromosome, start, end)

  # initialize result
  regions[,    region_id := .I]
  regions[, frac_patients := 0]

  ## Amplifications: logr > min_tCR
  seg_amp <- seg[logr >= min_tCR]
  if (nrow(seg_amp)) {
    hits_amp <- foverlaps(seg_amp, regions[type == "Amp"], nomatch = 0L)
    # hits_amp has both seg + region columns; count unique samples per region
    amp_counts <- hits_amp[, .(n_samples_hit = uniqueN(sample)), by = region_id]
    regions[amp_counts, frac_patients := n_samples_hit / n_samples, on = "region_id"]
  }

  ## Deletions: logr < -min_tCR  (adjust sign if your Seg.CN is coded differently)
  seg_del <- seg[logr <= -min_tCR]
  if (nrow(seg_del)) {
    hits_del <- foverlaps(seg_del, regions[type == "Del"], nomatch = 0L)
    del_counts <- hits_del[, .(n_samples_hit = uniqueN(sample)), by = region_id]
    regions[del_counts, frac_patients := n_samples_hit / n_samples, on = "region_id"]
  }

  regions[, region_id := NULL] # drop temp ID
  regions[]
}

regions.case <- annotate_frac_patients(regions.case, seg.case, min_tCR)
regions.control <- annotate_frac_patients(regions.control, seg.control, min_tCR)


# The "score" is a measure of confidence in this interval as a characteristic of the difference between those two cohorts.
    # Scales with 1 - fdr and with 1 - fc^k for fc < 1 and some power k. 
    # This score could be used e.g. for downstream pathway or over-representation analysis.
beta = -log(fdr_threshold)
x_abslog <- function(f, t = exp(1)) { return(abs(log(f) / log(t)) )}
sig <- function(x) { return( 1 - exp(-beta * x) ) }  # = 1 - exp(beta * log(f) / log(t)) = 1 - f ^ (beta / log(t))
def_score <- function(fdr, fc) { return(sig(x_abslog(fdr, t=fdr_threshold)) * sig(x_abslog(fc, t=fc_threshold)) ) }

regions.case <- regions.case %>%
  mutate(score = def_score(fdr=fdr, fc=exp(ldiff)), is_significant = "S") %>%
  mutate(IGV = paste0("chr", chromosome, ":", start, "-", end)) %>%
  arrange(desc(score))
idx <- with(
  regions.case,
  frac_patients >= frac_patients_threshold
  & score >= score_threshold
  & end - start + 1 > min_seg_size
  & n_obs >= n_obs_threshold
)
regions.case$is_significant[!idx] = "NS"

regions.control <- regions.control %>%
  mutate(score = def_score(fdr=fdr, fc=exp(ldiff)), is_significant = "S") %>%
  mutate(IGV = paste0("chr", chromosome, ":", start, "-", end)) %>%
  arrange(desc(score))
idx <- with(
  regions.control,
  frac_patients >= frac_patients_threshold
  & score >= score_threshold
  & end - start + 1 > min_seg_size
  & n_obs >= n_obs_threshold
)
regions.control$is_significant[!idx] = "NS"

io::qwrite(regions.case, file.path(outdir, paste0("cngpld_sig-regions_", case, ".tsv")))
io::qwrite(regions.control, file.path(outdir, paste0("cngpld_sig-regions_", control, ".tsv")))

cat(paste0("Saved outputs to: ", outdir, "\n"))