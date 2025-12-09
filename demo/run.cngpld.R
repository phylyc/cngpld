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
  library(cngpld)
  library(rjson)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- "."
study <- "CASE_vs_CONTROL"
use_cache <- TRUE
if (length(args) > 0) {
  for (arg in args) {
    if (grepl("^--dir=", arg)) { outdir <- sub("^--dir=", "", arg) }
    if (grepl("^--study=", arg)) { study <- sub("^--study=", "", arg) }
    if (grepl("^--no-cache", arg)) { use_cache <- FALSE }
  }
}
message("Output directory:", outdir)

dir.create(paste0(outdir, "/cngpld"), showWarnings = FALSE, recursive = TRUE)


# Configurations ##############################################################

case <- unlist(strsplit(study, "_vs_"))[1]
control <- unlist(strsplit(study, "_vs_"))[2]

genome <- "hg19"
fits.fn <- paste0(outdir, "/cngpld/", case, "-vs-", control, ".rds")

# lodds.cut - probability of discovery
# 3 <-> 0.953
# 4 <-> 0.982
# 5 <-> 0.993
lodds.cut <- 3
min_tCR <- 0.3  # absolute threshold for copy-number log ratio to be considered in summary statistics
min_nprobes <- 4  # minimum number of targets supporting a segment

# For interval annotations:
fdr_threshold <- 0.1  # defines where score == 0.9 for fc >> 1
fc_threshold <- 1.15  # defines where score == 0.9 for fdr << 0
frac_patients_threshold <- 0.01
score_threshold <- 0.5  # annotation threshold
min_seg_size <- 1e5
n_obs_threshold <- 5


# Run analysis ################################################################

seg.case <- cngpld::read_seg(
  paste0(outdir, "/tcga-", case, ".seg.gz")
) %>% mutate( logr = pmax(log(2) * logr, -3) ) %>% filter( nprobes >= min_nprobes )
seg.control <- cngpld::read_seg(
  paste0(outdir, "/tcga-", control, ".seg.gz")
) %>% mutate( logr = pmax(log(2) * logr, -3) ) %>% filter( nprobes >= min_nprobes )

if (file.exists(fits.fn) & use_cache) {

  # cat("Using cached model.")
  fits <- io::qread(fits.fn)

} else {

  options(mc.cores = 8)
  fits <- cngpld::compare_segs(
    seg.case, seg.control,
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

regions.case <- summary(fits, genome = genome, lodds.cut = lodds.cut)
regions.control <- summary(fits, direction = -1, genome = genome, lodds.cut = lodds.cut)
if (is.null(regions.control)) {
  regions.control <- regions.case[FALSE,]
}

plot_region <- function(chr_arm, profile, profile_str = "amp") {
  qdraw({
    with(profile[[chr_arm]], plot(model, data, which = c("response", "latent", "odds"), xlab = "position (Mbp)"))
  }, width = 5, height = 10, file = paste0(outdir, "/cngpld/chr/chr", chr_arm, "_", profile_str, ".pdf"))
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

  regions[, region_id := NULL]               # drop temp ID
  regions[]
}

regions.case    <- annotate_frac_patients(regions.case,    seg.case,    min_tCR)
regions.control <- annotate_frac_patients(regions.control, seg.control, min_tCR)


# The "score" is a measure of confidence in this interval as a characteristic of
# the difference between those two cohorts. It scales with 1 - fdr and with
# 1 - fc^k for fc < 1 and some power k. This score could be used e.g. for downstream
# pathway or over-representation analysis.

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

io::qwrite(regions.case, paste0(outdir, "/cngpld/cngpld_sig-regions_", case, ".tsv"))
io::qwrite(regions.control, paste0(outdir, "/cngpld/cngpld_sig-regions_", control, ".tsv"))
