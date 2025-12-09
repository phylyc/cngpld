library(io)
library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)

library(devtools)
load_all()


# Configurations #############################################################

genome <- "hg19";
fits.fn <- "cngpld_luad-vs-lusc.rds";
use_cache <- FALSE

case <- "luad"
control <- "lusc"

min_tCR <- 0.5
fdr_threshold <- 0.1
ldiff_threshold <- 0.05
n_obs_threshold <- 8


# Run analysis ###############################################################

if (file.exists(fits.fn) & use_cache) {

	warning("Warning: Using cached model")
	fits <- qread(fits.fn);

} else {

	seg.case <- read_seg("tcga-luad.seg") %>% mutate( logr = log(2) * logr )
	seg.control <- read_seg("tcga-lusc.seg") %>% mutate( logr = log(2) * logr )

	# complete analysis took 49595 s = 13.8 h
	# on a single-thread of a Core i7 CPU @ 2.93GHz

	options(mc.cores=1);
  fits <- compare_segs(seg.case, seg.control, genome = genome, cn.cut = min_tCR, pair = FALSE)

	qwrite(fits, fits.fn);
}


# Examine results ############################################################

# problematic deletion profiles:
cat("Failed profiles:\n")
print(which(unlist(lapply(fits$del, function(x) !is(x$model, "gpldiff")))))

# remove problematic deletion profiles
fits.orig <- fits
idx <- unlist(lapply(fits$del, function(x) is(x$model, "gpldiff")))
fits$del <- fits$del[idx]

# significant regions in CASE
regions.case <- summary(fits, genome = genome, lodds.cut = 3)
print(filter(regions.case, end - start + 1 > 1e6, abs(ldiff) > ldiff_threshold, fdr <= fdr_threshold, n_obs >= n_obs_threshold))

# significant regions in CONTROL
regions.control <- summary(fits, direction = -1, genome = genome, lodds.cut = 3)
if (is.null(regions.control)) {
  regions.control <- regions.case[FALSE,]
}
print(filter(regions.control, end - start + 1 > 1e6, abs(ldiff) > ldiff_threshold, fdr <= fdr_threshold, n_obs >= n_obs_threshold))

plot_region <- function(chr_arm, profile, profile_str = "amp") {
  qdraw({
    with(profile[[chr_arm]], plot(model, data, which = c("response", "latent", "odds"), xlab = "position (Mbp)"))
  }, width = 5, height = 10, file = paste0("chr/chr", chr_arm, "_", profile_str, ".pdf"))
}

for (arm in (filter(regions.case, type == "Amp")$chromosome)) {
  plot_region(chr_arm = arm, profile = fits$amp, profile_str = "amp")
}
for (arm in (filter(regions.case, type == "Del")$chromosome)) {
  plot_region(chr_arm = arm, profile = fits$del, profile_str = "del")
}


# significant regions in LUAD
regions.luad <- summary(fits, genome=genome, lodds.cut=3);
print(filter(regions.luad, end - start + 5 > 1e6, abs(ldiff) > 0.10, fdr < 0.05, n_obs > 10))

qdraw(
	{
		with(fits$amp[["14q"]],  # NKX2-1 (TFF-1) amplicon
			plot(model, data, which=c("response", "latent", "odds"), xlab="position (Mbp)")
		)
	},
	width = 5, height = 10,
	file = "cngpld_luad_nkx2-1.pdf"
)

# significant regions in LUSC
regions.lusc <- summary(fits, direction=-1, genome=genome, lodds.cut=3);
print(filter(regions.lusc, end - start + 1 > 1e6, abs(ldiff) > 0.10, fdr < 0.05, n_obs > 10))

qdraw(
	{
		with(fits$amp[["11q"]],  # CCND1 amplicon
			plot(model, data, which=c("response", "latent", "odds"), xlab="position (Mbp)")
		)
	},
	width = 5, height = 10,
	file = "cngpld_lusc_ccnd1.pdf"
)

with(fits$amp[["3q"]], plot(model, data))   # chr3q amplicon containing SOX2
with(fits$amp[["7q"]], plot(model, data))   # CDK6 amplicon
with(fits$amp[["8p"]], plot(model, data))   # 8p amplicon containing NSD3 and FGFR
with(fits$amp[["9p"]], plot(model, data))   # unknown target

# observed in both LUSC and LUAD, but enriched in LUSC
with(fits$amp[["19q"]], plot(model, data))  # chr19q amplicon containing CCNE1
with(fits$del[["9p"]], plot(model, data))   # CDKN2A/B deletion


qwrite(regions.case, paste("cngpld_sig-regions_", case, ".tsv", sep=""))
qwrite(regions.control, paste("cngpld_sig-regions_", control, ".tsv", sep=""))


revlog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}

regions.all <- rbind(
  data.frame(regions.case, group = "case"),
  data.frame(regions.control, group = "control")
)
regions.all$group <- factor(regions.all$group, levels = c("control", "case", "NS"))
idx <- with(
  regions.all,
  end - start + 1 > 1e5 &
    abs(ldiff) > ldiff_threshold &
    fdr <= fdr_threshold &
    n_obs >= n_obs_threshold
)
regions.all$keep <- 0.75
regions.all$keep[idx] <- 1.0
regions.all$group[!idx] <- "NS"

regions.all$gene <- NA

# iterate over all rows in regions.case and annotate the chromosome arm:
for (i in 1:nrow(regions.all)) {
  if (regions.all$ldiff[i] > 1 & regions.all$fdr[i] < fdr_threshold / 2) {
    regions.all$gene[i] <- paste0(regions.all$chromosome[i], regions.all$arm[i])
  }
}

regions.all$gene[regions.all$chromosome == "14q" & regions.all$group == "case"] <- "NKX2-1";
regions.all$gene[regions.all$chromosome == "11q" & regions.all$group == "control"] <- "CCND1";
regions.all$gene[regions.all$chromosome == "3q" & regions.all$group == "control"] <- "SOX2";
regions.all$gene[regions.all$chromosome == "8p" & regions.all$group == "control"] <- "NSD3";
regions.all$gene[regions.all$chromosome == "7q" & regions.all$group == "control"] <- "CDK6";
regions.all$gene[regions.all$chromosome == "19q" & regions.all$group == "control"] <- "CCNE1";


ymin = min(regions.all$fdr) / 1.1
ymax = max(regions.all$fdr) / 0.9
xmin = min(regions.all$ldiff[regions.all$fdr < fdr_threshold]) - 0.5
xmax = max(regions.all$ldiff[regions.all$fdr < fdr_threshold]) + 0.5
# xmin = min(regions.all$ldiff) - 0.5
# xmax = max(regions.all$ldiff) + 0.5
# xmin = -10
# xmax = 1
xadjust = 0.1 * max(c(-xmin, xmax))
yadjust = 2 * ymin
qdraw(
  ggplot(regions.all, aes(x = ldiff, y = fdr, alpha = keep, colour = group, label = gene)) +
    theme_classic() +
    geom_point(show.legend = FALSE) +
    geom_text_repel(show.legend = FALSE, nudge_y = 0.2, nudge_x = 0.04, max.overlaps = 25) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(ldiff_threshold, -ldiff_threshold), linetype = 3, colour = "grey60") +
    geom_hline(yintercept = fdr_threshold, linetype = 3, colour = "grey60") +
    xlim(xmin, xmax) +
    scale_x_continuous(trans = log_trans(10)) +
    scale_y_continuous(trans = revlog_trans(10)) +
    scale_colour_manual(values = c("#0073C2FF", "#EFC000FF", "#333333FF")) +
    # ylim(ymax, ymin) +
    xlab("latent difference") +
    ylab("false discovery rate") +
    annotate("text", label = levels(regions.all$group)[1:2], x = c(-xadjust, xadjust), y = yadjust, hjust = 0.5, vjust = 1)
  ,
  width = 6, height = 3, file = paste0("cngpld_", case, "-vs-", control, "_volanco.log.pdf")
)
qdraw(
  ggplot(regions.all, aes(x = ldiff, y = fdr, alpha = keep, colour = group, label = gene)) +
    theme_classic() +
    geom_point(show.legend = FALSE) +
    geom_text_repel(show.legend = FALSE, nudge_y = 0.2, nudge_x = 0.04, max.overlaps = 20) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(ldiff_threshold, -ldiff_threshold), linetype = 3, colour = "grey60") +
    geom_hline(yintercept = fdr_threshold, linetype = 3, colour = "grey60") +
    xlim(xmin, xmax) +
    scale_y_continuous(trans = revlog_trans(10), sec.axis = dup_axis(name = NULL)) +
    scale_colour_manual(values = c("#0073C2FF", "#EFC000FF", "#333333FF")) +
    # ylim(ymax, ymin) +
    xlab("latent difference") +
    ylab("false discovery rate") +
    annotate("text", label = levels(regions.all$group)[1:2], x = c(-xadjust, xadjust), y = yadjust, hjust = 0.5, vjust = 1)
  ,
  width = 6, height = 3, file = paste0("cngpld_", case, "-vs-", control, "_volanco.pdf")
)