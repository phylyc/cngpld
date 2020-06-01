library(io)
library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)

library(devtools)
load_all();


# Read data #################################################################

genome <- "hg19";
fits.fn <- "cngpld_luad-vs-lusc.rds";


# Run analysis ###############################################################

if (file.exists(fits.fn)) {

	warning("Warning: Using cached model")
	fits <- qread(fits.fn);

} else {

	seg.luad <- read_seg("tcga-luad.seg");
	seg.lusc <- read_seg("tcga-lusc.seg");

	# complete analysis took 49595 s = 13.8 h
	# on a single-thread of a Core i7 CPU @ 2.93GHz

	#options(mc.cores=1);
	fits <- compare_segs(seg.luad, seg.lusc);

	qwrite(fits, fits.fn);
}


# Examine results ############################################################

# problematic deletion profiles:
print(which(unlist(lapply(fits$del, function(x) ! is(x$model, "gpldiff")))))
# chroms 1 11 13 15 17 19 20 22  3  5  7  9

# remove problematic deletion profiles
fits.orig <- fits;
idx <- unlist(lapply(fits$del, function(x) is(x$model, "gpldiff")));
fits$del <- fits$del[idx];

# significant regions in LUAD
regions.luad <- summary(fits, genome=genome);
print(filter(regions.luad, end - start + 1 > 2e6, abs(ldiff) > 0.10, fdr < 0.05, n_obs > 10))

qdraw(
	{
		with(fits$amp[["14"]],  # NKX2-1 (TFF-1) amplicon
			plot(model, data, which=c("response", "latent", "odds"), xlab="position (Mbp)")
		)
	},
	width = 5, height = 10,
	file = "cngpld_luad_nkx2-1.pdf"
)

# significant regions in LUSC
regions.lusc <- summary(fits, direction=-1, genome=genome);
print(filter(regions.lusc, end - start + 1 > 2e6, abs(ldiff) > 0.10, fdr < 0.05, n_obs > 10))

qdraw(
	{
		with(fits$amp[["11"]],  # CCND1 amplicon
			plot(model, data, which=c("response", "latent", "odds"), xlab="position (Mbp)")
		)
	},
	width = 5, height = 10,
	file = "cngpld_lusc_ccnd1.pdf"
)

with(fits$amp[["6"]], plot(model, data))   # chr6p arm
with(fits$amp[["3"]], plot(model, data))   # chr3q arm
with(fits$amp[["2"]], plot(model, data))   # chr2q amplicon
with(fits$amp[["9"]], plot(model, data))   # chr9p amplicon
with(fits$del[["2"]], plot(model, data))   # LRP1B deletion
#with(fits$del[["3"]], plot(model, data))   # chr3q deletion
#with(fits$del[["5"]], plot(model, data))   # chr3p deletion

qwrite(regions.luad, "cngpld_sig-regions_luad.tsv");
qwrite(regions.lusc, "cngpld_sig-regions_lusc.tsv");


# plot results summary

revlog_trans <- function(base = exp(1)) {
	trans <- function(x) -log(x, base)
	inv <- function(x) base^(-x)
	trans_new(paste0("reverselog-", format(base)), trans, inv, 
						log_breaks(base = base), 
						domain = c(1e-100, Inf))
}

regions.all <- rbind(
	data.frame(regions.luad, group="case"),
	data.frame(regions.lusc, group="control")
);
regions.all$group <- factor(regions.all$group, levels=c("control", "case", "NS"));
idx <- with(regions.all, end - start + 1 > 2e6 & abs(ldiff) > 0.15 & fdr < 0.05 & n_obs > 10);
regions.all$keep <- 0.75;
regions.all$keep[idx] <- 1.0;
regions.all$group[!idx] <- "NS";

regions.all$gene <- NA;
regions.all$gene[regions.all$start_idx == 568] <- "CCND1";
regions.all$gene[regions.all$start_idx == 112] <- "NKX2-1";

qdraw(
	ggplot(regions.all, aes(x=ldiff, y=fdr, alpha=keep, colour=group, label=gene)) + theme_classic() +
		geom_point(show.legend=FALSE) +
		geom_text_repel(show.legend=FALSE, nudge_y=0.1, nudge_x=0.095) +
		geom_vline(xintercept=0) +
		geom_vline(xintercept=c(0.15, -0.15), linetype=3, colour="grey60") +
		geom_hline(yintercept=0.05, linetype=3, colour="grey60") +
		scale_y_continuous(trans=revlog_trans(10), sec.axis = dup_axis(name=NULL)) +
		scale_colour_manual(values=c("#0073C2FF", "#EFC000FF", "#333333FF")) +
		xlim(-0.6, 0.6) +
		xlab("latent difference") + ylab("false discovery rate") +
		annotate("text", label=levels(regions.all$group)[1:2], x=c(-0.5, 0.5), y=0.2)
	,
	width = 6, height = 3, 
	file = "cngpld_luad-vs-lusc_volanco.pdf"
);

