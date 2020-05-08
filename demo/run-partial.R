library(io)
library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)

library(devtools)
load_all();


# Read data #################################################################

genome <- "hg19";

seg.luad <- qread("tcga-luad.seg");
seg.lusc <- qread("tcga-lusc.seg");

chrno <- "14";
chrom <- paste0("chr", chrno);

seg.luad.chr <- seg.luad[seg.luad$chromosome %in% chrno, ];
seg.lusc.chr <- seg.lusc[seg.lusc$chromosome %in% chrno, ];


# Run analysis ###############################################################

# complete analysis took 10 min
# on a single-thread of a Core i7 CPU @ 2.93GHz

#options(mc.cores=1);
fits <- compare_segs(seg.luad.chr, seg.lusc.chr);


# Examine results ############################################################

# significant regions in LUAD
regions.luad <- summary(fits);
filter(regions.luad, end - start + 1 > 2e6, abs(ldiff) > 0.1, fdr < 0.05, n_obs > 10)

qdraw(
	{
		with(fits$amp[["14"]],  # NKX2-1 (TFF-1) amplicon
			plot(model, data, which=c("response", "latent", "odds"), xlab="position (Mbp)")
		)
	},
	width = 5, height = 10,
	file = "cngpld_luad_nkx2-1.pdf"
)

