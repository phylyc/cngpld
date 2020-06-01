library(io)

library(devtools)
load_all();


# Read data #################################################################

genome <- "hg19";

seg.luad <- read_seg("tcga-luad.seg");
seg.lusc <- read_seg("tcga-lusc.seg");

# NKX2-1 (TFF-1) amplicon is on chr14
chrno <- c("11", "14");

seg.luad.chr <- seg.luad[seg.luad$chromosome %in% chrno, ];
seg.lusc.chr <- seg.lusc[seg.lusc$chromosome %in% chrno, ];


# Run analysis ###############################################################

# complete analysis of chr14 took 15 min
# on a single-thread of a Core i7 CPU @ 2.93GHz

#options(mc.cores=1);
fits <- compare_segs(seg.luad.chr, seg.lusc.chr);


# Examine results ############################################################

# significant regions in LUAD
regions.luad <- summary(fits, genome=genome);
if (!is.null(regions.luad)) {
	idx <- with(regions.luad, end - start + 1 > 2e6 & abs(ldiff) > 0.10 & fdr < 0.05 & n_obs > 10);
	print(regions.luad[idx, ]);
}

# significant regions in LUSC
regions.lusc <- summary(fits, direction=-1, genome=genome);
if (!is.null(regions.lusc)) {
	idx <- with(regions.lusc, end - start + 1 > 2e6 & abs(ldiff) > 0.10 & fdr < 0.05 & n_obs > 10);
	print(regions.lusc[idx, ]);
}

for (ch in chrno) {
	qdraw(
		{
			with(fits$amp[[ch]],  
				plot(model, data, which=c("response", "latent", "odds"), xlab="position (Mbp)")
			)
		},
		width = 5, height = 10,
		file = sprintf("cngpld_luad_amp_chr%s.pdf", ch)
	)
}

