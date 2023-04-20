library(devtools)
load_all()

seg <- read_seg("tcga-pancan.seg");
tss <- read.table("tss.tsv", sep="\t", header=TRUE);

# rename chrX
seg$chromosome[seg$chromosome == "23"] <- "X";

code.kirc <- as.character(tss$code[tss$study_name == "Kidney renal clear cell carcinoma"]);
seg.kirc <- do.call(rbind, lapply(code.kirc,
	function(code) {
		seg[grep(paste0("TCGA-", code), seg$sample), ]
	}
));
write_seg(seg.kirc, "tcga-kirc.seg");

# collect LUAD samples
code.luad <- as.character(tss$code[tss$study_name == "Lung adenocarcinoma"]);
seg.luad <- do.call(rbind, lapply(code.luad,
	function(code) {
		seg[grep(paste0("TCGA-", code), seg$sample), ]
	}
));

# collect LUSC samples
code.lusc <- as.character(tss$code[tss$study_name == "Lung squamous cell carcinoma"]);
seg.lusc <- do.call(rbind, lapply(code.lusc,
	function(code) {
		seg[grep(paste0("TCGA-", code), seg$sample), ]
	}
));

samples.luad <- unique(seg.luad$sample);
samples.lusc <- unique(seg.lusc$sample);

pheno <- data.frame(
	sample = c(samples.luad, samples.lusc),
	cancer_subtype = c(
		rep("LUAD", length(samples.luad)),
		rep("LUSC", length(samples.lusc))
	)
);

write_seg(seg.luad, "tcga-luad.seg");
write_seg(seg.lusc, "tcga-lusc.seg");
write.table(pheno, "pheno.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE);

