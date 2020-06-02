library(devtools)
load_all()

seg <- read_seg("tcga-pancan.seg");
tss <- read.table("tss.tsv", sep="\t", header=TRUE);

# rename chrX
seg$chromosome[seg$chromosome == "23"] <- "X";

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

write_seg(seg.luad, "tcga-luad.seg");
write_seg(seg.lusc, "tcga-lusc.seg");

