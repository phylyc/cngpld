library(io)
library(GenomicRanges)

genome <- "hg19";
seg <- qread("tcga-pancan.seg");
colnames(seg) <- c("sample", "chromosome", "start", "end", "nprobes", "logr");

tss <- qread("tss.tsv");

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

qwrite(seg.luad, "tcga-luad.seg");
qwrite(seg.lusc, "tcga-lusc.seg");

