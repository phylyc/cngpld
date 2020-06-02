#!/usr/bin/env Rscript

library(argparser)

pr <- arg_parser("Convert cytoband to centromeres data table");
pr <- add_argument(pr, "input", help="UCSC cytoBand input file");
pr <- add_argument(pr, "output", help="output file name");

argv <- parse_args(pr);
input.fn <- argv$input;
output.fn <- argv$output;

d <- read.table(input.fn, header=FALSE, sep="\t", stringsAsFactors=FALSE);
names(d) <- c("chromosome", "start", "end", "cytoband", "stain");

# select centromeres
cen <- d[d$stain == "acen", ];

# remove chromosome prefix
cen$chromosome <- sub("chr", "", cen$chromosome);

# merge centromeric regions on the same chromosome
cens <- split(cen, cen$chromosome);

chroms <- names(cens);
starts <- unlist(lapply(cens, function(s) min(s$start)));
ends <- unlist(lapply(cens, function(s) max(s$end)));

y <- data.frame(
	chromosome = chroms,
	# ensure that coordinates are 1-based
	start = starts + 1,
	end = ends
);

write.table(y, output.fn, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t");

