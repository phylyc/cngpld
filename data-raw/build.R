#!/usr/bin/env Rscript

# Build package data

genomes <- list.dirs(full.names=FALSE, recursive=FALSE);
names(genomes) <- genomes;

# read centromeres
centromeres <- lapply(genomes,
	function(genome) {
		read.table(sprintf("%s/centromeres.tsv", genome), sep="\t", header=TRUE)
	}
);

save(list = c("centromeres"), file = "../R/sysdata.rda");
