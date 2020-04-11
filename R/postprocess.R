# process copy-number regions
process_cn_regions <- function(regions) {
	# convert position from Mbp back to bp
	regions$start <- regions$start * 1e6;
	regions$end <- regions$end * 1e6;

	filter_centromere_regions(regions)
}

# filter out regions that overlap with padded centromere regions
filter_centromere_regions <- function(regions, padding=10e6, genome="hg19") {
	cen_chroms <- centromeres[[genome]]$chromosome;
	cen_starts <- centromeres[[genome]]$start - padding;
	cen_ends <- centromeres[[genome]]$end + padding;

	idx <- match(regions$chromosome, cen_chroms);
	regions[!overlap(regions$start, regions$end, cen_starts[idx], cen_ends[idx]), ]
}

# combine regions from different chromosomes together
combine_regions <- function(regions) {
	combined <- do.call(rbind,
		mapply(
			function(d, chrom) {
				if (!is.null(d)) {
					data.frame(
						chromosome = as.integer(chrom),
						d
					)
				} else {
					NULL
				}
			},
			regions,
			names(regions),
			SIMPLIFY = FALSE
		)
	);
	rownames(combined) <- NULL;

	combined
}

