# process copy-number regions
process_cn_regions <- function(regions) {
	if (is.null(regions)) {
		NULL
	} else {
		# convert position from Mbp back to bp
		regions$start <- regions$start * 1e6;
		regions$end <- regions$end * 1e6;

		regions <- filter_centromere_regions(regions);

		if (is.null(regions)) {
			NULL
		} else {
			# sort by false discovery rate
			regions[order(regions$fdr), ]
		}
	}
}

# filter out regions that overlap with padded centromere regions
# @param padding  padding factor relative to the size of the centromere
filter_centromere_regions <- function(regions, padding=1, genome="hg19") {
	# centromere coordinates are stored in 1-based
	cens <- centromeres[[genome]];
	cen_chroms <- cens$chromosome;
	cen_sizes <- cens$end - cens$start;
	cen_starts <- (cens$start + 1) - padding * cen_sizes;
	cen_ends <- cens$end + padding * cen_sizes;

	idx <- match(regions$chromosome, cen_chroms);
	regions[!gpldiff:::overlap(regions$start, regions$end, cen_starts[idx], cen_ends[idx]), ]
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

