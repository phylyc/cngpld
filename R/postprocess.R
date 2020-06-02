# process copy-number regions
process_cn_regions <- function(regions) {
	if (is.null(regions)) {
		NULL
	} else {
		# convert position from Mbp back to bp
		regions$start <- regions$start * 1e6;
		regions$end <- regions$end * 1e6;

		if (is.null(regions)) {
			NULL
		} else {
			# sort by false discovery rate
			regions[order(regions$fdr), ]
		}
	}
}

# combine regions from different chromosomes together
combine_regions <- function(regions) {
	combined <- do.call(rbind,
		mapply(
			function(d, chrom) {
				if (!is.null(d)) {
					data.frame(
						chromosome = chrom,
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

