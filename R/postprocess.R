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

#' Filter out regions that overlap with padded centromere regions
#'
#' @param regions  \code{data.frame} of significant regions
#' @param padding  padding factor relative to the size of the centromere
#' @param genome   genome build
#' @return a \code{data.frame} of filtered regions
#' @export
filter_centromere_regions <- function(regions, padding=1, genome=c("hg19", "hg38")) {
	genome <- match.arg(genome);

	# centromere coordinates are stored in 1-based
	cens <- centromeres[[genome]];
	cen_chroms <- cens$chromosome;
	cen_sizes <- cens$end - cens$start;
	cen_starts <- (cens$start + 1) - padding * cen_sizes;
	cen_ends <- cens$end + padding * cen_sizes;

	idx <- match(regions$chromosome, cen_chroms);
	regions[!gpldiff:::overlap(regions$start, regions$end, cen_starts[idx], cen_ends[idx]), , drop=FALSE]
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

