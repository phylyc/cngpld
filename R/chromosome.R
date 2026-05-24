# Get centromere regions with padding
# @param padding  centromere padding factor relative to its size
get_padded_centromere_regions <- function(genome=names(centromeres), padding=1) {
	genome <- match.arg(genome)

	# centromere coordinates are stored in 1-based
	cens <- centromeres[[genome]]
	cen_sizes <- cens$end - cens$start + 1

	data.frame(
		chromosome = cens$chromosome,
		start = cens$start - padding * cen_sizes,
		end = cens$end + padding * cen_sizes
	)
}

#' Split chromosomes into chromosome arms.
#'
#' Assign segments to chromosome arms, splitting them as necessary.
#'
#' @param seg      a \code{data.frame} of segments
#' @param genome   genome build
#' @param padding  centromere padding factor relative to its size
#' @return a \code{data.frame} of modified segments
#' @export
split_chromosome_arm_seg <- function(seg, genome, padding=1) {
	if (is.null(seg) || nrow(seg) == 0) return(seg)

	cens <- get_padded_centromere_regions(genome, padding)

	# remove unknown chromosomes
	idx <- which(! seg$chromosome %in% as.character(cens$chromosome))
	if (length(idx) > 0) {
		warning("Unknown chromosomes are removed: ",
			paste(unique(seg$chromosome[idx]) ,sep=", "))
		seg <- seg[-idx, ]
	}

	idx <- match(seg$chromosome, cens$chromosome)

	p.arm <- seg$start < cens$start[idx]
	q.arm <- seg$end > cens$end[idx]

	p.arm.only <- p.arm & !q.arm
	q.arm.only <- q.arm & !p.arm
	
	# segments that only span p arms
	seg.p <- seg[p.arm & !q.arm, ]
	if (nrow(seg.p) > 0) {
		# right-truncate segments at the centromere
		seg.p$end <- pmin(seg.p$end, cens$start[match(seg.p$chromosome, cens$chromosome)] - 1)
		seg.p$chromosome <- paste0(seg.p$chromosome, "p")
	}

	# segments that only span q arms
	seg.q <- seg[q.arm & !p.arm, ]
	if (nrow(seg.q) > 0) {
		# left-truncate segments at the centromere
		seg.q$start <- pmax(seg.q$start, cens$end[match(seg.q$chromosome, cens$chromosome)] + 1)
		seg.q$chromosome <- paste0(seg.q$chromosome, "q")
	}

	# segments that span both p and q arms
	both.arm <- p.arm & q.arm

	seg.bp <- seg[both.arm, ]
	if (nrow(seg.bp) > 0) {
		# right-truncate segments at the centromere
		seg.bp$end <- pmin(seg.bp$end, cens$start[match(seg.bp$chromosome, cens$chromosome)] - 1)
		seg.bp$chromosome <- paste0(seg.bp$chromosome, "p")
	}

	seg.bq <- seg[both.arm, ]
	if (nrow(seg.bq) > 0) {
		# left-truncate segments at the centromere
		seg.bq$start <- pmax(seg.bq$start, cens$end[match(seg.bq$chromosome, cens$chromosome)] + 1)
		seg.bq$chromosome <- paste0(seg.bq$chromosome, "q")
	}

	seg2 <- rbind(seg.p, seg.bp, seg.q, seg.bq)
	seg2 <- with(seg2, seg2[order(sample, chromosome, start), ])

	seg2	
}

# Mark whether each segment spans each arm of the chromosome
mark_chromosome_arm_seg <- function(seg, genome, padding=1) {
	cens <- get_padded_centromere_regions(genome, padding)

	idx <- match(seg$chromosome, cens$chromosome)
	seg$p_arm <- seg$start < cens$start[idx]
	seg$q_arm <- seg$end > cens$end[idx]

	seg
}

#' Filter out regions that overlap with padded centromere regions
#'
#' @param regions  \code{data.frame} of significant regions
#' @param padding  centromere padding factor relative to its size
#' @param genome   genome build
#' @return a \code{data.frame} of filtered regions
#' @export
filter_centromere_regions <- function(regions, genome, padding=1) {
	cens <- get_padded_centromere_regions(genome, padding)

	idx <- match(sub("p|q", "", regions$chromosome), cens$chromosome)
	regions[!gpldiff:::overlap(regions$start, regions$end, cens$start[idx], cens$end[idx]), , drop=FALSE]
}

