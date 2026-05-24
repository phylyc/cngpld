#' Read a copy-number segmentation file
#'
#' @param file  path to the seg file containing log ratios
#' @return a \code{data.frame}
#' @export
read_seg <- function(file, ...) {

	# open a file connection
	if (is.character(file)) {
		f <- base::file(file, "rt")
	} else {
		f <- file
	}

	header <- scan(f, character(), sep="\t", nlines=1, comment.char="#", quiet=TRUE)

	# continue reading the file connection
	# assume first four columns are:
	# name (character), chromosome (character), start (integer), end (integer)
	x <- read.table(f, sep="\t", header=FALSE,
		colClasses=c("character", "character", "numeric", "numeric",
			rep("numeric", length(header)-4)), ...)
	colnames(x) <- header

	if (is.character(file)) {
		close(f)
	}

	colnames(x) <- c("sample", "chromosome", "start", "end", "nprobes", "logr")

	x
}

#' Write a copy-number segmentation file
#'
#' @param x     \code{data.frame} with segmentation data
#' @param file  path to the seg file containing log ratios
#' @return a \code{data.frame}
#' @export
write_seg <- function(x, file, ...) {
	write.table(x, file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t", ...)
}

# Median center each chromosome for each sample
median_center_seg <- function(seg) {
	segs <- split(seg, list(seg$sample, seg$chromosome))
	d <- do.call(rbind,
		lapply(segs,
			function(s) {
				s$logr <- s$logr - median(s$logr)
				s
			}
		)
	)
	rownames(d) <- NULL

	d
}

# center by substract overall mean of chromosome
wmean_center_seg <- function(seg) {
	segs <- split(seg, list(seg$sample, seg$chromosome))
	d <- do.call(rbind,
		lapply(segs,
			function(s) {
				w <- s$end - s$start + 1
				w <- w / sum(w)
				m <- sum(w * s$logr)
				# do not allow subtraction to induce a copy-number change
				s$logr <- s$logr - ifelse(abs(s$logr) > abs(m), m, s$logr)
				s
			}
		)
	)
	rownames(d) <- NULL

	d
}

# center by subtracting overall mean of chromosome arm
wmean_center_arm_seg <- function(seg, genome) {
	seg <- mark_chromosome_arm_seg(seg, genome)

	segs <- split(seg, list(seg$sample, seg$chromosome))
	d <- do.call(rbind,
		lapply(segs,
			function(s) {
				w <- s$end - s$start + 1
				x <- w * s$logr
				m.w <- sum(x) / sum(w)
				p.arm <- s$p_arm
				q.arm <- s$q_arm
				m.p <- sum(x[p.arm]) / sum(w[p.arm])
				m.q <- sum(x[q.arm]) / sum(w[q.arm])
				# use whole chromosome mean if segment spans both arms
				# use p arm if segment only spans p arm
				# use q arm if segment only spans q arm
				# use 0 if segment spans neither arm
				m <- ifelse(s$p_arm,
					ifelse(s$q_arm, m.w, m.p),
					ifelse(s$q_arm, m.q, 0)
				)
				s$logr <- s$logr - m
				s
			}
		)
	)
	rownames(d) <- NULL

	d
}

#' Convert data.frame of segments to GRanges
#' 
#' @import GenomicRanges IRanges GenomeInfoDb
#'
#' @param seg     \code{data.frame} containing log ratios
#' @param genome  genome build
#' @export
seg_to_gr <- function(seg, genome=NULL) {
	chroms <- as.character(seg$chromosome)
	#if (!grepl("^chr", chroms[1])) {
	#	chroms <- paste0("chr", chroms)
	#}

	gr <- GRanges(
		seqnames = chroms,
		ranges = IRanges(start = seg$start, end = seg$end),
		logr = seg$logr,
		sample = seg$sample,
		nprobes = seg$nprobes
	)
	if (!is.null(genome)) {
		genome(gr) <- genome
	}

	gr
}

# Compute mean copy-number at a genomic position
# @param gr  GRanges object
summarize_cn_at_position <- function(gr, pos, direction, cutoff) {
	# S4Vectors::from is not available in older versions (v0.8.11),
	# so we avoid using it
	ov <- findOverlaps(ranges(gr), IRanges(start=pos, end=pos))
	idx <- as.matrix(ov)[,1]
	logr <- direction * gr$logr[idx]
	idx2 <- logr > cutoff

	if (sum(idx2) > 0) {
		sum(exp(logr[idx2])) / length(logr)
	} else {
		0
	}
}

#' Summarize copy-number values
#' 
#' @param gr         GRanges object
#' @param direction  direction of change: 1 for amp, -1 for del
#' @param cutoff     absolute threshold for copy-number log ratio
#' @return  a \code{cn_summary} object
#' @export
summarize_cn <- function(gr, direction, cutoff, positions=NULL) {
	if (is.null(positions)) {
		positions <- sort(unique(c(start(gr), end(gr))))
	}
	values <- unlist(lapply(positions,
		function(pos) {
			summarize_cn_at_position(gr, pos, direction=direction, cutoff=cutoff)
		}
	))
	structure(
		data.frame(
			position = positions,
			value = values
		),
		class = "cn_summary"
	)
}

# Count copy-number alterations at a genomic position
count_cn_at_position <- function(gr, pos, direction, cutoff) {
	ov <- findOverlaps(ranges(gr), IRanges(start=pos, end=pos))
	idx <- as.matrix(ov)[,1]
	logr <- direction * gr$logr[idx]
	idx2 <- logr > cutoff

	sum(idx2)
}

# Count copy-number alterations
count_cn <- function(gr, direction, cutoff, positions=NULL) {
	if (is.null(positions)) {
		positions <- sort(unique(c(start(gr), end(gr))))
	}
	values <- unlist(lapply(positions,
		function(pos) {
			count_cn_at_position(gr, pos, direction=direction, cutoff=cutoff)
		}
	))
	structure(
		data.frame(
			position = positions,
			value = values
		),
		class = "cn_counts"
	)
}

#' Collapse repeated runs
#'
#' @param  d          \code{cn_summary} object
#' @param  res        resolution (higher order of magnitude leads to 
#'                    less coarsening)
#' @param  max.len    max length that can be collapsed into one data point;
#'                    all other runs will be collapsed to the two end points
#' @return \code{cn_summary} object
#' @export
collapse_runs <- function(d, res, max.len=2e6) {
	if (res <= 0) return(d)

	# get repeated runs
	r <- rle(round(d$value * res) / res)

	# start and end index of repeated runs
	starts <- cumsum(c(1, r$lengths[-length(r$lengths)]))
	ends <- cumsum(r$lengths)

	# genomic size of these runs
	gsizes <- d$position[ends] - d$position[starts] + 1

	ridx <- gsizes <= max.len

	if (sum(ridx) > 0) {
		# repeated runs that can be collapsed into a single data point (mid point)
		positions <- region_center(d$position[starts[ridx]], d$position[ends[ridx]])
		values <- mapply(
			function(s, e) {
				mean(d$value[s:e])
			},
			starts[ridx], ends[ridx]
		)
	} else {
		# no runs should be collapsed into a single data point
		positions <- NULL
		values <- NULL
	}

	# repeated runs that need to be collapsed into two data points (end points)
	large.starts <- starts[!ridx]
	large.ends <- ends[!ridx]
	positions.1 <- d$position[large.starts]
	positions.2 <- d$position[large.ends]
	values.1 <- d$value[large.starts]
	values.2 <- d$value[large.ends]

	dc <- data.frame(
		position = c(positions, positions.1, positions.2),
		value = c(values, values.1, values.2)
	)

	dc[order(dc$position), ]
}

# @param pair  whether to pair case and control seg for each chromosome
split_segs <- function(case, control, genome=NULL, pair=FALSE) {
	if (is.character(case)) {
		case <- read_seg(case)
	}

	if (is.character(control)) {
		control <- read_seg(control)
	}

	if (!is.null(genome)) {
		case <- split_chromosome_arm_seg(case, genome)
		control <- split_chromosome_arm_seg(control, genome)
	}

	# split by chromosome
	case.split <- split(case, list(chromosome = case$chromosome))
	control.split <- split(control, list(chromosome = control$chromosome))

	# look for common chromosomes
	chroms.common <- intersect(names(case.split), names(control.split))
	names(chroms.common) <- chroms.common
	if (length(chroms.common) == 0) {
		message("case chromosomes: ")
		cat(case.split, stderr())
		message("")
		message("control chromosomes: ")
		cat(case.split, stderr())
		message("")
		stop("Error: case and control samples have no chromosomes in common.")
	}

	if (!all(union(names(case.split), names(control.split)) %in% chroms.common)) {
		# some chromosomes are missing in case or control
		missing.in.case <- setdiff(names(control.split), names(case.split))
		missing.in.control <- setdiff(names(case.split), names(control.split))
		warning("Warning: cases and controls contain different chromosomes.\n",
			"Missing in case: ", paste(missing.in.case, collapse=", "), "\n",
			"Missing in control: ", paste(missing.in.control, collapse=", ")
		)
	}

	if (pair) {
		lapply(chroms.common,
			function(chrom) {
				list(
					case = case.split[[chrom]],
					control = control.split[[chrom]]
				)
			}
		)
	} else {
		# use common chromosomes and ensure that they are in same order
		list(
			case = case.split[chroms.common],
			control = control.split[chroms.common]
		)
	}
}

# collapse runs based on count difference in order to improve speed
collapse_runs_paired <- function(d, res=1, ...) {
	e <- data.frame(
		position = d$position,
		value = d$case - d$control
	)
	e <- collapse_runs(e, res=res, max.len=-1)

	d[match(e$position, d$position), ]
}

# Count copy-number alterations at breakpoints in seg data.frame objects
count_segs <- function(case, control,
	cn.cut=0.5, genome=NULL, verbose=1, ...
) {
	segs.split <- split_segs(case, control, genome=genome, pair=TRUE)

	# segs contain a list of case and control seg data.frames for one chromosome
	count_amp_del <- function(segs) {
		grs <- lapply(segs, function(seg) {
			seg_to_gr(wmean_center_seg(seg), genome)
		})

		# identify common positions across case and control
		positions <- sort(unique(c(start(grs$case), start(grs$control), end(grs$case), end(grs$control))))

		lapply(grs, function(gr) {
			list(
				amp = count_cn(gr, direction=1, cutoff=cn.cut, positions=positions),
				del = count_cn(gr, direction=-1, cutoff=cn.cut, positions=positions)
			)
		})
	}

	# summary s is organized by chromosomes, group, cna type
	s <- mclapply(segs.split, count_amp_del)
	
	# organize by cna type and chromosome
	list(
		amp = lapply(s, function(ss) {
			d <- data.frame(
				position = ss$case$amp$position,
				case = ss$case$amp$value,
				control = ss$control$amp$value
			)
			collapse_runs_paired(d)
		}),
		del = lapply(s, function(ss) {
			d <- data.frame(
				position = ss$case$del$position,
				case = ss$case$del$value,
				control = ss$control$del$value
			)
			collapse_runs_paired(d)
		})
	)
}

#' Compare copy-number segments using GPLDIFF.
#'
#' Assume that input seg file or \code{data.frame} contain copy-number
#' data as log ratios. The column names must be:
#' \code{sample, chromosome, start, end, nprobes, logr}.
#'
#' If \code{genome} build is provided,
#' results will be organized by chromosome arms,
#' broad copy-number events will be more effectively removed,
#' and computational speed will improve.
#' See \code{data(centromeres)} for supported genomes.
#' The chromosome nomenclature must match the target genome.
#'
#' @import parallel gpldiff
#'
#' @param case     seg file name for the case cohort,
#'                 or \code{data.frame}
#' @param control  seg file name for the control cohort,
#'                 or \code{data.frame}
#' @param param    initial parameter values to \code{gpldiff()}
#' @param hparams  hyperparameter values to \code{gpldiff()}
#' @param cn.cut   absolute threshold for copy-number log ratio
#' @param pair     analyze data in pair mode
#' @param smooth   whether to median smooth the copy-number data
#' @param cn.res   copy-number resolution; if value is positive, collapse 
#'                 runs of repeats (which improves speed) while preserving 
#'                 the target resolution
#' @param genome   genome build
#' @param verbose  verbosity level; none: 0, info: 1, debug: 2
#' @param ...      other parameters to \code{gpldiff()}
#' @return a list of \code{gpldiff} objects
#' @export
compare_segs <- function(case, control, params=NULL, hparams=NULL,
	cn.cut=0.5, pair=FALSE, smooth=TRUE, cn.res=100, weight.N.ref=100, genome=NULL, verbose=1, weight_by_cohort_size=FALSE, ...
) {

	n.case <- length(unique(case$sample))
	n.control <- length(unique(control$sample))

	if (is.null(hparams)) {
		hparams <- default_hparams()
	}

	add_n_samples <- function(d, n) {
		attr(d, "n_samples") <- n
		d
	}

	if (pair) {

		segs.split <- split_segs(case, control, genome=genome, pair=TRUE)

		# segs contain a list of case and control seg data.frames for one chromosome
		summarize_amp_del <- function(segs) {
			grs <- lapply(segs, function(seg) {
				seg_to_gr(wmean_center_seg(seg), genome)
			})

			# identify common positions across case and control
			positions <- sort(unique(c(start(grs$case), start(grs$control), end(grs$case), end(grs$control))))

			lapply(grs, function(gr) {
				d.amp <- summarize_cn(gr, direction=1, cutoff=cn.cut, positions=positions)
				d.del <- summarize_cn(gr, direction=-1, cutoff=cn.cut, positions=positions)
				if (smooth) {
					d.amp$value <- as.numeric(smooth(d.amp$value))
					d.del$value <- as.numeric(smooth(d.del$value))
				}
				list(
					amp = d.amp,
					del = d.del
				)
			})
		}

		# summary s is organized by chromosomes, group, cna type
		s <- mclapply(segs.split, summarize_amp_del)

		# organize by cna type and chromosome and collapse runs
		s.paired <- list(
			amp = lapply(s, function(ss) {
				d <- data.frame(
					position = ss$case$amp$position,
					case = ss$case$amp$value,
					control = ss$control$amp$value
				)
				collapse_runs_paired(d, res=cn.res)
			}),
			del = lapply(s, function(ss) {
				d <- data.frame(
					position = ss$case$del$position,
					case = ss$case$del$value,
					control = ss$control$del$value
				)
				collapse_runs_paired(d, res=cn.res)
			})
		)

		# organized by type, then by chromosome
		s.case <- lapply(s.paired, function(ss) {
			lapply(ss, function(d) {
				data.frame(
					position = d$position,
					value = d$case
				)
			})
		})
		s.control <- lapply(s.paired, function(ss) {
			lapply(ss, function(d) {
				data.frame(
					position = d$position,
					value = d$control
				)
			})
		})

		s.case$amp <- lapply(s.case$amp, add_n_samples, n = n.case)
		s.case$del <- lapply(s.case$del, add_n_samples, n = n.case)

		s.control$amp <- lapply(s.control$amp, add_n_samples, n = n.control)
		s.control$del <- lapply(s.control$del, add_n_samples, n = n.control)

		dsets <- list(
			amp = mcmapply(
				function(ca, co) {
					prepare_cn(ca, co, weight.N.ref = weight.N.ref, weight_by_cohort_size = weight_by_cohort_size)
				},
				s.case$amp,
				s.control$amp,
				SIMPLIFY = FALSE
			),
			del = mcmapply(
				function(ca, co) {
					prepare_cn(ca, co, weight.N.ref = weight.N.ref, weight_by_cohort_size = weight_by_cohort_size)
				},
				s.case$del,
				s.control$del,
				SIMPLIFY = FALSE
			)
		)

	} else {

		seg.split <- split_segs(case, control, genome=genome)

		# seg contains only segments from one chromosome
		summarize_amp_del <- function(seg) {
			gr <- seg_to_gr(wmean_center_seg(seg), genome)
			d.amp <- summarize_cn(gr, direction=1, cutoff=cn.cut)
			d.del <- summarize_cn(gr, direction=-1, cutoff=cn.cut)
			if (smooth) {
				d.amp$value <- as.numeric(smooth(d.amp$value))
				d.del$value <- as.numeric(smooth(d.del$value))
			}
			d.amp <- collapse_runs(d.amp, res=cn.res)
			d.del <- collapse_runs(d.del, res=cn.res)
			list(
				amp = d.amp,
				del = d.del
			)
		}

		# organized by chromosome, then by type
		s.case <- mclapply(seg.split$case, summarize_amp_del)
		s.control <- mclapply(seg.split$control, summarize_amp_del)

		amp.case <- lapply(s.case, function(x) x$amp)
		amp.control <- lapply(s.control, function(x) x$amp)

		del.case <- lapply(s.case, function(x) x$del)
		del.control <- lapply(s.control, function(x) x$del)

		amp.case <- lapply(amp.case, add_n_samples, n = n.case)
		del.case <- lapply(del.case, add_n_samples, n = n.case)

		amp.control <- lapply(amp.control, add_n_samples, n = n.control)
		del.control <- lapply(del.control, add_n_samples, n = n.control)

		dsets <- list(
			amp = mcmapply(
				function(ca, co) {
					prepare_cn(ca, co, weight.N.ref = weight.N.ref, weight_by_cohort_size = weight_by_cohort_size)
				},
				amp.case,
				amp.control,
				SIMPLIFY = FALSE
			),
			del = mcmapply(
				function(ca, co) {
					prepare_cn(ca, co, weight.N.ref = weight.N.ref, weight_by_cohort_size = weight_by_cohort_size)
				},
				del.case,
				del.control,
				SIMPLIFY = FALSE
			)
		)

	}

	# fit models for amp and del separately, and for each chromosome
	msets <- mapply(
		function(dset, type) {
			lapply(names(dset), function(i) {
					if (verbose >= 1) {
						message("Processing ", type, " ", i)
					}
				tryCatch(
					gpldiff::gpldiff(dset[[i]], params = params, hparams = hparams, verbose = verbose, ...),
					error = function(e) {
						message("FAILED: ", type, " ", i)
						message(conditionMessage(e))
						stop(e)
					}
				)
			})
		},
		dsets, names(dsets),
		SIMPLIFY = FALSE
	)

	# zip matching data and model together
	fits <- mapply(
		function(dset, mset) {
			mapply(function(d, m) list(data = d, model = m), dset, mset, SIMPLIFY=FALSE)
		},
		dsets, msets,
		SIMPLIFY=FALSE
	)
	attr(fits, "genome") <- genome

	structure(fits, class="cn_gpldiffs")
}


