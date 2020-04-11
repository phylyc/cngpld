#' Read a copy-number segmentation file
#'
#' @param file  path to the seg file containing log ratios
#' @return a \code{data.frame}
#' @export
read_seg <- function(file) {

	# open a file connection
	if (is.character(file)) {
		f <- base::file(file, "rt");
	} else {
		f <- file;
	}

	header <- scan(f, character(), nlines=1, comment.char="#", quiet=TRUE);

	# continue reading the file connection
	# assume first four columns are:
	# name (character), chromosome (character), start (integer), end (integer)
	x <- read.table(f, sep="\t", header=FALSE,
		colClasses=c("character", "character", "numeric", "numeric",
			rep("numeric", length(header)-4)), ...);
	colnames(x) <- header;

	if (is.character(file)) {
		close(f);
	}

	colnames(x) <- c("sample", "chromosome", "start", "end", "nprobes", "logr");

	x
}

#' Convert data.frame of segments to GRanges
#' 
#' @param seg     \code{data.frame} containing log ratios
#' @param genome  genome build (e.g. hg19)
#' @export
seg_to_gr <- function(seg, genome=NA) {
	chroms <- as.character(seg$chromosome);
	if (!grepl("^chr", chroms[1])) {
		chroms <- paste0("chr", chroms);
	}

	gr <- GRanges(
		seqnames = paste0("chr", seg$chroms),
		ranges = IRanges(start = seg$start, end = seg$end),
		logr = seg$logr,
		sample = seg$sample,
		nprobes = seg$nprobes
	);
	genome(gr) <- genome;

	gr
}

# Compute mean copy-number at a genomic position
# @param gr  GRanges object
mean_cn_at_position <- function(gr, pos, direction, cutoff) {
	idx <- from(findOverlaps(ranges(gr), IRanges(start=pos, end=pos)));
	logr <- direction * gr$logr[idx];
	idx2 <- logr > cutoff;
	sum(logr[idx2]) / sum(idx2)
}

#' Summarize copy-number values
#' 
#' @param gr        GRanges object
#' @param direction  direction of change: 1 for amp, -1 for del
#' @param cutoff     absolute threshold for copy-number log ratio
#' @return  a \code{cn_summary} object
#' @export
summarize_cn <- function(gr, direction=1, cutoff=0.1) {
	positions <- sort(unique(c(start(gr), end(gr))));
	y <- unlist(lapply(positions,
		function(pos) {
			mean_cn_at_position(gr, pos, direction=direction, cutoff=cutoff)
		}
	));
	structure(
		data.frame(
			position = positions,
			value = y
		),
		class = "cn_summary"
	)
}

#' Collapse repeated runs
#'
#' @param  d          \code{cn_summary} object
#' @param  digits     number of digits to keep after rounding
#' @param  max.len    max length that can be collapsed into one data point;
#'                    all other runs will be collapsed to the two end points
#' @return \code{cn_summary} object
#' @export
collapse_runs <- function(d, digits=2, max.len=2e6) {
	# get repeated runs
	r <- rle(round(d$value, digits=digits));

	# start and end index of repeated runs
	starts <- cumsum(c(1, r$lengths[-length(r$lengths)]));
	ends <- cumsum(r$lengths);

	# genomic size of these runs
	gsizes <- d$position[ends] - d$position[starts] + 1;

	ridx <- gsizes <= max.len;

	# repeated runs that can be collapsed into a single data point (mid point)
	positions <- region_center(d$position[starts[ridx]], d$position[ends[ridx]]);
	values <- mapply(
		function(s, e) {
			mean(d$value[s:e])
		},
		starts[ridx], ends[ridx]
	);

	# repeated runs that need to collapsed into two data points (end points)
	large.starts <- starts[!ridx];
	large.ends <- ends[!ridx];
	positions.1 <- d$position[large.starts];
	positions.2 <- d$position[large.ends];
	values.1 <- d$value[large.starts];
	values.2 <- d$value[large.ends];

	dc <- data.frame(
		position = c(positions, positions.1, positions.2),
		value = c(values, values.1, values.2)
	);

	dc[order(dc$position), ]
}

#' Compare copy-number segments using GPLDIFF.
#'
#' Assume that input seg file or \code{data.frame} contain copy-number
#' data as log ratios. The column names must be:
#' \code{sample, chromosome, start, end, nprobes, logr}.
#'
#' @param case     seg file name for the case cohort,
#'                 or \code{data.frame}
#' @param control  seg file name for the control cohort,
#'                 or \code{data.frame}
#' @param param    initial parameter values to \code{gpldiff()}
#' @param hparams  hyperparameter values to \code{gpldiff()}
#' @param collapse whether to collapse runs of repeats (improves speed)
#' @param cutoff   absolute threshold for copy-number log ratio
#' @param verbose  verbosity level; none: 0, info: 1, debug: 2
#' @param ...      other parameters to \code{gpldiff()}
#' @return a list of \code{gpldiff} objects
compare_seg <- function(case, control, params=NULL, hparams=NULL, collapse=TRUE, cutoff=0.1, verbose=1, ...) {
	
	if (is.character(case)) {
		case <- read_seg(case);
	}

	if (is.character(control)) {
		control <- read_seg(control);
	}

	if (is.null(hparams)) {
		hparams <- default_hparams();
	}

	# split by chromosome
	case.split <- split(case, list(chromosome = case$chromosome));
	control.split <- split(control, list(chromosome = control$chromosome));

	# seg contains only segments from one chromosome
	summarize_amp_del <- function(seg) {
		gr <- seg_to_gr(seg);
		d.amp <- summarize_cn(gr, direction=1, cutoff=cutoff);
		d.del <- summarize_cn(gr, direction=-1, cutoff=cutoff);
		if (collapse) {
			d.amp <- collapse_runs(d.amp);
			d.del <- collapse_runs(d.del);
		}
		list(
			amp = d.amp,
			del = d.del
		)
	}

	s.case <- mclapply(case.split, summarize_amp_del);
	s.control <- mclapply(control.split, summarize_amp_del);

	amp.case <- lapply(s.case, function(x) x$amp);
	amp.control <- lapply(s.control, function(x) x$amp);

	del.case <- lapply(s.case, function(x) x$del);
	del.control <- lapply(s.control, function(x) x$del);

	dsets <- list(
		amp = mcmapply(prepare_cn, amp.case, amp.control, SIMPLIFY=FALSE);
		del =  mcmapply(prepare_cn, del.case, del.control, SIMPLIFY=FALSE);
	);

	# fit models for amp and del separately, and for each chromosome
	msets <- mapply(
		function(dset, type) {
			mclapply(names(dset),
				function(i) {
					if (verbose >= 1) {
						message("Processing ", type, i)
					}
					gpldiff(dset[[i]], params=params, hparams=hparams, verbose=verbose, ...)
			});
		},
		dsets, names(dests),
		SIMPLIFY=FALSE
	);

	# zip matching data and model together
	fits <- mapply(
		function(dset, mset) {
			mapply(function(d, m) list(data = d, model = m), dset, mset, SIMPLIFY=FALSE);
		},
		dsets, msets,
		SIMPLIFY=FALSE
	);
	
	structure(fits, class="cn_gpldiffs")
}


#' Summarize \code{cn_gpldiffs} object.
#'
#' @param object  a \code{cn_gpldiffs} object
#' @return a \code{list} of \code{data.frame}
summary.cn_gpldiffs <- function(object, direction=1) {
	# cn_gpldiffs is organized as a list of lists,
	# with first level grouping by type (amp vs. del)

	rsets <- lapply(
		object,
		function(fset) {
			mclapply(fset,
				function(fit) {
					find_sig_regions(fit$model, fit$data, direction=direction, process=FALSE);
				}
			)
		}
	);

	regions <- mcmapply(
		function(regions, type) {
			d <- process_regions(combine_regions(rsets$amp), direction=direction);
			if (nrow(d) > 0) {
				data.frame(type = type, d)
			} else {
				NULL
			}
		},
		rsets, c("Amp", "Del"),
		SIMPLIFY=FALSE
	);

	process_cn_regions(rbind(regions$amp, regions$del))
}

