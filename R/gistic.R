#' Read a GISTIC summary file
#'
#' @param file  path to the file
#' @return a \code{data.frame}
#' @export
read_gistic <- function(file) {
	x <- read.table(file, sep="\t", header=TRUE);
	# rename colnames; in particular, G score is renamed to value
	colnames(x) <- c("type", "chromosome", "start", "end", "nlq", "value", "mean_amplitude", "frequency");
	# calculate position as the center of each segment
	x$position <- with(x, region_center(start, end));

	x
}

#' Compare GISTIC scores using GPLDIFF.
#'
#' @import parallel gpldiff
#'
#' @param case     file name of GISTIC scores table for case cohort,
#'                 or \code{data.frame}
#' @param control  file name of GISTIC scores table for control cohort,
#'                 or \code{data.frame}
#' @param param    initial parameter values to \code{gpldiff()}
#' @param hparams  hyperparameter values to \code{gpldiff()}
#' @param verbose  verbosity level; none: 0, info: 1, debug: 2
#' @param ...      other paramsters to \code{gpldiff()}
#' @return a \code{cn_gpldiff} object
#' @export
compare_gistics <- function(case, control, params=NULL, hparams=NULL, verbose=1, ...) {

	if (is.character(case)) {
		case <- read_gistic(case);
	}

	if (is.character(control)) {
		control <- read_gistic(control);
	}

	if (is.null(hparams)) {
		hparams <- default_hparams();
	}

	case.split <- split(case, list(type = case$type, chromosome = case$chromosome));
	control.split <- split(control, list(type = control$type, chromosome = control$chromosome));

	# prepare data and fit model
	datas <- mcmapply(prepare_cn, case.split, control.split, SIMPLIFY=FALSE);
	models <- mclapply(names(datas), function(i) { 
		if (verbose >= 1) {
			message("Processing ", i)
		}
		gpldiff::gpldiff(datas[[i]], params=params, hparams=hparams, ...)
	});

	# organize results by amp and del
	idx.amp <- grep("Amp", names(datas));
	idx.del <- grep("Del", names(datas));

	dsets <- list(
		amp = datas[idx.amp],
		del = datas[idx.del]
	);
	names(dsets$amp) <- sub("Amp.", "", names(dsets$amp));
	names(dsets$del) <- sub("Del.", "", names(dsets$del));

	msets <- list(
		amp = models[idx.amp],
		del =	models[idx.del] 
	);
	names(msets$amp) <- sub("Amp.", "", names(msets$amp));
	names(msets$del) <- sub("Del.", "", names(msets$del));

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

read_gistic_peaks <- function(fname) {
	x <- read.table(fname, sep="\t", header=TRUE, check.names=FALSE);	
	y <- data.frame(
		type = ifelse(grepl("Amplification", x[["Unique Name"]]), "Amp",
			ifelse(grepl("Deletion", x[["Unique Name"]]), "Del", NA)),
		string_to_coordinate(x[["Wide Peak Limits"]]),
		q = x[["q values"]]
	);
	y
}

summary_append_gistic_peaks <- function(x, peaks, jaccard.cut = 0) {
	y <- x;
	y$gistic_jaccard <- NA;
	y$gistic_q <- NA;
	for (i in 1:nrow(x)) {
		r <- x[i, ];
		s <- peaks[peaks$type == r$type & peaks$chromosome == r$chromosome, ];
		if (nrow(s) > 0) {
			jaccard <- jaccard_similarity(r$start, r$end, s$start, s$end);
			idx <- which.max(jaccard);
			if (jaccard[idx] > jaccard.cut) {
				y$gistic_jaccard[i] <- jaccard[idx];
				y$gistic_q[i] <- s$q[idx];
			}
		}
	}

	y
}

