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
#' @param case     file name of GISTIC scores table for case cohort
#' @param control  file name of GISTIC scores table for control cohort
#' @param param    initial parameter values to \code{gpldiff()}
#' @param hparams  hyperparameter values to \code{gpldiff()}
#' @param verbose  verbosity level; none: 0, info: 1, debug: 2
#' @param ...      other paramsters to \code{gpldiff()}
#' @return a list of \code{gpldiff} objects
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

	data.sets <- mcmapply(prepare_cn, case.split, control.split, SIMPLIFY=FALSE);
	models <- mclapply(names(data.sets), function(i) { 
		if (verbose >= 1) {
			message("Processing ", i)
		}
		gpldiff::gpldiff(data.sets[[i]], params=params, hparams=hparams, ...)
	});

	fits <- mapply(function(d, m) list(data = d, model = m), data.sets, models, SIMPLIFY=FALSE)
	
	structure(fits, class="gistic_gpldiffs")
}


#' Summarize \code{gistic_gpldiffs} object.
#'
#' @param object  a \code{gistic_gpldiffs} object
#' @return a \code{list} of \code{data.frame}
summary.gistic_gpldiffs <- function(object, direction=1) {
	# gistic_gpldiffs is organized as a flat list

	regions.all <- lapply(
		object,
		function(fit) {
			gpldiff::find_sig_regions(fit$model, fit$data, direction=direction, process=FALSE);
		}
	);

	regions.amp <- regions.all[grep("Amp", names(regions.all)), drop=FALSE];
	names(regions.amp) <- sub("Amp.", "", names(regions.amp));
	regions.amp <- gpldiff::process_regions(combine_regions(regions.amp), direction=direction);
	if (!is.null(regions.del) && nrow(regions.amp) > 0) {
		regions.amp <- data.frame(
			type = "Amp",
			regions.amp
		);
	} else {
		regions.amp <- NULL;
	}

	regions.del <- regions.all[grep("Del", names(regions.all)), drop=FALSE];
	names(regions.del) <- sub("Del.", "", names(regions.del));
	regions.del <- gpldiff::process_regions(combine_regions(regions.del), direction=direction);
	if (!is.null(regions.del) && nrow(regions.del) > 0) {
		regions.del <- data.frame(
			type = "Del",
			regions.del
		);
	} else {
		regions.del <- NULL;
	}

	process_cn_regions(
		rbind(regions.amp, regions.del)
	)
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

