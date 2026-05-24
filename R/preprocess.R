region_center <- function(start, end) {
	start + floor((end - start)/2)
}

#' Prepare copy-number data for GPLDIFF
#'
#' @param case     \code{cn_summary} for the case cohort
#' @param control  \code{cn_summary} for the control cohort
#' @param unit     unit size of the genomic coorindate
#' @return a \code{gpldiff_data} object
#' @export
prepare_cn <- function(case, control, unit = 1e6, n.case = NULL, n.control = NULL, weight.N.ref = 100, weight_by_cohort_size = FALSE) {
	if (is.null(n.case)) {
		n.case <- attr(case, "n_samples")
	}
	if (is.null(n.control)) {
		n.control <- attr(control, "n_samples")
	}

	if (is.null(n.case)) {
		n.case <- weight.N.ref
	}
	if (is.null(n.control)) {
		n.control <- weight.N.ref
	}

	# Backward-compatible default: original unweighted behavior
	if (weight_by_cohort_size) {
		w.case <- n.case / weight.N.ref
		w.control <- n.control / weight.N.ref
	} else {
		w.case <- 1
		w.control <- 1
	}

	# concatenate control data with case data
	data <- list(
		# total sample size
		J = nrow(control) + nrow(case),
		# convert position from bp to Mbp
		x = c(control$position, case$position) / unit,
		# group membership
		g = c(rep(-0.5, nrow(control)), rep(0.5, nrow(case))),
		y = c(control$value, case$value),
		w = c(rep(w.control, nrow(control)), rep(w.case, nrow(case)))
	)

	# data are currently sorted by cohort but need to be sorted by position
	idx <- order(data$x)
	data$x <- data$x[idx]
	data$g <- data$g[idx]
	data$y <- data$y[idx]
	data$w <- data$w[idx]

	data
}

