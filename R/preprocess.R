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
prepare_cn <- function(case, control, unit=1e6) {
	# concatenate control data with case data
	data <- list(
		# total sample size
		J = nrow(control) + nrow(case),
		# convert position from bp to Mbp
		x = c(control$position, case$position) / unit,
		# group membership
		g = c(rep(-0.5, nrow(control)), rep(0.5, nrow(case))),
		y = c(control$value, case$value)
	);

	# data are currently sorted by cohort but need to be sorted by position
	idx <- order(data$x);
	data$x <- data$x[idx];
	data$g <- data$g[idx];
	data$y <- data$y[idx];

	data
}

