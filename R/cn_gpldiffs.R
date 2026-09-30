#' Summarize \code{cn_gpldiffs} object.
#'
#' @import parallel
#'
#' @param object     a \code{cn_gpldiffs} object
#' @param direction  comparison direction:
#'                    1 indicates case vs. control,
#'                   -1 indicates control vs. case.
#' @param genome     genome build to use for filtering centromeres;
#'                   \code{NULL} suppresses filtering
#' @param lodds.cut  log posterior odds threshold defining candidate support regions
#' @param lodds.core optional higher log posterior odds threshold defining peak cores;
#'                   \code{NULL} preserves the original single-threshold behavior
#' @param peak.prominence minimum log-odds drop required to split adjacent cores
#' @param ...        other parameters passed to underlying functions
#' @return a \code{list} of \code{data.frame}
#' @method summary cn_gpldiffs
#' @export
summary.cn_gpldiffs <- function(object, direction=1, genome=NULL, lodds.cut=10,
	lodds.core=NULL, peak.prominence=1.5, ...) {
	# cn_gpldiffs is organized as a list of lists,
	# with first level grouping by type (amp vs. del)

	rsets <- lapply(
		object,
		function(fset) {
			mclapply(fset,
				function(fit) {
					gpldiff::find_sig_regions(
						fit$model, fit$data, direction=direction, process=FALSE,
						lodds.cut=lodds.cut, lodds.core=lodds.core,
						min.prominence=peak.prominence, ...
					)
				}
			)
		}
	)

	regions <- mcmapply(
		function(rset, type) {
			d <- gpldiff::process_regions(combine_regions(rset), direction=direction)
			if (!is.null(d) && nrow(d) > 0) {
				data.frame(type = type, d)
			} else {
				NULL
			}
		},
		rsets, c("Amp", "Del"),
		SIMPLIFY=FALSE
	)

	rs <- process_cn_regions(rbind(regions$amp, regions$del))

	if (is.null(genome)) {
		rs
	} else {
		filter_centromere_regions(rs, genome=genome, ...)
	}
}

