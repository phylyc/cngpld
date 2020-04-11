coordinate_to_string <- function(x) {
	sprintf("chr%s:%d-%d",
		x$chromosome,
		x$start,
		x$end
	)
}

# NB this is human specific
chromname_to_integer <- function(s, species="hsa") {
	s <- gsub("chr", "", as.character(s), fixed=TRUE);
	if (species == "hsa") {
		ifelse(s == "X", 23,
			ifelse(s == "Y", 24,
				ifelse(s == "MT", 25,
					suppressWarnings(as.integer(s))
				)
			)
		)
	} else {
		stop("Species ", species, " is not supported.");
	}
}

integer_to_chromname <- function(x, species="hsa") {
	if (species == "hsa") {
		ifelse(x == 23, "X",
			ifelse(x == 24, "Y",
				ifelse(s == 25, "MT",
					suppressWarnings(as.character(s))
				)
			)
		)
	} else {
		stop("Species ", species, " is not supported.");
	}
}

string_to_coordinate <- function(s) {
	# removed parenthesized tokens
	s <- gsub("\\(.*\\)", "", s);
	# split on : or -
	# illegal expression can be parsed here too
	tokenss <- strsplit(s, ":|-");
	do.call(rbind, lapply(
		tokenss,
		function(tokens) {
			data.frame(
				chromosome = chromname_to_integer(tokens[1]),
				start = as.integer(tokens[2]),
				end = as.integer(tokens[3])
			)
		}
	))
}
