library(io)
library(dplyr)

library(devtools)
load_all();

#options(mc.cores=1);

case.fn <- "gistic2/tcga-luad/scores.gistic";
control.fn <- "gistic2/tcga-lusc/scores.gistic";

out.fn <- filename("cngpld-gistic2", date=NA);
fits.fn <- insert(out.fn, tag="luad-vs-lusc", ext="rds");

#chroms <- c("11", "14");
chroms <- NULL;

case <- read_gistic(case.fn);
control <- read_gistic(control.fn);

if (!is.null(chroms)) {
	case <- filter(case, chromosome %in% chroms);
	control <- filter(control, chromosome %in% chroms);
}


fits <- compare_gistics(case, control);
qwrite(fits, fits.fn);


regions.case <- summary(fits);
regions.case.f <- filter(regions.case, end - start + 1 > 2e6, abs(ldiff) > 0.1, fdr < 0.05, n_obs > 10);
print(regions.case.f)

qdraw(
	{
		with(fits$amp[["14"]],  # NKX2-1 (TFF-1) amplicon
			plot(model, data, which=c("response", "latent", "odds"), xlab="position (Mbp)")
		)
	},
	width = 5, height = 10,
	file = insert(out.fn, tag=c("luad", "nkx2-1"), ext="pdf")
)


regions.control <- summary(fits, direction=-1);
regions.control.f <- filter(regions.control, end - start + 1 > 2e6, abs(ldiff) > 0.1, fdr < 0.05, n_obs > 10);
print(regions.control.f)

qdraw(
	{
		with(fits$amp[["11"]],  # CCND1 amplicon
			plot(model, data, which=c("response", "latent", "odds"), xlab="position (Mbp)")
		)
	},
	width = 5, height = 10,
	file = insert(out.fn, tag=c("lusc", "ccnd1"), ext="pdf")
)

qwrite(regions.case, insert(out.fn, tag=c("sig-regions", "luad"), ext="tsv"));
qwrite(regions.control, insert(out.fn, tag=c("sig-regions", "lusc"), ext="tsv"));

