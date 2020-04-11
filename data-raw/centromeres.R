centromeres <- list();
centromeres$hg19 <- read.table("centromeres_hg19.seg", sep="\t", header=TRUE);

save(list = c("centromeres"), file = "../R/sysdata.rda");
