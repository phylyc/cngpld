library(io)

library(devtools)
load_all();

#options(mc.cores=1);

control.fname <- "gistic2/tcga-luad/scores.gistic";
case.fname <- "gistic2/tcga-lusc/scores.gistic";
out.fname <- filename("bm-luad-bmet_tcga-luad", tag="gpldiff");

case <- read_gistic(case.fname);
control <- read_gistic(control.fname);

fits <- compare_gistics(case.fname, control.fname);
qwrite(fits, insert(out.fname, ext="rds"));

results <- summary(fits);
results.f <- results[results$diff > 0.5 & results$fdr < 0.01, ];
print(results.f)

results.down <- summary(fits, direction=-1);
results.down.f <- results.down[results.down$diff < -0.5 & results.down$fdr < 0.01, ];
print(results.down.f)

qwrite(results, insert(out.fname, ext="tsv"));
qwrite(results.down, insert(out.fname, tag="down", ext="tsv"));

