library(devtools)
load_all();


# Read data #################################################################

genome <- "hg19";

seg.luad <- read_seg("tcga-luad.seg");
seg.lusc <- read_seg("tcga-lusc.seg");

#chrom <- "3";
chrom <- "14";

seg.luad.chr <- seg.luad[seg.luad$chromosome == chrom, ];
seg.lusc.chr <- seg.lusc[seg.lusc$chromosome == chrom, ];


# Run analysis ###############################################################

gr.luad <- seg_to_gr(median_center_seg(seg.luad.chr));
d.amp.luad <- summarize_cn(gr.luad, direction=1, cutoff=0.1);
d.del.luad <- summarize_cn(gr.luad, direction=-1, cutoff=0.1);
summary(d.amp.luad$value)
summary(d.del.luad$value)

gr <- gr.luad;

direction <- 1;
cutoff <- 0.1;

positions <- sort(unique(c(start(gr), end(gr))));
vlist <- lapply(positions,
	function(pos) {
		ov <- findOverlaps(ranges(gr), IRanges(start=pos, end=pos));
		idx <- as.matrix(ov)[,1];
		logr <- direction * gr$logr[idx];
		ifelse(logr > cutoff, logr, 0)
	}
);
values <- unlist(vlist);
#values <- values[values > 0];

library(MASS)
alpha <- coef(fitdistr(values, "exponential"));

hist(values, breaks=100, freq=FALSE)
curve(dexp(x, rate=alpha), add=TRUE)


s.amp.luad <- d.amp.luad;
s.amp.luad$value <- as.numeric(smooth(smooth(s.amp.luad$value)));
s.del.luad <- d.del.luad;
s.del.luad$value <- as.numeric(smooth(smooth(s.del.luad$value)));

with(d.amp.luad, plot(position, value, type="l"))
with(s.amp.luad, lines(position, value, col="royalblue"))

with(d.del.luad, plot(position, value, type="l"))
with(s.del.luad, lines(position, value, col="royalblue"))

c.amp.luad <- collapse_runs(s.amp.luad);
c.del.luad <- collapse_runs(s.del.luad);

with(d.amp.luad, plot(position, value, type="l"))
with(c.amp.luad, lines(position, value, col="royalblue"))

with(d.del.luad, plot(position, value, type="l"))
with(c.del.luad, lines(position, value, col="royalblue"))

print(str(c.amp.luad))
print(str(c.del.luad))


gr.lusc <- seg_to_gr(median_center_seg(seg.lusc.chr));

d.amp.lusc <- summarize_cn(gr.lusc, cutoff=0.1);
summary(d.amp.lusc$value)

d.del.lusc <- summarize_cn(gr.lusc, direction=-1, cutoff=0.1);
summary(d.del.lusc$value)

s.amp.lusc <- d.amp.lusc
s.amp.lusc$value <- as.numeric(smooth(smooth(s.amp.lusc$value)));
c.amp.lusc <- collapse_runs(s.amp.lusc);

s.del.lusc <- d.del.lusc
s.del.lusc$value <- as.numeric(smooth(smooth(s.del.lusc$value)));
c.del.lusc <- collapse_runs(s.del.lusc);

with(d.amp.lusc, plot(position, value, type="l"))
with(c.amp.lusc, lines(position, value, col="royalblue"))

with(d.del.lusc, plot(position, value, type="l"))
with(c.del.lusc, lines(position, value, col="royalblue"))

print(str(c.amp.lusc))
print(str(c.del.lusc))


amp.set <- prepare_cn(c.amp.luad, c.amp.lusc);
amp.fit <- gpldiff(amp.set);

plot(amp.fit, amp.set);

del.set <- prepare_cn(c.del.luad, c.del.lusc);
del.fit <- gpldiff(del.set);

plot(del.fit, del.set);

