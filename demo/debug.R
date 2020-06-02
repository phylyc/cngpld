library(devtools)
load_all();


# Read data #################################################################

genome <- "hg19";

seg.luad <- read_seg("tcga-luad.seg");
seg.lusc <- read_seg("tcga-lusc.seg");

seg.luad$chromosome[seg.luad$chromosome == "23"] <- "X";
seg.lusc$chromosome[seg.lusc$chromosome == "23"] <- "X";

#chrom <- "3";
#chrom <- "5";

chrom <- "14";
#chrom <- "11";

seg.luad.chr <- seg.luad[seg.luad$chromosome == chrom, ];
seg.lusc.chr <- seg.lusc[seg.lusc$chromosome == chrom, ];

cutoff <- 0.5;
res <- 100;


# Run analysis ###############################################################


seg.luad.chr <- split_chromosome_arm_seg(seg.luad.chr, genome);
seg.lusc.chr <- split_chromosome_arm_seg(seg.lusc.chr, genome);


#gr.luad <- seg_to_gr(seg.luad.chr);
#gr.luad <- seg_to_gr(median_center_seg(seg.luad.chr));
gr.luad <- seg_to_gr(wmean_center_seg(seg.luad.chr));
#gr.luad <- seg_to_gr(wmean_center_arm_seg(seg.luad.chr, genome));
#gr.luad <- seg_to_gr(wmean_center_seg(split_chromosome_arm_seg(seg.luad.chr, genome)));

d.amp.luad <- summarize_cn(gr.luad, direction=1, cutoff=cutoff);
d.del.luad <- summarize_cn(gr.luad, direction=-1, cutoff=cutoff);
summary(d.amp.luad$value)
summary(d.del.luad$value)

gr <- gr.luad;

direction <- 1;

#positions <- sort(unique(c(start(gr), end(gr))));
#vlist <- lapply(positions,
#	function(pos) {
#		ov <- findOverlaps(ranges(gr), IRanges(start=pos, end=pos));
#		idx <- as.matrix(ov)[,1];
#		logr <- direction * gr$logr[idx];
#		ifelse(logr > cutoff, logr, 0)
#	}
#);
#values <- unlist(vlist);
#values <- values[values > 0];

#library(MASS)
#alpha <- coef(fitdistr(values[values < 0.1], "exponential"));

#hist(values, breaks=1000, freq=FALSE, xlim=c(-1, 1))
#h <- hist(values, breaks=1000, freq=FALSE);
#curve(dexp(x, rate=alpha), add=TRUE)

#plot(log(h$mids), log(h$density))
#idx <- h$density > 0 & h$mids > cutoff;
#f <- lm(log(h$density[idx]) ~ log(h$mids[idx]))
#abline(a=coef(f)[1], b=coef(f)[2])

#plot(h$mids, log(h$density))
#abline(a=log(alpha), b=-alpha)
#abline(h=log(alpha), )

#qqnorm(values, pch='.')


s.amp.luad <- d.amp.luad;
s.amp.luad$value <- as.numeric(smooth(smooth(s.amp.luad$value)));
s.del.luad <- d.del.luad;
s.del.luad$value <- as.numeric(smooth(smooth(s.del.luad$value)));

c.amp.luad <- collapse_runs(s.amp.luad, res);
c.del.luad <- collapse_runs(s.del.luad, res);

with(d.amp.luad, plot(position, value, type="l"))
with(s.amp.luad, lines(position, value, col="royalblue"))
with(c.amp.luad, lines(position, value, col="firebrick"))

with(d.del.luad, plot(position, value, type="l"))
with(s.del.luad, lines(position, value, col="royalblue"))
with(c.del.luad, lines(position, value, col="firebrick"))

print(str(c.amp.luad))
print(str(c.del.luad))


#gr.lusc <- seg_to_gr(seg.lusc.chr);
#gr.lusc <- seg_to_gr(median_center_seg(seg.lusc.chr));
gr.lusc <- seg_to_gr(wmean_center_seg(seg.lusc.chr));
#gr.lusc <- seg_to_gr(wmean_center_arm_seg(seg.lusc.chr, genome));
#gr.lusc <- seg_to_gr(wmean_center_seg(split_chromosome_arm_seg(seg.lusc.chr, genome)));

d.amp.lusc <- summarize_cn(gr.lusc, direction=1, cutoff=cutoff);
summary(d.amp.lusc$value)

d.del.lusc <- summarize_cn(gr.lusc, direction=-1, cutoff=cutoff);
summary(d.del.lusc$value)

s.amp.lusc <- d.amp.lusc
s.amp.lusc$value <- as.numeric(smooth(smooth(s.amp.lusc$value)));
c.amp.lusc <- collapse_runs(s.amp.lusc, res);

s.del.lusc <- d.del.lusc
s.del.lusc$value <- as.numeric(smooth(smooth(s.del.lusc$value)));
c.del.lusc <- collapse_runs(s.del.lusc, res);

with(d.amp.lusc, plot(position, value, type="l"))
with(s.amp.lusc, lines(position, value, col="royalblue"))
with(c.amp.lusc, lines(position, value, col="firebrick"))

with(d.del.lusc, plot(position, value, type="l"))
with(s.del.lusc, lines(position, value, col="royalblue"))
with(c.del.lusc, lines(position, value, col="firebrick"))

print(str(c.amp.lusc))
print(str(c.del.lusc))


amp.set <- prepare_cn(c.amp.luad, c.amp.lusc);
amp.fit <- gpldiff(amp.set);

plot(amp.fit, amp.set);

del.set <- prepare_cn(c.del.luad, c.del.lusc);
del.fit <- gpldiff(del.set);

plot(del.fit, del.set);

