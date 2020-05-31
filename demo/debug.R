library(devtools)
load_all();


# Read data #################################################################

genome <- "hg19";

seg.luad <- read_seg("tcga-luad.seg");
seg.lusc <- read_seg("tcga-lusc.seg");

chrom <- "3";

seg.luad.chr <- seg.luad[seg.luad$chromosome == chrom, ];
seg.lusc.chr <- seg.lusc[seg.lusc$chromosome == chrom, ];


# Run analysis ###############################################################

gr.luad <- seg_to_gr(seg.luad.chr);
d.amp.luad <- summarize_cn(gr.luad, direction=1, cutoff=0.1);
d.del.luad <- summarize_cn(gr.luad, direction=-1, cutoff=0.1);
summary(d.del.luad$value)

s.amp.luad <- d.amp.luad;
s.amp.luad$value <- as.numeric(smooth(s.amp.luad$value));
s.del.luad <- d.del.luad;
s.del.luad$value <- as.numeric(smooth(s.del.luad$value));

with(d.amp.luad, plot(position, value, type="l"))
with(s.amp.luad, lines(position, value, col="royalblue"))

with(d.del.luad, plot(position, value, type="l"))
with(s.del.luad, lines(position, value, col="royalblue"))

c.amp.luad <- collapse_runs(s.amp.luad);
c.del.luad <- collapse_runs(s.del.luad);

with(c.amp.luad, plot(position, value, type="l"))
with(c.del.luad, plot(position, value, type="l"))

print(str(c.amp.luad))
print(str(c.del.luad))


gr.lusc <- seg_to_gr(seg.lusc.chr);
d.del.lusc <- summarize_cn(gr.lusc, direction=-1, cutoff=0.1);
summary(d.del.lusc$value)
s.del.lusc <- d.del.lusc
s.del.lusc$value <- as.numeric(smooth(s.del.lusc$value));
c.del.lusc <- collapse_runs(s.del.lusc);

with(c.del.lusc, plot(position, value, type="l"))

str(c.del.luad)
str(c.del.lusc)

del.set <- prepare_cn(c.del.luad, c.del.lusc);

del.fit <- gpldiff(del.set);

plot(del.fit, del.set);

