genome ?= hg19
table ?= cytoBand

url = http://hgdownload.soe.ucsc.edu/goldenPath/$(genome)/database

all: $(genome)/$(table).txt $(genome)/centromeres.tsv
	

$(genome)/$(table).txt:
	mkdir -p $(genome)
	curl $(url)/$(table).txt.gz | gunzip > $(genome)/$(table).txt

$(genome)/centromeres.tsv: $(genome)/cytoBand.txt
	Rscript cytoband_to_centromeres.r $< $@
