
all: tcga-pancan.seg tcga-luad.seg tcga-lusc.seg
	

tcga-pancan.seg:
	curl -L -o tcga-pancan.seg http://api.gdc.cancer.gov/data/00a32f7a-c85f-4f86-850d-be53973cbc4d

tcga-luad.seg: tcga-pancan.seg
	Rscript prepare.R

tcga-lusc.seg: tcga-pancan.seg
	Rscript prepare.R

clean:
	rm -f *.rds *.pdf *.seg

