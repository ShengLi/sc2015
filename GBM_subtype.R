# https://tcga-data.nci.nih.gov/docs/publications/gbm_exp/TCGA_unified_CORE_ClaNC840.txt
source('convert.R')
tcga840=read.table("doc/TCGA_unified_CORE_ClaNC840.txt",header=T, sep="\t", stringsAsFactors = F)
tcga.genes=tcga840[-1,]
tcga.genes1=tcga.genes[tcga.genes[,2]!="" ,1:2]
tcga.genes1$mouse=foreach(g=tcga.genes1[,1], .combine=c) %do% simpleCap(g)
tcga.genes2=tcga.genes1[tcga.genes1$mouse %in% rownames(mydata),]
subtypes=unique(tcga.genes2[,2])

# subytpe index
subtype.dat=foreach(i=subtypes, .combine=rbind) %do%{
  subtypes.g=tcga.genes2[tcga.genes2[,2]==i,"mouse"]
  colMeans(mydata[subtypes.g,])-colMeans(mydata)
}
rownames(subtype.dat)=subtypes
