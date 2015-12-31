# use science 2014 paper for stemness in GBM: DOI: 10.1126/science.1254257
library(doMC)
source("lib/convert.R")
# gene list from Figure 3B: new_science2014_bradley_stemness.txt
sci.stem=unique(read.table("new_science2014_bradley_stemness.txt", header=F, sep="\t", stringsAsFactors = F)$V1)
# add stemness associated genes that are not in the main figure from Supplemental Figure 19
sci.stem1=unique(intersect(c("Sox4","Sox9", "Fzd3", "Ptn", "Nfib", "Kdm4A", "Rfx4",foreach(i = sci.stem1, .combine=c) %do% simpleCap(tolower(i))), rownames(mydata)))
stem.dat=data.frame(stage=splitn(colnames(mydata), "-",1),stemness=colMeans(mydata[sci.stem1,])-colMeans(mydata))
ggplot(stem.dat, aes(x=stemness, colour=stage)) + stat_ecdf()+stage.colour.scale+
  ylab("cumulative distribution function")  + gb2.theme("bottom")
ks.test(stem.dat[stem.dat$stage=="T1","stemness"],stem.dat[stem.dat$stage=="T2","stemness"])
ks.test(stem.dat[stem.dat$stage=="T1","stemness"],stem.dat[stem.dat$stage=="T3","stemness"])
ks.test(stem.dat[stem.dat$stage=="T2","stemness"],stem.dat[stem.dat$stage=="T3","stemness"])

stem.dat.sd=ddply(stem.dat, .(stage), summarize, sd.stem=sd(stemness))
pheatmap(mydata[sci.stem1,rownames(stem.dat)[with(stem.dat, order(stage, stemness))]], cluster_cols = F, annotation_col = ann[,c(1,4)],color = colorRampPalette(c("white","brown"))(100))
# simulated variance 
s=splitn(colnames(mydata), "-",1)
cm.mydata=colMeans(mydata)
s.non.var=foreach(i=1:1000, .combine=rbind) %do% {
  g=sample(1:nrow(mydata), size=121)
  nonstem.dat=data.frame(stage=s,simulate=i, nonstemness=colMeans(mydata[g,])-cm.mydata)
  nonstem.dat.var=ddply(nonstem.dat, .(stage), summarize, var.nonstem=sd(nonstemness))
}
