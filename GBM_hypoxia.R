# source: http://www.genome.jp/dbget-bin/www_bget?pathway+mmu04066
kegg.hypoxia=read.table("kegg_mmu04066_hypoxia.txt", header=F, sep="\t", stringsAsFactors = F)$V1
hypoxiag=c("Loc102642819", grep("transcription", invert=T, value=T,splitn(grep(";", value=T, kegg.hypoxia),";",1)))
hypoxiag1=intersect(hypoxiag, rownames(mydata))
# calcualte hypoxialevel 
hypoxia.dat=data.frame(stage=splitn(colnames(mydata), "-",1),hypoxia=colMeans(mydata[hypoxiag1,])-colMeans(mydata))
ggplot(hypoxia.dat, aes(x=stage, y=hypoxia)) + geom_boxplot(notch=T)
ggplot(hypoxia.dat, aes(x=hypoxia, colour=stage)) + stat_ecdf()+stage.colour.scale+
  ylab("cumulative distribution function")  + gb2.theme("bottom")
ks.test(hypoxia.dat[hypoxia.dat$stage=="T1","hypoxia"],hypoxia.dat[hypoxia.dat$stage=="T3","hypoxia"])
ks.test(hypoxia.dat[hypoxia.dat$stage=="T2","hypoxia"],hypoxia.dat[hypoxia.dat$stage=="T3","hypoxia"])
ks.test(hypoxia.dat[hypoxia.dat$stage=="T1","hypoxia"],hypoxia.dat[hypoxia.dat$stage=="T2","hypoxia"])
hypoxia.dat.sd=ddply(hypoxia.dat, .(stage), summarize, sd.stem=sd(hypoxia))
# simulated variance 
s=splitn(colnames(mydata), "-",1)
cm.mydata=colMeans(mydata)
h.non.var=foreach(i=1:1000, .combine=rbind) %do% {
  g=sample(1:nrow(mydata), size=length(hypoxiag1))
  nonstem.dat=data.frame(stage=s,simulate=i, nonhypoxia=colMeans(mydata[g,])-cm.mydata)
  nonstem.dat.var=ddply(nonstem.dat, .(stage), summarize, var.nonhypoxia=sd(nonhypoxia))
}
hdat=ddply(h.non.var, .(stage), summarize, mean=mean(var.nonhypoxia), se=se(var.nonhypoxia), type="controls")
hdat=rbind(hdat, data.frame(hypoxia.dat.sd, se=0, type="hypoxia"))
