# source:http://www.genome.jp/dbget-bin/www_bget?mmu:102642920
kegg.cycle=read.table("kegg_mmu04110_cellcycle.txt", header=F, sep="\t", stringsAsFactors = F)$V1
cycleg=c("Loc102642920", grep("cut9", invert=T, value=T,splitn(grep(";", value=T, kegg.cycle),";",1)))
cycleg1=intersect(cycleg, rownames(mydata))
# calcualte cyclelevel 
cycle.dat=data.frame(stage=splitn(colnames(mydata), "-",1),cellcycle=colMeans(mydata[cycleg1,])-colMeans(mydata))
ggplot(cycle.dat, aes(x=cellcycle, colour=stage)) + stat_ecdf()+stage.colour.scale+
  ylab("cumulative distribution function") + xlab("cell cycle") + gb2.theme("bottom")
ks.test(cycle.dat[cycle.dat$stage=="T1","cellcycle"],cycle.dat[cycle.dat$stage=="T3","cellcycle"])
ks.test(cycle.dat[cycle.dat$stage=="T1","cellcycle"],cycle.dat[cycle.dat$stage=="T2","cellcycle"])
ks.test(cycle.dat[cycle.dat$stage=="T2","cellcycle"],cycle.dat[cycle.dat$stage=="T3","cellcycle"])
cellcycle.dat.sd=ddply(cycle.dat, .(stage), summarize, sd.stem=sd(cellcycle))
# simulated variance 
s=splitn(colnames(mydata), "-",1)
cm.mydata=colMeans(mydata)
c.non.var=foreach(i=1:1000, .combine=rbind) %do% {
  g=sample(1:nrow(mydata), size=50)
  nonstem.dat=data.frame(stage=s,simulate=i, noncycle=colMeans(mydata[g,])-cm.mydata)
  nonstem.dat.var=ddply(nonstem.dat, .(stage), summarize, var.noncycle=sd(noncycle))
}
colnames(cellcycle.dat.sd)[2]="mean"
cdat=ddply(c.non.var, .(stage), summarize, mean=mean(var.noncycle), se=se(var.noncycle), type="controls")
cdat=rbind(cdat, data.frame(cellcycle.dat.sd, se=0, type="cell cycle"))
cdat$type=factor(cdat$type, levels=c("controls", "cell cycle"))
