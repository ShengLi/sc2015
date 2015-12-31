source("monoclewrap.R")
#get TPM from 115 cells
# filter genes by sd
tpm.a.sd=apply(tpm.a, 1, sd)
mytpm = tpm.a[tpm.a.sd > quantile(tpm.a.sd, 0.2),]

# get stage specific DEGs

# get groups vector
getgroups=function(mydata, include){
  stages=splitn(colnames(mydata),"-",1)
  idx=which(stages %in% include)
  mydata1=mydata[,idx]
  stages1=stages[idx]
  groups=ifelse(stages1 == sort(unique(stages1), decreasing = T)[1], 1, 0)
  names(groups)=colnames(mydata1)
  groups
}
# comparing  T1 and T2
groups=getgroups(mytpm, c("T1","T2"))
st12.deg=monoclewrap(mytpm[,names(groups)], groups)

# comparing  T1 and T3
groups=getgroups(mytpm, c("T1","T3"))
st13.deg=monoclewrap(mytpm[,names(groups)], groups)

# comparing  T2 and T3
groups=getgroups(mytpm, c("T2","T3"))
st23.deg=monoclewrap(mytpm[,names(groups)], groups)

# process degs
# log fold change cutoff
cutoff=log2(1.5)
# T1 vs T2
st12.degs.sig=as.character(st12.deg[which(st12.deg$qval <=0.1 & st12.deg$log2foldchange>cutoff),"variable"])
# T1 vs T3
st13.degs.sig=as.character(st13.deg[which(st13.deg$qval <=0.1 &st13.deg$log2foldchange>cutoff),"variable"])
# T2 vs T3
st23.degs.sig=as.character(st23.deg[which(st13.deg$qval <=0.1  & st23.deg$log2foldchange>cutoff),"variable"])

