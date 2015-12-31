monoclewrap=function(mytpm, groups, cores=5){
  library(monocle)
  library(reshape)
  library(plyr)
  # gene name
  gene_annotation=data.frame(gene_short_name=rownames(mytpm), 
                             biotype="protein_coding",num_cells_expressed=0, use_for_ordering=FALSE)
  gene_annotation$num_cells_expressed=rowSums(mytpm>0.1)
  #gene_annotation$use_for_ordering[which(gene_annotation$gene_short_name %in% ordering_genes)]=rep(TRUE, 3)
  rownames(gene_annotation)=gene_annotation$gene_short_name
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  
  
  pdm=data.frame(Library=colnames(mytpm), group=groups)
  pd <- new("AnnotatedDataFrame", data = pdm)
  
  HSMM <- newCellDataSet(as.matrix(mytpm), phenoData = pd, featureData = fd)
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
  diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr="expression~group", cores = cores)
  
  #sig_degs=diff_test_res[diff_test_res$qval<0.05,]
  sig_degs1=diff_test_res[with(diff_test_res, order(pval)),]
  log2tpm=log2(1+mytpm[rownames(sig_degs1),])
  pdm=data.frame(Library=colnames(log2tpm), group=groups)
  dat=cbind(data.frame(group=as.character(pdm$group), data.frame(t(log2tpm))))
  mdat=melt(dat)
  mdat1=ddply(mdat, .( variable, group), summarize, mean=mean(value))
  mdat2=cast(mdat1, variable~group, values=mean)
  mdat3=cbind(mdat2, sig_degs1[gsub("\\.", "-", as.character(mdat2$variable)),2:3])
  mdat3$log2foldchange=mdat3[,"1"] - mdat3[,"0"]
  mdat3
  #mdat4=mdat3[mdat3$log2foldchange >= 1.5,]
  #mdat5=mdat4[order(mdat4$log2foldchange, decreasing = T),]
  #rownames(mdat5)=1:nrow(mdat5)
  
}

# without cutoff
getlogfc=function(mytpm, groups, cores=5){
  # filter by exprs and cell number
  idx1=which(rowMeans(mytpm) >= 0.1)
  mytpm1=mytpm[idx1,]
  idx2=which(rowSums(mytpm1 >= 0.1) >= 20)
  mytpm2=mytpm1[idx2,]
  pdm=data.frame(Library=colnames(mytpm2), group=groups)
  dat=cbind(data.frame(group=as.character(pdm$group), data.frame(t(log2tpm))))
  mdat=melt(dat)
  mdat1=ddply(mdat, .( variable, group), summarize, mean=mean(value))
  mdat2=cast(mdat1, variable~group, values=mean)
  mdat3=cbind(mdat2, sig_degs1[as.character(mdat2$variable),2:3])
  mdat3$log2foldchange=mdat3[,"1"] - mdat3[,"0"]
  mdat3
  #mdat4=mdat3[mdat3$log2foldchange >= 1.5,]
  #mdat5=mdat4[order(mdat4$log2foldchange, decreasing = T),]
  #rownames(mdat5)=1:nrow(mdat5)
  
}