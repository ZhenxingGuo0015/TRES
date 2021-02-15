### RNASeq tools: DESeq2
m6Apeak.RNAseq <- function(counts, sf = NULL) {
  ### DESeq2
  require(DESeq2)
  require(stats)
  nreps = ncol(counts)/2
  design = data.frame(Trt = rep(c("Input", "IP"), nreps),
                      Reps = c(rep("Rep1", nreps), rep("Rep2", nreps)))
  model = ~ Reps + Trt
  model.matrix(model, design)

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = design,
                                design = model)

  if(length(sf) ==  0){
    ## compute size factor from all counts
    sf = colSums(counts)
    sf = sf/median(sf)
  }
  counts.norm = sweep(counts, 2, sf, FUN = "/")
  sizeFactors(dds) = sf
  dds <- DESeq(dds)

  resultsNames(dds) # name of each variable in the design matrix used for DESeq
  res = results(dds, name= "Trt_IP_vs_Input")
  # res$pvalue[res$log2FoldChange < 0] = 1
  # ix <- order(res$pvalue)
  # res1 = res[ix, ]
  # return(res1$pvalue)
  dat = data.frame(pval = res$pvalue, stat = res$stat)
  return(dat)
}
