getGeneID <- function(x, genome){
  ### x: a list of peak positions
  ### genome: reference genome
  ### this function is used to annotate which gene the peak is located on
  require(GenomicRanges)
  require(genomeutils)
  require(GenomicFeatures)
  if(genome == "mm9"){
    require(TxDb.Mmusculus.UCSC.mm9.knownGene)
    require(org.Mm.eg.db)
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
    annodb = org.Mm.eg.db
  }else if(genome == "mm8"){
    require(TxDb.Mmusculus.UCSC.mm8.knownGene)
    require(org.Mm.eg.db)
    txdb <- TxDb.Mmusculus.UCSC.mm8.knownGene
    annodb = org.Mm.eg.db
  }else if(genome == "hg18"){
    require(TxDb.Hsapiens.UCSC.hg18.knownGene)
    require(org.Hs.eg.db)
    txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
    annodb = org.Hs.eg.db
  }else if(genome == "hg19"){
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    require(org.Hs.eg.db)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    annodb = org.Hs.eg.db
  }

  allgene = genes(txdb)
  columns(annodb)
  features <- select(annodb, keys = allgene@elementMetadata@listData$gene_id,
                     columns = c("SYMBOL"), keytype = "ENTREZID")
  genes.GR = GRanges(seqnames = Rle(allgene@seqnames), allgene@ranges, strand = Rle(allgene@strand),
                     gene_id = Rle(allgene$gene_id), geneSymbol = Rle(features$SYMBOL))

  x.GR = GRanges(seqnames = Rle(x$chr), IRanges(x$start, x$end))
  hits = findOverlaps(x.GR, genes.GR)
  query = queryHits(hits)
  subject = subjectHits(hits)
  head(hits)
  ENTREZID = NULL
  Symbol = NULL
  for(i in 1:length(x.GR)){
    idx = which(query==i)
    if (length(idx) == 0){
      ENTREZID[i] = NA
      Symbol[i] = NA
    }else{
      ENTREZID[i] = paste0(genes.GR$gene_id[subject[idx]], collapse = " / ")
      Symbol[i] = paste0(genes.GR$geneSymbol[subject[idx]], collapse = " / ")
    }
  }

  res = data.frame(gene_ENTREZID = ENTREZID, gene_Symbol = Symbol,
                   chr = x.GR@seqnames, start = x.GR@ranges@start,
                   end = x.GR@ranges@start + x.GR@ranges@width -1
  )
  return(res)
}



getPeakSeq <- function(peaks, genome){
  require(GenomicFeatures)
  require(rlist)
  if(genome == "mm9"){
    require(BSgenome.Mmusculus.UCSC.mm9)
    genome = BSgenome.Mmusculus.UCSC.mm9
    allseq = Mmusculus
  }else if(genome == "mm8"){
    require(BSgenome.Mmusculus.UCSC.mm8)
    genome = BSgenome.Mmusculus.UCSC.mm8
    allseq = Mmusculus
  }else if(genome == "hg18"){
    require(BSgenome.Hsapiens.UCSC.hg18)
    genome = BSgenome.Hsapiens.UCSC.hg18
    allseq = Hsapiens
  }else if(genome == "hg19"){
    require(BSgenome.Hsapiens.UCSC.hg19)
    genome = BSgenome.Hsapiens.UCSC.hg19
    allseq = Hsapiens
  }

  ### this function is used to extract the sequence for each peak
  allidx = split(1:nrow(peaks),peaks$chr)
  length(allidx)
  allidx
  allidx[sapply(allidx, length) == 0] <- NULL
  length(allidx)
  allchr = names(allidx)

  DNASets = NULL
  for(ichr in 1:length(allidx)){
    thischr = allchr[ichr]
    thisidx = allidx[[thischr]]
    thisSeq = allseq[[thischr]]
    tmp.set = DNAStringSet(thisSeq, start = peaks$start[thisidx],
                           end = peaks$end[thisidx])
    tmp.ids <- sprintf(paste0(thischr, "_ID%06d"),  thisidx)

    tmp.set = data.frame(name = tmp.ids, seq = as.character(tmp.set))
    DNASets[[ichr]] = tmp.set
  }
  ### combine all sets into one sets
  AllSets = list.rbind(DNASets)
  rownames(AllSets) = AllSets$name
  AllSets$name = NULL
  return(AllSets)
}

getPeakCounts <- function(peaks, allCounts, allBins){
  ### grasp reads counts for a list of peaks
  ### allCounts: counts for all sites across the whole genome
  ### all bins cut from the whole genome
  require(GenomicFeatures)
  peak.GR = GRanges(Rle(peaks$chr), IRanges(peaks$start, peaks$end))
  iii = findOverlaps(peak.GR, allBins)
  query = queryHits(iii)
  subject = subjectHits(iii)
  peakCounts = matrix(0, nrow = length(peak.GR), ncol = ncol(allCounts))
  for(i in 1:length(peak.GR)) {
    ix = query == i
    if(sum(ix) == 0) next
    peakCounts[i,] = colSums(allCounts[subject[ix],,drop = FALSE])
  }
  colnames(peakCounts) = colnames(allCounts)
  return(peakCounts)
}

lfc.peak <- function(peak, sf, binCounts, bins.GR){
  ### peak: genomic sites
  ### sf: size factor for samples that peaks belong to
  ### binCounts: reads counts in all bins across the whloe genome in samples that peaks belong to
  ### bin.GR: Grange object to store the genomic location of each bin
  thiscount = getPeakCounts(peaks = peak, allCounts = binCounts, allBins = bins.GR)
  blocks = seq(1, ncol(binCounts), 2)
  c0 = rep(0, length(blocks))
  for (j in 1:length(blocks)) {
    id = blocks[j]
    dat = binCounts[, id:(id+1)]
    c0[j] = mean(as.matrix(dat), na.rm = TRUE)
  }
  idx.x = seq(1, ncol(binCounts), 2)
  idx.y = seq(2, ncol(binCounts), 2)
  lfc = log(rowMeans((thiscount[, idx.y]/sf[idx.y] + c0)/(thiscount[, idx.x]/sf[idx.x] + c0), na.rm = TRUE))
  smooth.lfc <- mySmooth(lfc, windlen = 5)
  ###
  smooth.lfc
}


mySmooth <- function(vec, windlen = 3){
  ### smooth the value of vec
  cum = cumsum(vec)
  n = length(cum)
  a = c(cum[((windlen - 1)/2 + 1):n], rep(cum[n], (windlen - 1)/2))
  b = c(rep(0, (windlen - 1)/2 + 1), cum[1:(n - (windlen - 1)/2 -1)])
  res = (a-b)/windlen
  return(res)
}









getPeakStrand <- function(bins, genomeBuild){
  ### get bin strand given bin sites and genome
  if(genomeBuild == "mm9"){
    require(TxDb.Mmusculus.UCSC.mm9.knownGene)
    txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
  }else if(genomeBuild == "mm10"){
    require(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  }else if(genomeBuild == "hg19"){
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
  }else if(genomeBuild == "hg18"){
    require(TxDb.Hsapiens.UCSC.hg18.knownGene)
    txdb = TxDb.Hsapiens.UCSC.hg18.knownGene
  }else if(genomeBuild == "rn4"){
    require(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
    txdb = TxDb.Rnorvegicus.UCSC.rn4.ensGene
  }else if(genomeBuild == "dm3"){
    require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
    txdb = TxDb.Dmelanogaster.UCSC.dm3.ensGene
  }

  allGenes = genes(txdb)
  allGenes.strand = as.character(strand(allGenes))

  tmp = findOverlaps(bins, allGenes)
  length(tmp) ## there are bins overlapping more than one gene. Careful
  binIdx = queryHits(tmp)
  geneIdx = subjectHits(tmp)

  tt = table(binIdx)
  binNames = as.integer(names(tt))

  ### start to get strands
  binStrand = rep(".", length(bins))

  ## unique ones
  ix.unique = binNames[which(tt==1)]
  tmpidx = which(binIdx %in% ix.unique)
  geneidx.unique = geneIdx[tmpidx]
  binStrand[ix.unique] = allGenes.strand[geneidx.unique]

  ## deal with bins overlapping multiple genes
  ix.dup = binNames[which(tt>1)]
  tmpidx = which(binIdx %in% ix.dup)
  geneidx.dup = geneIdx[tmpidx]
  strands.dup = allGenes.strand[geneidx.dup]
  binIdx.dup = binIdx[tmpidx]
  res <- tapply(strands.dup, binIdx.dup,
                function(x) {
                  if(all(x == "+")) return("+")
                  else if(all(x == "-")) return("-")
                  else return("*")
                })
  binStrand[ix.dup] = res
  return(binStrand)
}




