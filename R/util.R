## function to read BAM file and create GRanges for the reads
### read in the BAM file and create GRanges for all reads
read.BAM <- function(fn){
  require(Rsamtools)
  param = ScanBamParam(what=c("rname","strand","pos","qwidth"))
  bam = scanBam(fn, param=param)[[1]]
  ix = !is.na(bam$rname) & !is.na(bam$pos)
  qwidth = bam$qwidth[ix]
  IRange.reads <- GRanges(seqnames=Rle(bam$rname[ix]),
                          ranges=IRanges(bam$pos[ix], width=bam$qwidth[ix]),
                          strand=Rle(bam$strand[ix]))
  IRange.reads
}


## get read counts for given genomic windows
getWinCounts <- function(files, wins, filetype = c("bed", "bam", "GRanges")) {
  if(class(wins)!="data.frame" & class(wins)!="GRanges")
    stop("Input genomic intervals must be a GRanges or data frame!")
  if(class(wins)=="data.frame") wins = import(wins)
  counts = matrix(0, nrow = length(wins), ncol = length(files))

  for(i in 1:length(files)) {
    if(filetype!= "GRanges"){
      if(filetype == "bam")
        reads = read.BAM(files[i])
      else if(filetype == "bed")
        reads = import(files[i])
    }else if(filetype == "GRanges"){
      ### data in datasetTRES package
      reads = get(files[i])
    }
    ### added my Zhenxing to avoid repeated counting
    # width(reads) = 2
    #### ended
    counts[,i] = countOverlaps(wins,reads)
  }
  colnames(counts) = files
  counts
}


exonBins <- function(genomeBuild = "mm9", anno_TXDB = NULL, binsize = 50,
                     IncludeIntron = FALSE) {

  ## get exons and introns
  if(genomeBuild == "bosTau8"){
    require(BSgenome.Btaurus.UCSC.bosTau8)
    genome = BSgenome.Btaurus.UCSC.bosTau8
    if(length(anno_TXDB) == 0){
      require(TxDb.Btaurus.UCSC.bosTau8.refGene)
      txdb = TxDb.Btaurus.UCSC.bosTau8.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "bosTau9"){
    require(BSgenome.Btaurus.UCSC.bosTau9)
    genome = BSgenome.Btaurus.UCSC.bosTau9
    if(length(anno_TXDB) == 0){
      require(TxDb.Btaurus.UCSC.bosTau9.refGene)
      txdb = TxDb.Btaurus.UCSC.bosTau9.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "ce6"){
    require(BSgenome.Celegans.UCSC.ce6)
    genome = BSgenome.Celegans.UCSC.ce6
    if(length(anno_TXDB) == 0){
      require(TxDb.Celegans.UCSC.ce6.ensGene)
      txdb = TxDb.Celegans.UCSC.ce6.ensGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "ce11"){
    require(BSgenome.Celegans.UCSC.ce11)
    genome = BSgenome.Celegans.UCSC.ce11
    if(length(anno_TXDB) == 0){
      require(TxDb.Celegans.UCSC.ce11.ensGene)
      txdb = TxDb.Celegans.UCSC.ce11.ensGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "canFam3"){
    require(BSgenome.Cfamiliaris.UCSC.canFam3)
    genome = BSgenome.Cfamiliaris.UCSC.canFam3
    if(length(anno_TXDB) == 0){
      require(TxDb.Cfamiliaris.UCSC.canFam3.refGene)
      txdb = TxDb.Cfamiliaris.UCSC.canFam3.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "dm3"){
    require(BSgenome.Dmelanogaster.UCSC.dm3)
    genome = BSgenome.Dmelanogaster.UCSC.dm3
    if(length(anno_TXDB) == 0){
      require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
      txdb = TxDb.Dmelanogaster.UCSC.dm3.ensGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "dm6"){
    require(BSgenome.Dmelanogaster.UCSC.dm6)
    genome = BSgenome.Dmelanogaster.UCSC.dm6
    if(length(anno_TXDB) == 0){
      require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
      txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "danRer10"){
    require(BSgenome.Drerio.UCSC.danRer10)
    genome = BSgenome.Drerio.UCSC.danRer10
    if(length(anno_TXDB) == 0){
      require(TxDb.Drerio.UCSC.danRer10.refGene)
      txdb = TxDb.Drerio.UCSC.danRer10.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "danRer11"){
    require(BSgenome.Drerio.UCSC.danRer11)
    genome = BSgenome.Drerio.UCSC.danRer11
    if(length(anno_TXDB) == 0){
      require(TxDb.Drerio.UCSC.danRer11.refGene)
      txdb = TxDb.Drerio.UCSC.danRer11.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "galGal4"){
    require(BSgenome.Ggallus.UCSC.galGal4)
    genome = BSgenome.Ggallus.UCSC.galGal4
    if(length(anno_TXDB) == 0){
      require(TxDb.Ggallus.UCSC.galGal4.refGene)
      txdb = TxDb.Ggallus.UCSC.galGal4.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "galGal5"){
    require(BSgenome.Ggallus.UCSC.galGal5)
    genome = BSgenome.Ggallus.UCSC.galGal5
    if(length(anno_TXDB) == 0){
      require(TxDb.Ggallus.UCSC.galGal5.refGene)
      txdb = TxDb.Ggallus.UCSC.galGal5.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "galGal6"){
    require(BSgenome.Ggallus.UCSC.galGal6)
    genome = BSgenome.Ggallus.UCSC.galGal6
    if(length(anno_TXDB) == 0){
      require(TxDb.Ggallus.UCSC.galGal6.refGene)
      txdb = TxDb.Ggallus.UCSC.galGal6.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "hg17"){
    require(BSgenome.Hsapiens.UCSC.hg17)
    genome = BSgenome.Hsapiens.UCSC.hg17
    if(length(anno_TXDB) == 0){
      require(TxDb.Hsapiens.UCSC.hg17.knownGene)
      txdb = TxDb.Hsapiens.UCSC.hg17.knownGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "hg18"){
    require(BSgenome.Hsapiens.UCSC.hg18)
    genome = BSgenome.Hsapiens.UCSC.hg18
    if(length(anno_TXDB) == 0){
      require(TxDb.Hsapiens.UCSC.hg18.knownGene)
      txdb = TxDb.Hsapiens.UCSC.hg18.knownGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "hg19"){
    require(BSgenome.Hsapiens.UCSC.hg19)
    genome = BSgenome.Hsapiens.UCSC.hg19
    if(length(anno_TXDB) == 0){
      require(TxDb.Hsapiens.UCSC.hg19.knownGene)
      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "hg38"){
    require(BSgenome.Hsapiens.UCSC.hg38)
    genome = BSgenome.Hsapiens.UCSC.hg38
    if(length(anno_TXDB) == 0){
      require(TxDb.Hsapiens.UCSC.hg38.knownGene)
      txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "mm8"){
    require(BSgenome.Mmusculus.UCSC.mm8)
    genome = BSgenome.Mmusculus.UCSC.mm8
    if(length(anno_TXDB) == 0){
      require(TxDb.Mmusculus.UCSC.mm8.knownGene)
      txdb = TxDb.Mmusculus.UCSC.mm8.knownGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "mm9") {
    require(BSgenome.Mmusculus.UCSC.mm9)
    genome = BSgenome.Mmusculus.UCSC.mm9
    if(length(anno_TXDB) == 0){
      require(TxDb.Mmusculus.UCSC.mm9.knownGene)
      txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "mm10"){
    require(BSgenome.Mmusculus.UCSC.mm10)
    genome = BSgenome.Mmusculus.UCSC.mm10
    if(length(anno_TXDB) == 0){
      require(TxDb.Mmusculus.UCSC.mm10.knownGene)
      txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "panTro5"){
    require(BSgenome.Ptroglodytes.UCSC.panTro5)
    genome = BSgenome.Ptroglodytes.UCSC.panTro5
    if(length(anno_TXDB) == 0){
      require(TxDb.Ptroglodytes.UCSC.panTro5.refGene)
      txdb = TxDb.Ptroglodytes.UCSC.panTro5.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "panTro6"){
    require(BSgenome.Ptroglodytes.UCSC.panTro6)
    genome = BSgenome.Ptroglodytes.UCSC.panTro6
    if(length(anno_TXDB) == 0){
      require(TxDb.Ptroglodytes.UCSC.panTro6.refGene)
      txdb = TxDb.Ptroglodytes.UCSC.panTro6.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "rn4"){
    require(BSgenome.Rnorvegicus.UCSC.rn4)
    genome = BSgenome.Rnorvegicus.UCSC.rn4
    if(length(anno_TXDB) == 0){
      require(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
      txdb = TxDb.Rnorvegicus.UCSC.rn4.ensGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "rn5"){
    require(BSgenome.Rnorvegicus.UCSC.rn5)
    genome = BSgenome.Rnorvegicus.UCSC.rn5
    if(length(anno_TXDB) == 0){
      require(TxDb.Rnorvegicus.UCSC.rn5.refGene)
      txdb = TxDb.Rnorvegicus.UCSC.rn5.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "rn6"){
    require(BSgenome.Rnorvegicus.UCSC.rn6)
    genome = BSgenome.Rnorvegicus.UCSC.rn6
    if(length(anno_TXDB) == 0){
      require(TxDb.Rnorvegicus.UCSC.rn6.refGene)
      txdb = TxDb.Rnorvegicus.UCSC.rn6.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "rheMac3"){
    require(BSgenome.Mmulatta.UCSC.rheMac3)
    genome = BSgenome.Mmulatta.UCSC.rheMac3
    if(length(anno_TXDB) == 0){
      require(TxDb.Mmulatta.UCSC.rheMac3.refGene)
      txdb = TxDb.Mmulatta.UCSC.rheMac3.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "rheMac8"){
    require(BSgenome.Mmulatta.UCSC.rheMac8)
    genome = BSgenome.Mmulatta.UCSC.rheMac8
    if(length(anno_TXDB) == 0){
      require(TxDb.Mmulatta.UCSC.rheMac8.refGene)
      txdb = TxDb.Mmulatta.UCSC.rheMac8.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "rheMac10"){
    require(BSgenome.Mmulatta.UCSC.rheMac10)
    genome = BSgenome.Mmulatta.UCSC.rheMac10
    if(length(anno_TXDB) == 0){
      require(TxDb.Mmulatta.UCSC.rheMac10.refGene)
      txdb = TxDb.Mmulatta.UCSC.rheMac10.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "sacCer2"){
    require(BSgenome.Scerevisiae.UCSC.sacCer2)
    genome = BSgenome.Scerevisiae.UCSC.sacCer2
    if(length(anno_TXDB) == 0){
      require(TxDb.Scerevisiae.UCSC.sacCer2.sgdGene)
      txdb = TxDb.Scerevisiae.UCSC.sacCer2.sgdGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "sacCer3"){
    require(BSgenome.Scerevisiae.UCSC.sacCer3)
    genome = BSgenome.Scerevisiae.UCSC.sacCer3
    if(length(anno_TXDB) == 0){
      require(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
      txdb = TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "susScr3"){
    require(BSgenome.Sscrofa.UCSC.susScr3)
    genome = BSgenome.Scerevisiae.UCSC.sacCer3
    if(length(anno_TXDB) == 0){
      require(TxDb.Sscrofa.UCSC.susScr3.refGene)
      txdb = TxDb.Sscrofa.UCSC.susScr3.refGene
    }else{
      txdb = anno_TXDB
    }
  }else if(genomeBuild == "susScr11"){
    require(BSgenome.Sscrofa.UCSC.susScr11)
    genome = BSgenome.Scerevisiae.UCSC.sacCer11
    if(length(anno_TXDB) == 0){
      require(TxDb.Sscrofa.UCSC.susScr11.refGene)
      txdb = TxDb.Sscrofa.UCSC.susScr11.refGene
    }else{
      txdb = anno_TXDB
    }
  }

  allExons = exonsBy(txdb, by="gene")
  allIntron = intronsByTranscript(txdb)
  ## tile up the whole genome
  require(GenomeInfoDb)
  info = seqinfo(genome)
  ## exclude random chr
  ii = grep("random", seqnames(info))
  if(length(ii) > 0){
    info = info[seqnames(info)[-ii]]
  }else{
    info = info[seqnames(info)]
  }
  wins = tileGenome(info, tilewidth=binsize, cut.last.tile.in.chrom=TRUE)

  ###
  if(!IncludeIntron){
    ## get bins on the exons
    ix = wins %over% allExons
    wins.exons = wins[ix]
    ## return
    invisible(return(list(bins=wins.exons, keep.bin = ix, allExons=allExons)))
  }else{
    ix = (wins %over% allExons) | (wins %over% allIntron)
    wins.exIn = wins[ix]
    invisible(return(list(bins=wins.exIn,
                          keep.bin = ix,
                          allExons=allExons,
                          allIntrons = allIntron)))
  }
}


getStrand <- function(bins, genomeBuild, anno_TXDB){
  ### get bin strand given bin sites and genome
  if(length(anno_TXDB) > 0){
    txdb = anno_TXDB
  }else{
    if(genomeBuild == "bosTau8"){
      require(TxDb.Btaurus.UCSC.bosTau8.refGene)
      txdb = TxDb.Btaurus.UCSC.bosTau8.refGene
    }else if(genomeBuild == "bosTau9"){
      require(TxDb.Btaurus.UCSC.bosTau9.refGene)
      txdb = TxDb.Btaurus.UCSC.bosTau9.refGene
    }else if(genomeBuild == "ce6"){
      require(TxDb.Celegans.UCSC.ce6.ensGene)
      txdb = TxDb.Celegans.UCSC.ce6.ensGene
    }else if(genomeBuild == "ce11"){
      require(TxDb.Celegans.UCSC.ce11.ensGene)
      txdb = TxDb.Celegans.UCSC.ce11.ensGene
    }else if(genomeBuild == "canFam3"){
      require(TxDb.Cfamiliaris.UCSC.canFam3.refGene)
      txdb = TxDb.Cfamiliaris.UCSC.canFam3.refGene
    }else if(genomeBuild == "dm3"){
      require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
      txdb = TxDb.Dmelanogaster.UCSC.dm3.ensGene
    }else if(genomeBuild == "dm6"){
      require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
      txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene
    }else if(genomeBuild == "danRer10"){
      require(TxDb.Drerio.UCSC.danRer10.refGene)
      txdb = TxDb.Drerio.UCSC.danRer10.refGene
    }else if(genomeBuild == "danRer11"){
      require(TxDb.Drerio.UCSC.danRer11.refGene)
      txdb = TxDb.Drerio.UCSC.danRer11.refGene
    }else if(genomeBuild == "galGal4"){
      require(TxDb.Ggallus.UCSC.galGal4.refGene)
      txdb = TxDb.Ggallus.UCSC.galGal4.refGene
    }else if(genomeBuild == "galGal5"){
      require(TxDb.Ggallus.UCSC.galGal5.refGene)
      txdb = TxDb.Ggallus.UCSC.galGal5.refGene
    }else if(genomeBuild == "galGal6"){
      require(TxDb.Ggallus.UCSC.galGal6.refGene)
      txdb = TxDb.Ggallus.UCSC.galGal6.refGene
    }else if(genomeBuild == "hg17"){
      require(TxDb.Hsapiens.UCSC.hg17.knownGene)
      txdb = TxDb.Hsapiens.UCSC.hg17.knownGene
    }else if(genomeBuild == "hg18"){
      require(TxDb.Hsapiens.UCSC.hg18.knownGene)
      txdb = TxDb.Hsapiens.UCSC.hg18.knownGene
    }else if(genomeBuild == "hg19"){
      require(TxDb.Hsapiens.UCSC.hg19.knownGene)
      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
    }else if(genomeBuild == "hg38"){
      require(TxDb.Hsapiens.UCSC.hg38.knownGene)
      txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
    }else if(genomeBuild == "mm8"){
      require(TxDb.Mmusculus.UCSC.mm8.knownGene)
      txdb = TxDb.Mmusculus.UCSC.mm8.knownGene
    }else if(genomeBuild == "mm9") {
      require(TxDb.Mmusculus.UCSC.mm9.knownGene)
      txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
    }else if(genomeBuild == "mm10"){
      require(TxDb.Mmusculus.UCSC.mm10.knownGene)
      txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
    }else if(genomeBuild == "panTro5"){
      require(TxDb.Ptroglodytes.UCSC.panTro5.refGene)
      txdb = TxDb.Ptroglodytes.UCSC.panTro5.refGene
    }else if(genomeBuild == "panTro6"){
      require(TxDb.Ptroglodytes.UCSC.panTro6.refGene)
      txdb = TxDb.Ptroglodytes.UCSC.panTro6.refGene
    }else if(genomeBuild == "rn4"){
      require(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
      txdb = TxDb.Rnorvegicus.UCSC.rn4.ensGene
    }else if(genomeBuild == "rn5"){
      require(TxDb.Rnorvegicus.UCSC.rn5.refGene)
      txdb = TxDb.Rnorvegicus.UCSC.rn5.refGene
    }else if(genomeBuild == "rn6"){
      require(TxDb.Rnorvegicus.UCSC.rn6.refGene)
      txdb = TxDb.Rnorvegicus.UCSC.rn6.refGene
    }else if(genomeBuild == "rheMac3"){
      require(TxDb.Mmulatta.UCSC.rheMac3.refGene)
      txdb = TxDb.Mmulatta.UCSC.rheMac3.refGene
    }else if(genomeBuild == "rheMac8"){
      require(TxDb.Mmulatta.UCSC.rheMac8.refGene)
      txdb = TxDb.Mmulatta.UCSC.rheMac8.refGene
    }else if(genomeBuild == "rheMac10"){
      require(TxDb.Mmulatta.UCSC.rheMac10.refGene)
      txdb = TxDb.Mmulatta.UCSC.rheMac10.refGene
    }else if(genomeBuild == "sacCer2"){
      require(TxDb.Scerevisiae.UCSC.sacCer2.sgdGene)
      txdb = TxDb.Scerevisiae.UCSC.sacCer2.sgdGene
    }else if(genomeBuild == "sacCer3"){
      require(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
      txdb = TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
    }else if(genomeBuild == "susScr3"){
      require(TxDb.Sscrofa.UCSC.susScr3.refGene)
      txdb = TxDb.Sscrofa.UCSC.susScr3.refGene
    }else if(genomeBuild == "susScr11"){
      require(TxDb.Sscrofa.UCSC.susScr11.refGene)
      txdb = TxDb.Sscrofa.UCSC.susScr11.refGene
    }
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


divid.into.bins <- function(TXDB, binsize = 50, IncludeIntron = FALSE){
  #### each transcript will be divided into bins of length "binsize"
  ### when the annotation file are provided by user instead of downloading from internet
  tmp = unique(transcripts(TXDB))
  idx <- split(1:length(as.data.frame(tmp)$seqnames), as.data.frame(tmp)$seqnames)
  iii = (grepl("random", names(idx)) | grepl("M", names(idx)))
  idx[iii] = NULL

  allchr = names(idx)
  allbins = NULL
  for (ichr in 1:length(allchr)) {
    cat(allchr[ichr], sep = "\t")
    thischr = allchr[ichr]
    thiset  = idx[[ichr]]
    thiset.strand = as.character(strand(tmp[thiset]))
    ### split negative and positive bins
    ### positive
    ii.pos = which(thiset.strand == "+")
    pos.min = min(tmp[thiset[ii.pos]]@ranges@start)
    pos.max = max(tmp[thiset[ii.pos]]@ranges@start + tmp[thiset[ii.pos]]@ranges@width -1)

    posbins.start = seq(pos.min, pos.max, by = binsize)
    posbins.end = c(posbins.start[-1] - 1,pos.max)
    posbin.strand = unique(tmp@strand[thiset[ii.pos]])
    pos.bin = GRanges(Rle(thischr), IRanges(posbins.start, posbins.end), Rle("+"))

    ### negative
    ii.neg = which(thiset.strand == "-")
    neg.min = min(tmp[thiset[ii.neg]]@ranges@start)
    neg.max = max(tmp[thiset[ii.neg]]@ranges@start + tmp[thiset[ii.neg]]@ranges@width - 1)

    negbins.start = seq(neg.min, neg.max, by = binsize)
    negbins.end = c(negbins.start[-1] - 1, neg.max)
    negbin.strand = unique(tmp@strand[thiset[ii.neg]])
    neg.bin = GRanges(Rle(thischr), IRanges(negbins.start, negbins.end), Rle("-"))

    #### combine positive and negtive bins
    thisbins = c(pos.bin, neg.bin)
    thisbins = as.data.frame(thisbins)

    if(length(allbins) == 0){
      allbins = thisbins
    }else{
      allbins = rbind(allbins, thisbins)
    }

  }

  allbins = GRanges(Rle(allbins$seqnames), IRanges(allbins$start, allbins$end), Rle(allbins$strand))

  ####
  allExons = unlist(exonsBy(TXDB))
  allIntron = unlist(intronsByTranscript(TXDB))

  if(!IncludeIntron){
    ## get bins on the exons
    ix = allbins %over% allExons
    bins.exons = allbins[ix]
    ## return
    #res = list(bins=bins.exons, keep.bin = ix, allExons=allExons)
    invisible(return(list(bins=bins.exons, keep.bin = ix, allExons=allExons)))
  }else{
    ix = (allbins %over% allExons) | (allbins %over% allIntron)
    bins.exIn = allbins[ix]
    invisible(return(list(bins=bins.exIn,
                          keep.bin = ix,
                          allExons=allExons,
                          allIntrons = allIntron)))
  }

}
