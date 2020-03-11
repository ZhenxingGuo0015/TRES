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
  if(filetype == "GRanges")
    require(datasetTRES)

  for(i in 1:length(files)) {
    if(filetype!= "GRanges"){
      if(filetype == "bam")
        reads = read.BAM(files[i])
      else if(filetype == "bed")
        reads = import(files[i])
    }else if(filetype == "GRanges"){
      ### data in datasetTRES package
      #reads = get(paste0(getwd(), "/", files[i]))
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

### get all exons given a genome build
### I will not use this for modeling. Only use for peak calling and annotation
getExons <- function(genomeBuild="mm9", binsize=50) {
    if(genomeBuild == "mm9") { ## need to expand this part to include multiple genome builds
        genome = BSgenome.Mmusculus.UCSC.mm9
        txdb = loadDb("mm9_knownGenes.sqlite")
        allExons = exonsBy(txdb, by="gene")
        ## add in other information like GC content

    }
    ## Maybe should add in other information like GC content, etc.
    allExons
}


### prepare equal sized bins on exons
### Note this only needs to be run once for each genome build
### This will return binned windows on each exon
exonBins <- function(genomeBuild = "mm9", binsize = 50) {

  ## get exons
  if(genomeBuild == "mm9") { ## need to expand this part to include multiple genome builds
    require(BSgenome.Mmusculus.UCSC.mm9)
    require(TxDb.Mmusculus.UCSC.mm9.knownGene)
    genome = BSgenome.Mmusculus.UCSC.mm9
    txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
   # txdb = loadDb("mm9_knownGenes.sqlite")
    allExons = exonsBy(txdb, by="gene")
  }else if(genomeBuild == "hg18"){
    require(BSgenome.Hsapiens.UCSC.hg18)
    genome = BSgenome.Hsapiens.UCSC.hg18
    txdb = TxDb.Hsapiens.UCSC.hg18.knownGene
    #txdb = loadDb("hg18_knownGenes.sqlite")
    allExons = exonsBy(txdb, by="gene")
  }else if(genomeBuild == "rn4"){
    require(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
    require(BSgenome.Rnorvegicus.UCSC.rn4)
    genome = BSgenome.Rnorvegicus.UCSC.rn4
    txdb = TxDb.Rnorvegicus.UCSC.rn4.ensGene
    allExons = exonsBy(txdb, by="gene")
  }else if(genomeBuild == "dm3"){
    require(BSgenome.Dmelanogaster.UCSC.dm3)
    require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
    genome = BSgenome.Dmelanogaster.UCSC.dm3
    txdb = TxDb.Dmelanogaster.UCSC.dm3.ensGene
    allExons = exonsBy(txdb, by="gene")
  }

  ## tile up the whole genome
  info = seqinfo(genome)
  ## exclude random chr
  ii = grep("random", seqnames(info))
  if(length(ii) > 0){
    info = info[seqnames(info)[-ii]]
  }else{
    info = info[seqnames(info)]
  }
  wins = tileGenome(info, tilewidth=binsize, cut.last.tile.in.chrom=TRUE)

  ## get bins on the exons
  ix = wins %over% allExons
  wins.exons = wins[ix]

  ## return
  invisible(return(list(bins=wins.exons, allExons=allExons)))
}


getStrand <- function(bins, genomeBuild){
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

