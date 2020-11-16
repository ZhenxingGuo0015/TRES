
TRES_peak <- function(IP.file, Input.file,
                     genomeBuild = NULL,
                     binsize = 50,
                     sf0 = NULL,
                     WhichThreshold = "fdr_lfc",
                     pval.cutoff0 = 1e-5,
                     fdr.cutoff0 = 0.05,
                     lfc.cutoff0 = 0.7,
                     InputDir = NULL,
                     OutputDir = NULL,
                     experiment_name,
                     filetype = "bam"){
  ### 1. divide the genome into bins and get bin counts, bin directions
  t.0 = Sys.time()
  t1 = Sys.time()
  cat(paste0("Divid genome of ", genomeBuild, " into bins..." ), sep = "\n")
  bins.exons = exonBins(genomeBuild = genomeBuild, binsize = binsize)
  binStrand = getStrand(genomeBuild = genomeBuild, bins = bins.exons$bins)
  t2 = Sys.time()
  t2 - t1
  cat("Time used for dividing genome is: ", t2 - t1, sep = "\n")

  datafiles = rep(NA, 2*length(IP.file))
  datafiles[seq(1, length(datafiles), 2)] = Input.file
  datafiles[seq(2, length(datafiles), 2)] = IP.file

  cat("Get bin counts...", sep = "\n")
  t1 =  Sys.time()
  allCounts = getWinCounts(files = file.path(paste0(InputDir, datafiles)),
                           wins = bins.exons$bins, filetype = filetype)
  t2 - t1
  t2 =  Sys.time()
  cat("Time used for obtaining bin counts is: ", t2 - t1, sep = "\n")


  ### 2. peak calling
  if(is.null(sf0))
    sf0 = colSums(allCounts)/median(colSums(allCounts))
  allBins = as.data.frame(bins.exons$bins)
  colnames(allBins)[1] = "chr"
  allBins$strand = binStrand

  cat("Start to call peaks...", sep = "\n")
  if(length(IP.file) > 1){
    ##### two-step procedure

    ### step 1
    cat("###### Step 1:...", sep = "\n")
    t1 = Sys.time()
    Peak.candidates = M6Apeak.MultiRep.step1(Counts = allCounts,
                                             bins = allBins,
                                             sf = sf0,
                                             WhichThreshold = WhichThreshold,
                                             fdr.cutoff = fdr.cutoff0,
                                             lfc.cutoff = lfc.cutoff0)
    t2 = Sys.time()
    t2 - t1
    cat("Time used in Step 1 is: ", t2 - t1, sep = "\n")

    if(nrow(Peak.candidates) >= 2){
      ###  step 2
      cat("###### Step 2:...", sep = "\n")
      t1 = Sys.time()

      ### added on April 14, 2020: estimate background methylation level bg.mu
      idx = which(grepl("rep", colnames(Peak.candidates)) | grepl("bam", colnames(Peak.candidates)))
      PeakCount = Peak.candidates[, idx]
      bgCount = colSums(allCounts) - colSums(PeakCount)
      bg.Input = bgCount[seq(1, length(bgCount), 2)]
      bg.IP = bgCount[seq(2, length(bgCount), 2)]
      bg.mu = mean((bg.IP/sf0[seq(2, length(bgCount), 2)] )/(bg.IP/sf0[seq(2, length(bgCount), 2)]
                                                        + bg.Input/sf0[seq(1, length(bgCount), 2)]),
                   na.rm = TRUE)
      #### stop editing

      Peaks = M6Apeak.MultiRep.step2(Candidates = Peak.candidates,
                                     sf = sf0,
                                     mu.cutoff = bg.mu)
      t2 = Sys.time()
      t2 - t1
      cat("Time used in Step 2 is: ", t2 - t1, sep = "\n")
    }else{
      Peaks = Peak.candidates
      cat("Less than 2 candidates!", sep = "\n")
    }
    ### save results
    save(bins.exons, binStrand, allCounts,
         file = paste0(OutputDir, experiment_name, ".rda"))
    if(!is.null(Peaks$width)){
      Peaks$width = NULL
    }
    write.table(Peaks, file = paste0(OutputDir, experiment_name, "_peaks.xls"),
                row.names = FALSE, sep = "\t")
  }else if(length(IP.file) == 1){
    Peaks = M6Apeak.oneRep(Counts = allCounts,
                           bins = allBins,
                           sf = sf0,
                           WhichThreshold = WhichThreshold,
                           pval.cutoff = pval.cutoff0,
                           fdr.cutoff = fdr.cutoff0,
                           lfc.cutoff = lfc.cutoff0,
                           lowCount = 10)
    ### save results
    save(bins.exons, binStrand, allCounts,
         file = paste0(OutputDir, experiment_name, ".rda"))
    if(!is.null(Peaks$width)){
      Peaks$width = NULL
    }
    write.table(Peaks, file = paste0(OutputDir, experiment_name, "_peaks.xls"),
                row.names = FALSE, sep = "\t")
  }
  cat("###### Done!", sep = "\n")
  ###
  t.1 = Sys.time()
  cat("Total time used for peak calling is: ", sep = "\n")
  t.1 - t.0
}


M6Apeak.MultiRep.step1 <- function(Counts, sf, bins,
                                   WhichThreshold = "fdr",
                                   pval.cutoff = 1e-5,
                                   fdr.cutoff = 1e-05,
                                   lfc.cutoff = 0.7,
                                   windlen = 5,
                                   lowcount = 50){
  ### peak calling when there are more than one replicate
  require(GenomicFeatures)
  require(stats)
  ### 1. find bumps for each replicate based on binomial test
  sx = sf[seq(1, length(sf), 2)]
  sy = sf[seq(2, length(sf), 2)]
  bn.prob = sy/(sx + sy)

  blocks = seq(1, ncol(Counts), 2)
  c0 = rep(0, length(blocks))
  Pvals = matrix(0, nrow = nrow(Counts), ncol = length(blocks))
  Bumps = vector("list", length = length(blocks))
  for (j in 1:length(blocks)) {
    cat(j, sep = "\n")
    id = blocks[j]
    dat = Counts[, id:(id+1)]
    thissf = sf[id:(id+1)]
    ### pvals based on binomial test
    idx = rowSums(dat) > 0
    Pvals[idx, j] = 1 - pbinom(dat[idx, 2], rowSums(dat[idx, ]), prob = bn.prob[j])
    Pvals[!idx, j] = 1

    ### lfc
    c0[j] = mean(as.matrix(dat), na.rm = TRUE)  ### pseudocount
    lfc = log((dat[, 2]/thissf[2] + c0[j])/(dat[, 1]/thissf[1] + c0[j]))   ### modified as this on Nov 19, 2019
    smooth.lfc <- mySmooth(lfc, windlen = windlen)
    x.vals = data.frame(pvals = Pvals[, j],
                        fdr = p.adjust(Pvals[, j], method = "fdr"),
                        lfc = lfc)

    ### find bumps based on pvals, fdr or lfc
    tmp = findBumps(chr = bins$chr, pos = bins$start, strand = bins$strand,
                    x = x.vals,
                    use = WhichThreshold,
                    pval.cutoff = pval.cutoff,
                    fdr.cutoff = fdr.cutoff,
                    lfc.cutoff = lfc.cutoff,
                    count = dat)
    Bumps[[j]] = tmp[tmp$counts > 10, ] ### remove very low count
  }

  ### 2. merge bumps from different replicates
  cat("Merge bumps from different replicates...", sep = "\n")
  thispeak = NULL
  for (i in 1:(length(Bumps))) {
    cat(paste0("Bumps ", i), sep = "\n")
    thisBump = Bumps[[i]]
    thisBump$strand[thisBump$strand=="."] ="*"
    thisBump.GR = GRanges(Rle(thisBump$chr), IRanges(thisBump$start, thisBump$end), Rle(thisBump$strand))
    if(is.null(thispeak)){
      thispeak = thisBump.GR
    }else{
     # require(GenomicRanges)
      thispeak = GenomicRanges::union(thispeak, thisBump.GR)
    }
  }
  thispeak = as.data.frame(thispeak)
  colnames(thispeak) = c("chr", "start", "end", "width", "strand")

  if(nrow(thispeak) >= 1){
    ### add summit in bumps to peak
    peak.summit = rep(NA, nrow(thispeak))
    for (i in 1:length(Bumps)) {
      thisBump = Bumps[[i]]
      thisBump$strand[thisBump$strand=="."] = "*"
      bump.GR = GRanges(Rle(thisBump$chr), IRanges(thisBump$start, thisBump$end), Rle(thisBump$strand))
      peak.GR = GRanges(Rle(thispeak$chr), IRanges(thispeak$start, thispeak$end), Rle(thispeak$strand))
      hits = findOverlaps(peak.GR, bump.GR)
      if(length(hits) > 0){
        idx.query = unique(queryHits(hits))
        idx.subject = subjectHits(hits)
        for (j in 1:length(idx.query)) {
          thisquery = idx.query[j]
          thissubject = idx.subject[which(queryHits(hits) == thisquery)]
          if(is.na(peak.summit[idx.query[j]]) ){
            # peak.summit[idx.query[j]] = Bumps[[i]]$summit[thissubject][1] ### only take one of the multiple summits
            peak.summit[idx.query[j]] = paste0(Bumps[[i]]$summit[thissubject], collapse = "_") ### record summits of all bumps
          }
        }
      }
    }
    thispeak$summit = peak.summit

    #### 3. get peak count and then assign a score for each candidate peak
    if(nrow(thispeak) > 1){
      mySums = rowSums
      myMeans = rowMeans
    }else{
      mySums = sum
      myMeans = mean
    }
    bins.GR = GRanges(Rle(bins$chr), IRanges(bins$start, bins$end))
    thiscount = getPeakCounts(peaks = thispeak, allCounts = Counts, allBins = bins.GR)
    mergecount = cbind(mySums(thiscount[, seq(1, ncol(thiscount), 2)]),
                       mySums(thiscount[, seq(2, ncol(thiscount), 2)]))
    y = mergecount[, 2]
    tol = mySums(mergecount)
    score = 1 - pbinom(y, tol, prob = mean(bn.prob))

    ### add lfc on sept 25, 2019
    idx.x = seq(1, ncol(thiscount), 2)
    idx.y = seq(2, ncol(thiscount), 2)
    lfc = log(myMeans((thiscount[, idx.y]/sf[idx.y] + c0)/(thiscount[, idx.x]/sf[idx.x] + c0), na.rm = TRUE))
    if(length(lfc) > windlen){
      smooth.lfc <- mySmooth(lfc, windlen = windlen)
    }else{
      smooth.lfc = lfc
    }
    ###
    Peaks = data.frame(thispeak, score = score, lfc = smooth.lfc, thiscount)
    colnames(Peaks) = c(colnames(thispeak), "score", "lg.fc", colnames(thiscount))

    if(nrow(thiscount) >= 2){
      #### remove lowcount peaks
      idx = which(rowSums(thiscount[, seq(2, ncol(thiscount), 2)])  > lowcount)
      tmp = Peaks[idx, ]
      rownames(tmp) = as.factor(1:nrow(tmp))
      Peaks = tmp
    }
    cat("The number of peaks in this sample is: ", nrow(Peaks), sep = "\n")
  }else{
    Peaks = data.frame()
    cat("No peaks in this sample!", sep = "\n")
  }

  Peaks$width = NULL
  return(Peaks)
}

M6Apeak.MultiRep.step2 <- function(Candidates, sf, mu.cutoff){
  require(stats)
  ### Order candidate peaks from step 1 with more sophisticated statistical model
  Candidates$score = NULL   ### first remove score from step 1
  thispeak = Candidates
  idx = which(grepl("rep", colnames(thispeak)) | grepl("bam", colnames(thispeak)))
  thiscount = thispeak[, idx]

  res = M6Apeak(mat = as.matrix(thiscount), sf = sf, cutoff = mu.cutoff)
  thispeak = cbind(Candidates, res)

  #### added on March 11, 2020: rank the peak list based on their statistics
  thispeak = thispeak[order(thispeak$rScore, decreasing = TRUE), ]
  #### ended editing

  return(thispeak)
}


M6Apeak.oneRep <- function(Counts, sf = NULL, bins,
                           WhichThreshold = "fdr_lfc",
                           pval.cutoff = 1e-05,
                           fdr.cutoff = 0.05,
                           lfc.cutoff = 0.7,
                           windlen = 5,
                           lowCount = 10){
  ### peak calling for real data when there are only one replicate
  require(GenomicFeatures)
  require(GenomicRanges)
  require(stats)
  ### step 1: grasp bumps based on lfc or binomial test
  if(length(sf) == 0){
    sf = colSums(Counts)/median(colSums(Counts))
  }
  Pvals = rep(NA, nrow(Counts))
  idx = rowSums(Counts) > 0
  Pvals[idx] = 1 - pbinom(Counts[idx, 2], rowSums(Counts[idx, ]), prob = sf[2]/sum(sf))
  Pvals[!idx] = 1

  ### lfc
  c0 = mean(as.matrix(Counts), na.rm = TRUE)  ### pseudocount
  lfc = log((Counts[, 2]/sf[2] + c0)/(Counts[, 1]/sf[1] + c0))
  smooth.lfc <- mySmooth(lfc, windlen = windlen)
  x.vals = data.frame(pvals = Pvals,
                      fdr = p.adjust(Pvals, method = "fdr"),
                      lfc = lfc)

  ### find bumps based on pvals, fdr or lfc
  tmp = findBumps(chr = bins$chr, pos = bins$start, strand = bins$strand,
                  x = x.vals,
                  use = WhichThreshold,
                  pval.cutoff = pval.cutoff,
                  fdr.cutoff = fdr.cutoff,
                  lfc.cutoff = lfc.cutoff,
                  count = Counts)
  Bumps = tmp[tmp$counts > lowCount, ] ### remove low count

  if(nrow(Bumps) >= 2){
    ### step 2: binomial test based on the counts of each bumps
    peaks = Bumps[, c("chr", "start",  "end", "strand", "summit")]
    bins.GR = GRanges(Rle(bins$chr), IRanges(bins$start, bins$end))
    count = getPeakCounts(peaks = peaks, allCounts = Counts, allBins = bins.GR)
    pvals = 1 - pbinom(count[, 2], rowSums(count), prob = sf[2]/sum(sf))
    fdr = p.adjust(pvals, method = "fdr")
    peaks = cbind(peaks, count)
    peaks$pvals = pvals
    peaks$score = fdr

    ### calculate log fold change for peak regions
    thiscount = peaks[, which(grepl("bam", colnames(peaks)) | grepl("rep", colnames(peaks)) ) ]
    c0 = mean(as.matrix(thiscount), na.rm = TRUE)
    lfc = log((thiscount[, 2]/sf[2] + c0)/(thiscount[, 1]/sf[1] + c0))
    if(length(lfc) > windlen){
      smooth.lfc <- mySmooth(lfc, windlen = windlen)
    }else{
      smooth.lfc = lfc
    }
    peaks$lg.fc = smooth.lfc

    ### Added on March 11, 2020: rank peaks based on the log fc
    peaks = peaks[order(peaks$lg.fc, decreasing = TRUE), ]
    ### eneded editing

    return(peaks)
  }else{
    cat("Less than 2 peaks!", sep = "\n")
    peaks = Bumps
    return(peaks)
    }

}



