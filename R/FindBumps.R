## should consider chr, but this should be okay most of the time
findRegion <- function(chr, pos, sep=1000) {
  pos.diff <- abs(c(as.integer(0), diff(pos)))
  idx.jump <- which(pos.diff>sep)
  regions <- rbind(c(1, idx.jump), c(idx.jump-1, length(pos)))
  regions
}

findBumps <- function(chr, pos, strand, x, count,
                      use = "pval", 
                      pval.cutoff, 
                      fdr.cutoff,
                      lfc.cutoff,
                      sep = 2000, minlen=100, minCount=3, dis.merge=100,
                      scorefun = mean, sort=TRUE) {
  ### This function is used to combine significant bins to form bumps
  ### given the location and p-value of each bin
  
  if(sep < dis.merge)
    sep = dis.merge + 1
  
  if(use == "pval"){
    flag <- as.numeric(x$pvals < pval.cutoff)
    flag[is.na(flag)]=FALSE
  }else if(use == "fdr"){
    flag <- as.numeric(x$fdr < fdr.cutoff)
    flag[is.na(flag)]=FALSE
  }else if(use == "lfc"){
    flag <- as.numeric(x$lfc > lfc.cutoff)
    flag[is.na(flag)]=FALSE
  } else if(use == "pval_lfc"){
    flag <- as.numeric(x$pvals < pval.cutoff & x$lfc > lfc.cutoff)
    flag[is.na(flag)]=FALSE
  }else if(use == "fdr_lfc"){
    flag <- as.numeric(x$fdr < fdr.cutoff & x$lfc > lfc.cutoff)
    flag[is.na(flag)]=FALSE
  }
  
  ## divide the whole genome into consecutive regions, which were sequenced
  regions <- findRegion(chr, pos, sep)
  ## loop on regions: within each region, combine significant bins to form a bump. 
  ## there may be multiple bumps (candidate peaks) within one region
  initn <- 100000
  result <- data.frame(chr=rep("chr1",initn), start=rep(0, initn),
                       end=rep(0, initn), length=rep(0, initn), 
                       strand = rep(NA, initn),
                       summit = rep(0, initn), counts = rep(0, initn), 
                       score=rep(0, initn))
  levels(result[,1]) <- unique(chr)
  result.idx <- 0
  Peakstrand <- NULL
  for(i in 1:ncol(regions)) {
    idx <- regions[1,i]:regions[2, i]
    pos.region <- pos[idx]
    strand.region <- strand[idx]
    # count (only this count calculation was added by Zhenxing)
    count.region <- count[idx, ]
    #
    if(length(idx)<minCount) next
    nn <- length(idx)
    flag.region <- flag[idx]
    ## get start/end position
    startidx <- which(flag.region[-nn]==0 & flag.region[-1]==1)+1
    if(flag.region[1]==1)
      startidx <- c(1, startidx)
    if(length(startidx)==0)
      next
    endidx <- which(flag.region[-nn]==1 & flag.region[-1]==0)
    if(flag.region[nn]==1)
      endidx <- c(endidx, nn)
    
    ## remove if there are less than minCount probes
    idx.keep <- (endidx-startidx+1)>=minCount
    startidx <- startidx[idx.keep]
    endidx <- endidx[idx.keep]
    if(length(endidx)==0) next
    
    # ## merge if they are really close
    nbump <- length(startidx)
    if(nbump>1) {
      bumppos <- cbind(pos[idx][startidx], pos[idx][endidx])
      dis <- bumppos[-1,1]>(bumppos[-nbump,2]+dis.merge)
      idx.start <- which(c(1,dis)==1)
      idx.end <- which(c(dis,1)==1)
      ## merged
      startidx <- startidx[idx.start]
      endidx <- endidx[idx.end]
    }
    nbump <- length(startidx)
    ll <- pos.region[endidx] - pos.region[startidx] + 1 ### length of each bump
    tmpn <- length(ll)
    # ## make bump scores
    x.thisregion <- x[idx, ]
    scores.thisregion <- rep(0, nbump)
    counts.thisregion <- rep(0, nbump)
    summit.thisregion <- rep(0, nbump)
    strand.thisregion <- rep(0, nbump)
    for(ibump in 1:nbump){
      thisrange <- startidx[ibump]:endidx[ibump]
      scores.thisregion[ibump] <- scorefun(x.thisregion$pval[thisrange])
      counts.thisregion[ibump] <- sum(count.region[thisrange,])
      thispos <- pos.region[thisrange]
      thistrand <- strand.region[thisrange]
      strand.thisregion[ibump] <- defineStrand(thistrand)
      summit.thisregion[ibump] <- thispos[which.max(count.region[thisrange, 2])]+24
      # summit.thisregion[ibump] <- thispos[which.max(rowSums(count.region[thisrange, ]) - 
      #                                                 mean(rowSums(count.region[thisrange,])))]
    }
    
    result[result.idx+(1:tmpn),] <- data.frame(chr = as.character(chr[idx][startidx]),
                                               start = pos[idx][startidx],
                                               end = pos[idx][endidx]+ 49, 
                                               length=ll+49, 
                                               strand = as.character(strand.thisregion),
                                               summit = summit.thisregion,
                                               count = counts.thisregion,
                                               score = scores.thisregion)
    Peakstrand = c(Peakstrand, strand.thisregion)
    
    result.idx <- result.idx + tmpn
  }
  
  result <- result[1:result.idx,]
  result$strand <- Peakstrand
  ## remove really short ones
  result <- result[result[,4]>minlen,]
  ## sort according to score

  ii <- sort(result$score, decreasing = FALSE, index=TRUE)
  result <- result[ii$ix,]
  result
}



makeBumps <- function(chr, pos, x, usestat = FALSE, cutoff, sep=2000, ext=0,
                      minlen=100, minCount=3, dis.merge=100,
                      index.return=FALSE, summit=TRUE) {
  if(usestat){
    flag <- as.numeric(x<cutoff)
    flag[is.na(flag)]=FALSE 
  }else if(!usestat){
    flag <- as.numeric(x>cutoff)
    flag[is.na(flag)]=FALSE 
  }
  
  
  ## loop on regions
  initn <- 20000
  if(summit)
    bumps <- data.frame(chr=rep("chr1",initn), start=rep(0,initn),
                        end=rep(0, initn), peak=rep(0, initn), length=rep(0, initn))
  else 
    bumps <- data.frame(chr=rep("chr1",initn), start=rep(0,initn),
                        end=rep(0, initn), length=rep(0, initn))
  if(index.return)
    res.idx <- matrix(0, nrow=initn, ncol=2)
  levels(bumps[,1]) <- unique(chr)
  
  counter <- 0
  idx <- split(1:length(chr), chr)
  for(ichr in seq(along=idx)) {
    thispos=pos[idx[[ichr]]]
    tmp <- .Call("findBumps", as.integer(thispos), x[idx[[ichr]]], as.double(cutoff), as.double(sep),
                 as.double(minlen), as.integer(minCount), as.double(dis.merge), as.integer(ext),
                 as.integer(summit))
    
    if(is.null(tmp))
      next
    if(!summit) {
      tmp=matrix(tmp, byrow=TRUE, ncol=2) + 1
      nn=nrow(tmp)
      if(nn>0) {
        if(index.return) 
          res.idx[counter+1:nn,] <- matrix(idx[[ichr]][tmp], ncol=2)
        bumps[counter+1:nn,] <- data.frame(chr=names(idx)[ichr], start=thispos[tmp[,1]],end=thispos[tmp[,2]])
        counter=counter+nn
      }
    } else { ## report summit of the peaks
      tmp2=matrix(tmp[[1]], byrow=TRUE, ncol=2) + 1
      nn=nrow(tmp2)
      if(nn>0) {
        if(index.return)
          res.idx[counter+1:nn,] <- matrix(idx[[ichr]][tmp2], ncol=2)
        bumps[counter+1:nn,] <- data.frame(chr=names(idx)[ichr], start=thispos[tmp2[,1]],
                                           end=thispos[tmp2[,2]], peak=thispos[tmp[[2]]+1])
        counter=counter+nn
      }
    }
  }
  
  if(counter>0) {
    bumps <- bumps[1:counter,]
    bumps$length <- bumps$end - bumps$start + 1
    if(index.return) {
      res.idx <- res.idx[1:counter,,drop=FALSE]
      return(list(bumps=bumps, ix=res.idx))
    }
    else
      bumps
  }
  else {
    return(NULL)
  }
  
}



defineStrand <- function(strand){
  if(all(strand== "+")){
    peakStrand = "+"
  }else if (all(strand == "-")) {
    peakStrand = "-"
  } else if (sum(strand == "*")>0 | (sum(strand == "+") >0 && sum(strand == "-") >0)){
    peakStrand = "*"
  }else if(sum(strand == ".") >0){
    peakStrand = "."
  }
  return(peakStrand)
}


addStrand <- function(peak, Bininf){
  ### add strand informaiton into peaks
  bin.GR = GRanges(seqnames = Rle(Bininf$seq), 
                   ranges = IRanges(start = Bininf$position.start, end = Bininf$position.end))
  peak.GR = GRanges(seqnames = Rle(peak$chr), ranges = IRanges(peak$start, peak$end))
  strand = Bininf$strand
  hits = findOverlaps(peak.GR, bin.GR)
  query = queryHits(hits)
  subject = subjectHits(hits)
  peakStrand = NULL
  
  for( i in 1:nrow(peak)){
    idx = which(query == i)
    if(all(strand[subject[idx]] == "+")){
      peakStrand[i] = "+"
    }else if (all(strand[subject[idx]] == "-")) {
      peakStrand[i] = "-"
    } else if (sum(strand[subject[idx]]  == "*")>0 | (sum(strand[subject[idx]] == "+") >0 && sum(strand[subject[idx]] == "-") >0)){
      peakStrand[i] = "*"
    }else if(sum(strand[subject[idx]]  == ".") >0){
      peakStrand[i] = "."
    }
  }
  
  peak$strand = peakStrand
  return(peak)
}
