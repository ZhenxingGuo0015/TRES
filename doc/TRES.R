## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("devtools") # if you have not installed "devtools" package
#  library(devtools)
#  install_github("https://github.com/ZhenxingGuo0015/TRES")

## ---- eval=FALSE---------------------------------------------------------
#  IP.file = c("ip1.bam", "ip2.bam")
#  Input.file = c("input1.bam", "input2.bam")
#  MyInputDir = "/Users/zhenxingguo/Documents/research/m6a/Me/working/RealdataAnalysis/data/Mouse/bam/"
#  MyOutDir = "/Users/zhenxingguo/Documents/research/m6a/Me/working/RealdataAnalysis/Results/Mouse/PeakList/MyPeak/"
#  TRES_peak(IP.file = IP.file,
#           Input.file = Input.file,
#           genomeBuild = "mm9",
#           InputDir = MyInputDir,
#           OutputDir = MyOutDir,
#           experiment_name = "example")
#  
#  peaks = read.table(paste0(MyOutDir, "example_peaks.xls"), sep = "\t", header = TRUE)
#  head(peaks)

## ---- eval=FALSE---------------------------------------------------------
#  data("Basal_binlevel")
#  ### The first 14 columns are from basal samples of mouse cortex
#  sf0 = colSums(allbincount)/median(colSums(allbincount))
#  peak_step1 = M6Apeak.MultiRep.step1(Counts = allbincount[, 1:14],
#                                      bins = allbins,
#                                      sf = sf0[1:14],
#                                      WhichThreshold = "lfc",
#                                      lfc.cutoff = 0.5)
#  Peaks = M6Apeak.MultiRep.step2(Candidates = peak_step1,
#                                 sf = sf0[1:14])
#  head(Peaks)

## ---- eval = FALSE-------------------------------------------------------
#  # A toy example
#  data("Basal_binlevel")
#  peaks = M6Apeak.oneRep(Counts = allbincount[, 1:2], sf = sf0[1:2], bins = allbins)
#  head(peaks)

## ---- eval=FALSE---------------------------------------------------------
#  data("Basal_regionlevel")
#  data("Basal_binlevel") ###  sf0 estimated from bin-level count
#  peaks = M6Apeak.MultiRep.step2(Candidates = Basal, sf = sf0, addDEseq2 = TRUE)
#  head(peaks)

