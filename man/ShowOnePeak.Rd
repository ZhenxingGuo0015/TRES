\name{ShowOnePeak}
\alias{ShowOnePeak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visulization of a single peak along the genome.
}
\description{
This function plots the estimated methylation level (as bars) of each bin within a peak for each replicate, and the corresponding normalized input read depth (grey curve).
}
\usage{
ShowOnePeak(onePeak, allBins, binCounts, ext = 500, ylim = c(0, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{onePeak}{
A one-row dataframe containing the genomic position of a single peak: chr, start, end, strand.
}
  \item{allBins}{
  A dataframe contaning the genomic position of all bins used to call peaks: chr, start, end, strand.
}
  \item{binCounts}{
  A dataframe contaning the read counts of all bins for each replicate. The sample order is: input1, ip1, input2, ip2, ...
}
  \item{ext}{
  An integer indicating the length of base paris to extend the peak on both sides: (start - ext, end + ext). Default is 500.
}
  \item{ylim}{
  The range of y-axis to plot. Default is c(0, 1)
}
}

\value{
It only generates a plot. No specific output.
}



%% ~Make other sections like Warning with \section{Warning}{....} ~

\seealso{
ShowOneDMR from "DSS" package.
}
\examples{
library(TRES)
library(datasetTRES)
IP.file = c("cb_ip_rep1_chr19.bam", "cb_ip_rep2_chr19.bam")
Input.file = c("cb_input_rep1_chr19.bam", "cb_input_rep2_chr19.bam")
BamDir = file.path(system.file(package = "datasetTRES"), "extdata/")
## need to specify your own output dir if you want to run these codes.
OutDir = "/Users/zhenxingguo/Dropbox/m6A/Zhenxing/analysis/working/package/results"

TRES_peak(IP.file = IP.file,
          Input.file = Input.file,
          genomeBuild = "mm9",
          InputDir = BamDir,
          OutputDir = OutDir,
          experiment_name = "examplebyBam",
          filetype = "bam")
peaks = read.table(paste0(OutDir, "/", "examplebyBam_peaks.xls"), sep = "\t", header = TRUE)
load(paste0(OutDir, "/", "examplebyBam.rda"))
allBins = as.data.frame(bins$bins)
colnames(allBins)[1] = "chr"
allBins$strand = binStrand

for (i in 1:4) {
  ShowOnePeak(onePeak = peaks[i,], allBins = allBins, binCounts = allCounts)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
