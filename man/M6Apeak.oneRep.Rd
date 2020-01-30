\name{M6Apeak.oneRep}
\alias{M6Apeak.oneRep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
m6A peak calling.
}
\description{
This function conducts peak calling for in data where there is only one biological replicate.
}
\usage{
M6Apeak.oneRep(Counts, sf, bins, WhichThreshold = "fdr",
pval.cutoff = 1e-05, fdr.cutoff = 1e-05,
lfc.cutoff = 0.7, windlen = 5, lowCount = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Counts}{
  A two-column data matrix containing bin-level read counts for both IP and input samples.
}
  \item{sf}{
 A numerical vector containg size factors of both IP and input samples. If can be provided by the user, or automatically estimated using "Counts".
}
  \item{bins}{
  A dataframe containing the genomic locations (chr, start, end, strand) of each bin.
}
  \item{WhichThreshold}{
   A character specifying a threshold for significant bins in bump finding using an ad hoc algorithm. There are three options: "pval" (only use p-values), "fdr" (only use FDR), "lfc" (only use log fold change), "pval_lfc" (use both p-values and log fold changes) and "fdr_lfc" (use FDR and log fold changes). Default is "fdr".
}
  \item{pval.cutoff}{
A constant indicating the cutoff for p-value. Default is 1e-05.
}
  \item{fdr.cutoff}{
A constant indicating the cutoff for FDR. Default is 1e-05.
}
  \item{lfc.cutoff}{
A constant indicating the cutoff for log fold change. Default is 0.7 for fold change of 2.
}
  \item{windlen}{
An integer specifying the length of consecutive bins used in simple moving average smooth of log fold change.
}
  \item{lowCount}{
  An iteger to filter out m6A regions with lower read counts. Default is 50.
}
}
\details{
When there is only one replicate, it assigns a p-value based on the binomial test, calculates a log fold change for each bin. Then it combines significant bins to form m6A regions using an ad hoc algorithm. This is very similiar to the first step of peak calling when there are more than two replicates.
}
\value{
It returns an excel containing the information for each peak:
\item{chr }{Chromosome number of that peak.}
\item{start}{The start of genomic position of that peak.}
\item{end}{The end of genomic position of that peak.}
\item{strand}{The strand of genomic position of that peak.}
\item{summit}{The summit position of that peak.}
\item{pvals}{The p-value for that peak calculated based on the binomial test.}
\item{score}{Adjusted p-values using FDR.}
\item{lg.fc}{The log foldchange between normalized IP and normalized input read counts.}
Note, there are additional columns whose name involve the character ".bam". These columns contain the read counts in ip and input samples.
}

\examples{
## A toy example
data("Basal_binlevel")
M6Apeak.oneRep(Counts = allbincount[, 1:2], sf = sf0[1:2], bins = allbins)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
