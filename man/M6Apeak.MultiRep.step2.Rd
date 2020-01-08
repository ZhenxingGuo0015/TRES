\name{M6Apeak.MultiRep.step2}
\alias{M6Apeak.MultiRep.step2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rank m6A candiates.
}
\description{
This is a wrap function, which calls another function to rank candidate m6A regions using Bayesian hierarchical negative binomial models.
}
\usage{
M6Apeak.MultiRep.step2(Candidates, sf, mu.cutoff, addDEseq2 = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Candidates}{
  A dataframe which contain genomic location and read counts of candidate m6A regions.
}
  \item{sf}{
  A vector of size factors for each sample, which is the same as that used in the first step of identifying candidate regions.
}
  \item{mu.cutoff}{
  A constant specifying the background methylation levels. Default is 0.5.
}
  \item{addDEseq2}{
  A logical value, which indicate whether also ranking candidate regions using DESeq2 by treat the IP and input as paired data. Default is FALSE.
}
}
\details{
This functions maninly call funciton "M6APeak" to rank candidate regions.
}
\value{
The output is a dataframe whose columns are:
The columns of the peak excel files are:
\item{chr}{Chromosome number of the peak.}
\item{start}{The start of genomic position of that peak.}
\item{end}{The end of genomic position of that peak.}
\item{width}{The genomic length in base pair of that peak.}
\item{strand}{The strand of genomic position of that peak.}
\item{summit}{The summit position of that peak.}
\item{lg.fc}{The log foldchange between normalized IP and normalized input read counts.}
\item{mu}{The methylation level of that peak if there are more than one replicate.}
\item{mu.var}{The estimated variance of for methylation level of that peak, if there are more than one replicate}
\item{stats}{The Wald test statistics of that peak, if there are more than one replicate}
\item{shrkPhi}{The shrinkage estimation of dispersion for mehtylation levels of that peak, if there are more than one replicate}
\item{shrkTheta}{The shrinkage estimation for scale parameter theta in the gamma distribution, if there are more than one replicate}
\item{pvals}{ P-value calculated based on the Wald-test.}
\item{score}{Adjusted p-values using FDR}
Note, there are additional columns whose name involve the character ".bam". These columns contain the read counts in respective samples.
}

\examples{
### A toy example
data("Basal_regionlevel")
data("Basal_binlevel") ### only use sf0
M6Apeak.MultiRep.step2(Candidates = Basal, sf = sf0, mu.cutoff = 0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
