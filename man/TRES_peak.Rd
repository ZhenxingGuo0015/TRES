\name{TRES_peak}
\alias{TRES_peak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An Empirical Bayesian Method for MeRIP-seq Peak Detection.
}
\description{
This function implements a two-step procedure to conduct peak calling for MeRIP-seq data.
In the first step, it loosely indentify candidate peak regions using an ad hoc procedure. In the second step, it ranks the candidate regions with more rigorous statistical modeling based on an empirical Bayesian hierarchical model.
}
\usage{
TRES_peak(IP.file, Input.file, genomeBuild, binsize,
sf0, WhichThreshold, pval.cutoff0,
fdr.cutoff0, lfc.cutoff0, mu.cutoff0,
addDEseq2, InputDir, OutputDir, experiment_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IP.file}{
  A vector of characters containing the name of bam files for IP samples.
}
\item{Input.file}{
  A vector of characters containing the name of bam files for input control samples.
}
\item{genomeBuild}{
  A character to specify the reference genome of the sequencing data, which takes value in ("mm9", "mm10", "hg18", "hg19", "dm3", "rn4").
}

\item{binsize}{
  A numerical value to specify the size of window to bin the genome and get bin-level read counts. Default value is 50.}

  \item{sf0}{
  Numerical vectors to specify the size factors of each sample. Default is NULL.
}

\item{WhichThreshold}{
  A character to specify the name of threshod to screen candidate regions in the first step. It takes among "pval", "fdr", "lfc", "pval_lfc" and "fdr_lfc".
  "pval": screen canidates only based on P-values from Binomial tests;
  "fdr": screen canidates only based on FDR;
  "lfc": screen canidates only based on log fold changes between normalized IP and normalized input read counts;
  "pval_lfc": screen canidates based on both P-values and log fold changes;
  "fdr_lfc": screen canidates based on both FDR and log fold changes.
  Default value of ``WhichThreshold" is "fdr".
}
  \item{pval.cutoff0}{
  A numerical value to specify the cutoff of p-values if identifying the candidate regions based on P-values. Default is 1e-5.
}
  \item{fdr.cutoff0}{
  A numerical value to specify the cutoff of fdr if identifying the candidate regions based on FDR. Default is 1e-5.
}
 \item{lfc.cutoff0}{
  A numerical value to specify the cutoff of log foldchanges between normalized IP and input read counts if identifying the candidate regions based on log fold changes. Default is 0.7.
}

\item{mu.cutoff0}{
  A numerical value to specify the background methylation level. Default is 0.5.
}

\item{addDEseq2}{
A logical value to indicate whether or not include DESeq2 in the second step of ranking candidate regions. Default is FALSE.
}

\item{InputDir}{
 A character to specify the input directory of all bam files.
}

\item{OutputDir}{
 A character to specify the output directory of all results.
}

\item{experiment_name}{
 A character to specify the name of results folder.
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output of this function contains two results. One is saved as ".rda", which contains all the bins and theirs counts. The other is one ".xls" file, which contains a list of peaks called under the two-step procedure.
The columns of the peak excel files are:
\item{chr}{Chromosome number of that peak.}
\item{start}{The start of genomic position of that peak.}
\item{end}{The end of genomic position of that peak.}
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
Zhenxing Guo
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
bamDir = "Your bam file dierctory"
peakDir = "Your output directory"
Sample = "The folder to store your results"

IP.bam = c("ip1.bam", "ip2.bam")
Input.bam = c("input1.bam", "input2.bam")
TRES_peak(IP.file = IP.bam, Input.file = Input.bam,
          InputDir = bamDir,
          OutputDir = peakDir,
          experiment_name = Sample)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
