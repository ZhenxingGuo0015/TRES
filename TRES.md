TRES: Detecting m6A methylation regions from Methylated RNA Immunoprecipitation Sequencing
================

## Introduction
The post-transcriptional epigenetic modiﬁcation on mRNA is an emerging ﬁeld to study the gene regulatory mechanism and their association with diseases. Recently developed high-throughput sequencing technology named Methylated RNA Immunoprecipitation Sequencing (MeRIP-seq) enables one to proﬁle mRNA epigenetic modiﬁcation transcriptome-wide. A basic task in the analysis of MeRIP-seq data is to identify transcriptome-wide m6A regions (namely "peak calling"). The package TRES provides methods for peak calling of MeRIP-seq data, based on an empirical Bayesian hierarchical model. The method accounts for various sources of variations in the data through rigorous modeling, and achieves shrinkage estimation by borrowing informations from transcriptome-wide data to stabilize the parameter estimation. This vignette explains the use of the package by introducing typical workflows. TRES package version: 0.1.0.

## Running TRES for peak calling
### 1. Installation
From GitHub: 
```{r, eval = FALSE}
install.packages("devtools") # if you have not installed "devtools" package
library(devtools)
install_github("https://github.com/ZhenxingGuo0015/TRES")
```
There are one main function and multiple subfunctions in TRES for different usages.

## 2. Main function of peak calling: TRES_peak
It is very convenient to conduct peak calling using this function. It starts from BAM files, so all you need are BAM files of IP and input samples output from sequence alignment tools like Bowtie2. It requires you know (of course) the reference genome ( like "mm9", "mm10", "hg18", "hg19",...) of your sequencing data. "InputDir" is the directory where you store the BAM files for both IP ("ip1.bam", "ip2.bam", ...) and input ("input1.bam", "input2.bam", ...) samples. "OutDir" and "experiment_name" are your output directory and the name you want for the peak excel file respectively.

A quick example with setup of required parameters, which uses the example dataset (only chr19 from the cerebellum sample of **Young Mouse data**) saved in data package **datasetTRES** (please install it first: install_github("https://github.com/ZhenxingGuo0015/datasetTRES").):
```{r, eval= FALSE}
library(TRES)
## Use data in the form of GRanges in package "datasetTRES"
require(datasetTRES)
IP.file = c("cb_6wk_ip_rep1_chr19", "cb_6wk_ip_rep2_chr19")
Input.file = c("cb_6wk_input_rep1_chr19", "cb_6wk_input_rep2_chr19")
OutDir = YourOutputDir
TRES_peak(IP.file = IP.file,
          Input.file = Input.file,
          genomeBuild = "mm9",
          InputDir = NULL,
          OutputDir = OutDir,
          experiment_name = "example",
          filetype = "GRanges")
peaks = read.table(paste0(OutDir, "example_peaks.xls"), sep = "\t", header = TRUE)
head(peaks)

```

```{r, eval= FALSE}
library(TRES)
## or directly take bam files in "datasetTRES"
require(datasetTRES)
IP.file = c("cb_ip_rep1_chr19.bam", "cb_ip_rep2_chr19.bam")
Input.file = c("cb_input_rep1_chr19.bam", "cb_input_rep2_chr19.bam")
BamDir = file.path(system.file(package = "datasetTRES"), "extdata/")
OutDir = YourOutputDir
TRES_peak(IP.file = IP.file,
          Input.file = Input.file,
          genomeBuild = "mm9",
          InputDir = BamDir,
          OutputDir = OutDir,
          experiment_name = "examplebyBam",
          filetype = "bam")
peaks = read.table(paste0(OutDir, "/", "examplebyBam_peaks.xls"), sep = "\t", header = TRUE)
head(peaks)
```

```{r, eval=FALSE}
## or, directly take bam files from your path
IP.file = c("ip_rep1.bam", "ip_rep2.bam")
Input.file = c("input_rep1.bam", "input_rep2.bam")
BamDir = YourInputDir
OutDir = YourOutputDir
TRES_peak(IP.file = IP.file,
          Input.file = Input.file,
          genomeBuild = YourReferenceGenome,
          InputDir = BamDir,
          OutputDir = OutDir,
          experiment_name = "example",
          filetype = "bam")
peaks = read.table(paste0(OutDir, "/", "example_peaks.xls"), sep = "\t", header = TRUE)
head(peaks)

```


### TRES_peak Parameters
**IP.file**: A vector of characters containing the name of bam files for IP samples.

**Input.file**: A vector of characters containing the name of bam files for input control samples.

**genomeBuild**: A character to specify the reference genome of the sequencing data, which takes value in ("bosTau8", "bosTau9", "ce6", "ce11","canFam3", "dm3", "dm6", "danRer10", "danRer11", "galGal4", "galGal5", "galGal6", "hg17", "hg18", "hg19", "hg38", "mm8", "mm9", "mm10", "panTro5", "panTro6", "rn4", "rn5", "rn6", "rheMac3", "rheMac8", "rheMac10",
 "sacCer2", "sacCer3",  "susScr3", "susScr11").

**Path_To_AnnoSqlite**: A character to specify the path to a ".sqlite" file used for genome annotation. If it's NULL, TRES will automatically load the annotation package starting with "TxDb." that corresponds to the provided genome. Default is NULL.

**binsize**: A numerical value to specify the size of window to bin the genome and get bin-level read counts. Default value is 50.

**sf0**: Numerical vectors to specify the size factors of each sample. Default is NULL.

**WhichThreshold**: A character to specify the name of threshod to select significant bins in the first step and significant peaks regions in the second step. It takes among "pval", "fdr", "lfc", "pval_lfc" and "fdr_lfc".
  "pval": The inference is only based on P-values;
  "fdr": The inference is only based on FDR;
  "lfc": The inference is only based on log fold changes between normalized IP and normalized input read counts;
  "pval_lfc": The inference is based on both P-values and log fold changes;
  "fdr_lfc": The inference is based on both FDR and log fold changes.
  Default value of ``WhichThreshold" is "fdr_lfc".
  
**pval.cutoff0**: A numerical value to specify the cutoff of p-values. Default is 1e-5.
  
**fdr.cutoff0**: A numerical value to specify the cutoff of fdr. Default is 0.05.
  
**lfc.cutoff0**: A numerical value to specify the cutoff of log foldchanges between normalized IP and input read counts. Default is 0.7.

**lowcount**: An iteger to filter out regions with total input counts < ``lowcount". Default is 30.

**InputDir**: A character to specify the input directory of all bam files.

**OutputDir**: A character to specify the output directory of all results. Default is NULL, which will output results under current working directory.

**experiment_name**:
 A character to specify the name of results folder.
 
**filetype**:
A character to specify the format of input data. Possible choices are: "bam", "bed" and "GRanges". Note, "GRanges" only works for the example data saved in the data package ''datasetTRES''. Default is "bam".

**IncludeIntron**: A logical value indicating whether to include (TRUE) intronic regions or not (False). Default is FALSE.


### TRES_peak Outputs
The output of this function contains two set of results. One is saved as ".rda", which contains all the genomic bins from the reference genome and theirs read counts. The other is one ".xls" file, which contains a list of called peaks. The columns of the peak excel files are:

**chr**: Chromosome number.

**start**: The start of genomic position of each peak.

**end**: The end of genomic position of each peak.

**strand**: The strand of genomic position of each peak.

**summit**: The summit position of each peak.

**lg.fc**: The log foldchange between normalized IP and normalized input read counts.

**mu**: The methylation level of each peak if there are more than one replicate.

**mu.var**: The estimated variance of methylation level for each peak, if there are more than one replicate.

**stats**: The Wald test statistics of each peak, if there are more than one replicate.

**shrkPhi**: The shrinkage estimation of dispersion for mehtylation levels of each peak, if there are more than one replicate.

**shrkTheta**: The shrinkage estimation for scale parameter theta in the gamma distribution, if there are more than one replicate.

**pvals**: P-value calculated based on the Wald-test.

**p.adj**: Adjusted p-values using FDR.

**rSocre**: A score used to ranke each site. The higher the score, the higher the rank.

Note, there are additional columns whose name involves the character ".bam". These columns contain the peak-level read counts in respective samples.

## 3. Usage of different sub-functions in TRES.
### 3.1 Call peaks starting from bin-level count matrix input 
If you have more than 1 replicates and already have bin-level counts from BAM files. You can use functions "M6Apeak.MultiRep.step1" and "M6Apeak.MultiRep.step2" togther for peak calling. "M6Apeak.MultiRep.step1" detects and combine significant bins to form candidate regions, and "M6Apeak.MultiRep.step2" identifies and ranks significant candidate regions by the proposed empirical Bayesian hirarchical negative binomial model and by DESeq2 (if required).

A quick example:
```{r, eval = FALSE, message=FALSE, warning=FALSE}
library(TRES)
data("Basal_binlevel") ### The first 14 columns are from basal samples of mouse cortex
sf0 = colSums(Basal_binlevel$bincounts)/median(colSums(Basal_binlevel$bincounts))
### first step to call candidate regions
peak_step1 = M6Apeak.MultiRep.step1(Counts = Basal_binlevel$bincounts[,
grepl("basal",colnames(Basal_binlevel$bincounts))],
                                    bins = Basal_binlevel$bins,
                                    sf = sf0[grepl("basal",names(sf0))],
                                    WhichThreshold = "lfc",
                                    lfc.cutoff = 0.5)
### estimate background methylation level
idx = which(grepl("rep", colnames(peak_step1)) | 
              grepl("bam", colnames(peak_step1)))
PeakCount = peak_step1[, idx]
bgCount = colSums(Basal_binlevel$bincounts[, grepl("basal",colnames(Basal_binlevel$bincounts))]) - colSums(PeakCount)
bg.Input = bgCount[seq(1, length(bgCount), 2)]
bg.IP = bgCount[seq(2, length(bgCount), 2)]
bg.mu = mean((bg.IP/sf0[seq(2, length(bgCount), 2)])/(bg.IP/sf0[seq(2, 
                                                                    length(bgCount), 2)] + bg.Input/sf0[seq(1, length(bgCount),                                                                                                            2)]), na.rm = TRUE)
### second step to detect and rank significant m6a regions among candidates
Peaks = M6Apeak.MultiRep.step2(Candidates = peak_step1, mu.cutoff = bg.mu,sf =sf0[grepl("basal",names(sf0))])

head(Peaks)
```

#### M6Apeak.MultiRep.step1
This function assumes you have already binned the whole reference transcriptome into equal sized bins (such as 50 basepair long) and obtained the corresponding read counts within each bin for all samples.

##### Parameters
**Counts**: A data matrix containing bin-level read counts in both IP and input samples. The number of column depends on the number of replicates, where the sample order is: input1, ip1, input2, ip2, ...

**sf**: A vector specifying the size factor of each sample, which is estimated using "Counts": colSums(Counts)/median(colSums(Counts)), or provided by the user.

**bins**: A data frame containing the genomic position of equal sized bins.

**WhichThreshold**: A character specifying a threshold for significant bins in bump finding using an ad hoc algorithm. There are three options: "pval" (only use p-values), "fdr" (only use FDR), "lfc" (only use log fold change), "pval_lfc" (use both p-values and log fold changes) and "fdr_lfc" (use FDR and log fold changes). Default is "fdr_lfc".
  
**pval.cutoff**: A constant indicating the cutoff for p-value. Default is 1e-05.

**fdr.cutoff**: A constant indicating the cutoff for FDR. Default is 0.05.

**lfc.cutoff**: A constant indicating the cutoff for log fold change. Default is 0.7 for fold change of 2.

**windlen**: An integer specifying the length of consecutive bins used in simple moving average smooth of log fold change. Default is 5.

**lowcount**: An iteger to filter out candidate regions with lower read counts. Default is 30.


##### Outputs
The output of this function is a dataframe containing candidate peak information denoted by the following columns:

**chr**: Chromosome number.

**start**: Start position.

**end**: End position.

**summit**: The summit of each peak.

**score**: The adjusted P-value of each candidate region. The P-value is derived using on the binomial test based on the region-level read counts.

**lg.fc**: The log fold changes between IP and input read counts within each m6A region.

**Others**: All the other columns corresponds to read counts of IP/input within each region from each replicate.


#### M6Apeak.MultiRep.step2
##### Parameters
**Candidates**: A dataframe which contains genomic locations and read counts of m6A regions. All detailed columns can be found from the output of function "M6Apeak.MultiRep.step1()".

**sf**: A vector of size factors for each sample, which is the same as that used in the first step of identifying candidate regions.

**mu.cutoff**:
  A constant specifying the background methylation levels. This is estimated automatically based on the first step of peak calling.

**WhichThreshold**: A character specifying a threshold for significant peaks. There are three options: "pval" (only use p-values), "fdr" (only use FDR), "lfc" (only use log fold change), "pval_lfc" (use both p-values and log fold changes) and "fdr_lfc" (use FDR and log fold changes). Default is "fdr_lfc".

**pval.cutoff**: A constant indicating the cutoff for p-value. Default is 1e-05.

**fdr.cutoff**: A constant indicating the cutoff for FDR. Default is 0.05.

**lfc.cutoff**: A constant indicating the cutoff for log fold change. Default is 0.7 for fold change of 2.


##### Outputs
The output is a dataframe. In addition to the genomic locations, read counts and log foldchanges that are output by "M6Apeak.MultiRep.step1", it also contains:

**mu**: The methylation level of that peak if there are more than one replicate.

**mu.var**: The estimated variance of for methylation level of that peak, if there are more than one replicate.

**stats**: The Wald test statistics of that peak, if there are more than one replicate.

**shrkPhi**: The shrinkage estimation of dispersion for mehtylation levels of that peak, if there are more than one replicate.

**shrkTheta**: The shrinkage estimation for scale parameter theta in the gamma distribution, if there are more than one replicate.

**pvals**: P-value calculated based on the Wald-test.

**p.adj**: Adjusted p-values using FDR.

**rSocre**: A score used to ranke each region. The higher the score, the higher the rank would be.


#### M6Apeak.oneRep
If there is only one biological replicate from your experiment and you already obtain the bin-level read counts and annotation. You can use function **M6Apeak.oneRep** to conduct peak calling, which also starts with bin-level read counts as input. The input and output of this function are pretty much similar to the input of "M6Apeak.MultiRep.step1".
```{r, eval = FALSE, message= FALSE, warning= FALSE}
# A toy example
library(TRES)
data("Basal_binlevel")
bincounts = Basal_binlevel$bincounts[, 1:2]
sf0 = Basal_binlevel$sf[1:2]
bins = Basal_binlevel$bins
peaks = M6Apeak.oneRep(Counts = bincounts, sf = sf0, bins = bins)
head(peaks)
```
##### Parameters

**Counts**: A two-column data matrix containing bin-level read counts for both IP and input samples.

**sf**: A numerical vector containing size factors of both IP and input samples. If can be provided by the user, or automatically estimated using "Counts".

**bins**: A dataframe containing the genomic locations (chr, start, end, strand) of each bin.

**WhichThreshold**: A character specifying a threshold for significant bins in bump finding using an ad hoc algorithm. There are three options: "pval" (only use p-values), "fdr" (only use FDR), "lfc" (only use log fold change), "pval_lfc" (use both p-values and log fold changes) and "fdr_lfc" (use FDR and log fold changes). Default is "fdr_lfc".

**pval.cutoff**: A constant indicating the cutoff for p-value. Default is 1e-05.

**fdr.cutoff**: A constant indicating the cutoff for FDR. Default is 0.05.

**lfc.cutoff**: A constant indicating the cutoff for log fold change. Default is 0.7 for fold change of 2.

**windlen**: An integer specifying the length of consecutive bins used in simple moving average smooth of log fold change. Default is 5.

**lowCount**: An iteger to filter out m6A regions with lower read counts. Default is 10.

##### outputs
It returns an excel containing the information for each peak:

**chr**: Chromosome number of that peak.

**start**: The start of genomic position of that peak.

**end**: The end of genomic position of that peak.

**strand**:The strand of genomic position of that peak.

**summit**: The summit position of that peak.

**pvals**: The p-value for that peak calculated based on the binomial test.

**score**: Adjusted p-values using FDR.

**lg.fc**: The log foldchange between normalized IP and normalized input read counts.

Note, there are additional columns whose name involve the character ".bam". These columns contain the read counts in ip and input samples.



### 3.2 Re-rank existing peaks
If you already have a list of peaks and the read counts but you want to re-rank them using TRES. You can run function "M6Apeak.MultiRep.step2" to do this. This may perform bad if you didn't appropriately estiamte size factors for each sample. Based on our experience, the estimation of size factor should be based on the bin-level counts across the whole transcriptome, not the  region-level counts. For background methylation level, you can use 0.5 but it would be informative if you can estimate it from your data. 

```{r, eval=FALSE, message= FALSE, warning= FALSE}
library(TRES)
data("Basal_regionlevel") ### load candidate regions
data("Basal_binlevel") ### load sf0 estimated from bin-level count
head(Basal_regionlevel)
peaks = M6Apeak.MultiRep.step2(Candidates = Basal_regionlevel,sf = Basal_binlevel$sf,mu.cutoff = 0.5,WhichThreshold = "fdr")
head(peaks)
```


### 3.3 Visualization of individual peaks
In addition to a set of sub-functions for m6A regions identification, TRES also contains one function "ShowOnePeak" to visualize individual peaks. The usage of this function is
```{r, eval=FALSE, message= FALSE, warning= FALSE}
ShowOnePeak(onePeak, allBins, binCounts, ext = 500, ylim = c(0,1))
```
In order to run this function, you need to have in hand for: 1) "onePeak": a pre-called peak saved as a dataframe, which contains genomic positions for that peak: "chr", "start", "end", "strand"; 2) "allBins": genomic positions ("chr", "start", "end", "strand") of all bins you used to call peaks; 3) "binCounts": the corresponding bin-level read counts in each replicate. This function will plot for each replicate: the methylation level of bins (blue bars) within the target peak(shade region in pink), and the normalized sequencing depth for input samples (curves in grey).
We show some example plots here:

```{r, eval=FALSE, message= FALSE, warning= FALSE}
library(datasetTRES)
library(TRES)
IP.file = c("cb_ip_rep1_chr19.bam", "cb_ip_rep2_chr19.bam")
Input.file = c("cb_input_rep1_chr19.bam", "cb_input_rep2_chr19.bam")
BamDir = file.path(system.file(package = "datasetTRES"), "extdata/")
OutDir = "/Users/zhenxingguo/Dropbox/m6A/Zhenxing/package/results"
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

# #for (i in 1:4) {
 #png(file = paste0("ExamplePeaks.png"), bg = "transparent")
   ShowOnePeak(onePeak = peaks[2,], allBins = allBins, binCounts = allCounts)
 #dev.off()
# #}
```

```{r echo=FALSE, out.height= "80%", out.width = "80%", include=TRUE}
knitr::include_graphics("ExamplePeaks.png")
```
