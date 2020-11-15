# TRES
Suppose you have a set of bam files for both input and IP from "mm9" genome, then a quick start of using TRES is:

############
IP.file = c("ip_rep1.bam", "ip_rep2.bam")
Input.file = c("input_rep1.bam", "input_rep2.bam")
BamDir = YourInputDir
OutDir = YourOutputDir
TRES_peak(IP.file = IP.file,
          Input.file = Input.file,
          genomeBuild = "mm9",
          InputDir = BamDir,
          OutputDir = OutDir,
          experiment_name = "example",
          filetype = "bam")
peaks = read.table(paste0(OutDir, "example_peaks.xls"), sep = "\t", header = TRUE)
head(peaks)
