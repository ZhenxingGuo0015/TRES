# MeRIPeak
A quick start to use this package:


bamDir = "Your bam file dierctory"

peakDir = "Your output directory"

Sample = "The folder to store your results"

IP.bam = c("ip1.bam", "ip2.bam")

Input.bam = c("input1.bam", "input2.bam")

MeRiPeak(IP.file = IP.bam,  
         Input.file = Input.bam, 
         InputDir = bamDir,         
         OutputDir = peakDir,         
         experiment_name = Sample)
