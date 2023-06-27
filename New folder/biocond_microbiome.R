library("knitr")
library("BiocStyle")
opts_chunk$set(cache = FALSE,fig.path="dadafigure/")
read_chunk(file.path("src", "bioinformatics.R"))

.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

set.seed(100)

setwd("C:/Users/Ronen/Documents/R/Tutorials/Tutorials")
miseq_path <- file.path("data", "MiSeq_SOP")
filt_path <- file.path("data", "filtered")

#get a list of my file names and lists of reverse and forward reads/sequences
fns <- sort(list.files(miseq_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

#quality control before trimming . Here, the forward reads maintain high quality throughout, while the quality of the reverse reads
#drops significantly at about position 160. Therefore, we choose to truncate the forward reads at position 245,
#and the reverse reads at position 160. We also choose to trim the first 10 nucleotides of each read based on empirical
#observations across many Illumina datasets that these base positions are particularly likely to contain pathological errors.

ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }      

# trimming
#We combine these trimming parameters with standard filtering parameters, the most important being the enforcement of a maximum of 2 expected
#errors per-read.
#Trimming and filtering is performed on paired reads jointly, i.e. both reads must pass the filter for the pair to pass.

if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
                    c(filtFs[[i]], filtRs[[i]]),
                    trimLeft=10, truncLen=c(245, 160),
                    maxN=0, maxEE=2, truncQ=2,
                    compress=TRUE)
}

