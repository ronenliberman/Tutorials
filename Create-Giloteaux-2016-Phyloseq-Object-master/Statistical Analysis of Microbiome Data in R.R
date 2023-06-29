install.packages("packrat")
packrat::init("~/R/Tutorials/Tutorials")

cran_packages <- c("tidyverse", "cowplot", "picante", "vegan", "HMP", "dendextend", "rms", "devtools")
.bioc_packages <- c("phyloseq", "DESeq2", "microbiome", "metagenomeSeq", "ALDEx2")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(.bioc_packages, version = "3.9")
devtools::install_github("adw96/breakaway")
# Installing the files in the repository
devtools::install_github(repo = "malucalle/selbal")

install.packages("ALDEx2","metagenomeSeq","HMP","rms")
# Loading the library
library("selbal")
library(tidyverse); packageVersion("tidyverse")                 
## [1] '2.0.0'
library(phyloseq); packageVersion("phyloseq")                    
## [1] '1.38.0'
library(DESeq2); packageVersion("DESeq2")                        
## [1] '1.34.0'
library(microbiome); packageVersion("microbiome")               
## [1] '1.16.0'
library(vegan); packageVersion("vegan")                          
## [1] '2.6.4'
library(picante); packageVersion("picante")                       
## [1] '1.8.2'
library(ALDEx2); packageVersion("ALDEx2")                        
## [1] '1.21.1'
library(metagenomeSeq); packageVersion("metagenomeSeq")          
## [1] '1.30.0'
library(HMP); packageVersion("HMP")                              
## [1] '2.0.1'
library(dendextend); packageVersion("dendextend")                
## [1] '1.17.0'
library(selbal); packageVersion("selbal")                          
## [1] '0.1.0'
library(rms); packageVersion("rms")
## [1] '6.7.0'
library(breakaway); packageVersion("breakaway")                  
## [1] '4.8.4'

#Reading in the Giloteaux data

#loading data of this project - ps project

(ps <- readRDS("Create-Giloteaux-2016-Phyloseq-Object-master/ps_giloteaux_2016.rds"))

#Sort samples on total read count, remove <5k reads, remove any OTUs seen in only those samples
sort(phyloseq::sample_sums(ps))

(ps <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 5000)) # removed 3 samples

(ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)) 

#Assign new sample metadata field
phyloseq::sample_data(ps)$Status <- ifelse(phyloseq::sample_data(ps)$Subject == "Patient", "Chronic Fatigue", "Control")
phyloseq::sample_data(ps)$Status <- factor(phyloseq::sample_data(ps)$Status, levels = c("Control", "Chronic Fatigue"))
ps %>% 
  sample_data %>%
  dplyr::count(Status)

#Visualizing relative abundance


