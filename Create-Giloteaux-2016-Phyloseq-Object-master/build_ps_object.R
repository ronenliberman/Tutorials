


##########################################################################

#Converting the data provided in the "Differential abundance analysis with gneiss"
#tutorial into a phyloseq object for analysis in R


# Data provided at: https://docs.qiime2.org/2019.4/tutorials/gneiss/
# Microbiome publication 2016: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0171-4  
# qiime2R at:  https://rdrr.io/github/jbisanz/qiime2R/f/README.md


#Created by: Nicholas Ollberding (for CCHMC/UC Summer Microbiome Course Session)

#On: 6/13/19

#Data accessed: 06/13/19


#########################################################################



#Installing qiime2R: https://rdrr.io/github/jbisanz/qiime2R/f/README.md
#install.packages("remotes")
#remotes::install_github("jbisanz/qiime2R")


#Loading qiime2R
library(qiime2R); packageVersion("qiime2R")      #version: 0.99.11
library(phyloseq); packageVersion("phyloseq")    #version: 1.28.0



#Generating phyloseq object from table and taxa files
dir.create("tmp")
(ps <- qza_to_phyloseq(features = "table.qza", taxonomy = "taxa.qza", tmp = "tmp"))



#Reading in metadata
meta <- read.delim("metadata.txt", header = TRUE, sep = "\t")
meta <- data.frame(meta, row.names = 1)
(ps <- merge_phyloseq(ps, sample_data(meta)))
rm(meta)



#Generating de-novo tree
library(DECIPHER); packageVersion("DECIPHER")      #version: 2.12.0
library(phangorn); packageVersion("phangorn")      #version: 2.5.3
library(dada2); packageVersion("dada2")            #version: 1.12.1

otu <- data.frame(otu_table(ps))
otu <- t(otu)

seqs <- getSequences(otu)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) 
fit = pml(treeNJ, data=phang.align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))


(ps <- merge_phyloseq(ps, phy_tree(fitGTR$tree)))



#Adding reference sequences
refs <- DNAStringSet(seqs)
(ps <- merge_phyloseq(ps, refseq(refs)))



#Renaming taxa names to OTU1, OTU2, etc. 
taxa_names(ps) <- paste0("OTU", seq(ntaxa(ps)))



#Fixing taxonomy
library(tidyverse); packageVersion("tidyverse")            #version: 1.2.1
taxa_tab <- data.frame(tax_table(ps)[, 1])

taxa_tab <- taxa_tab %>%
  rownames_to_column(var = "OTU") %>%
  separate(Kingdom, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  mutate(Kingdom = gsub("k__", "", Kingdom),
         Phylum = gsub("p__", "", Phylum),
         Class = gsub("c__", "", Class),
         Order = gsub("o__", "", Order),
         Family = gsub("f__", "", Family),
         Genus = gsub("g__", "", Genus),
         Species = gsub("s__", "", Species))

empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

taxa_tab <- taxa_tab %>% mutate_each(list(empty_as_na)) 
taxa_tab <- data.frame(taxa_tab, row.names = 1)
taxa_tab <- as.matrix(taxa_tab)

tax_table(ps) <- NULL
(ps <- merge_phyloseq(ps, tax_table(taxa_tab)))



#Saving ps object
saveRDS(ps, "ps_giloteaux_2016.rds")








