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

#this si taking a lot of CPU
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
                    c(filtFs[[i]], filtRs[[i]]),
                    trimLeft=10, truncLen=c(245, 160),
                    maxN=0, maxEE=2, truncQ=2,
                    compress=TRUE)
}
#Generating an error model of our data
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names

#Alternative - Generating an error model of our data
err_forward_reads <- learnErrors(filtFs,multithread=TRUE)
err_reverse_reads <- learnErrors(filtRs, multithread=TRUE)

#if the black line is in line with the black dots than i am ok

plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)

#Inferring ASVs
dada_forward <- dada(derepFs, err=err_forward_reads, pool=TRUE)
dada_reverse <- dada(derepRs, err=err_reverse_reads, pool=TRUE)
rm(dada_forward)

#Merging 
mergers <- mergePairs(dada_forward, derepFs, dada_reverse, derepRs)

#Construct sequence table and remove chimeras
seqtab.all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
seqtab <- removeBimeraDenovo(seqtab.all)

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)
sum(seqtab.nochim)/sum(seqtab)

#Assign taxonomy
## loading reference taxonomy object >> doing it manually 
library(DECIPHER)
packageVersion("DECIPHER")

## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

## and classifying
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)
View(tax_info)

#Extracting the standard goods from DADA2
#The typical standard outputs from amplicon processing are a fasta file, a count table, and a taxonomy table. So here’s one way we can generate those files from your DADA2 objects in R:
  # giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "New folder/ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "New folder/ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "New folder/ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

#Taxonomic filtering
rank_names(asv_tax)
table(tax_table(asv_tax)[, "phylum"], exclude = NULL) #for example by phylum , here we see that 16 were found to be NA's 

##### STUCK HERE #####

table(is.na(tax_table(asv_tax)$phylum))
table(tax_table(asv_tax)$phylum %in% c("", "uncharacterized"))

asv_tax0 <- subset_taxa(asv_tax, !is.na(phylum) & !phylum %in% c("", "uncharacterized")) #removed the unidentified ones


#Making Phylogenetic tree 
seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

mimarks_path <- "data/MIMARKS_Data_combined.csv"
samdf <- read.csv(mimarks_path, header=TRUE)
samdf$SampleID <- paste0(gsub("00", "", samdf$host_subject_id), "D", samdf$age-21)
samdf <- samdf[!duplicated(samdf$SampleID),] # Remove dupicate entries for reverse reads
rownames(seqtab) <- gsub("124", "125", rownames(seqtab)) # Fixing an odd discrepancy
all(rownames(seqtab) %in% samdf$SampleID) # TRUE

## [1] TRUE

rownames(samdf) <- samdf$SampleID
keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment",
               "host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass",
               "diet", "family_relationship", "genotype", "SampleID") 
samdf <- samdf[rownames(seqtab), keep.cols]
