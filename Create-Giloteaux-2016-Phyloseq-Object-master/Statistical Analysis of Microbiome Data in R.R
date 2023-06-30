#This is a tutorial about statistical analyis of Microbiome data
#This tutorial follows this link
#https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
#I will also implement examples from other sources

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
#Get count of phyla
table(phyloseq::tax_table(ps)[, "Phylum"])
#Convert to relative abundance
ps_rel_abund = phyloseq::transform_sample_counts(ps, function(x){x / sum(x)})
phyloseq::otu_table(ps)[1:5, 1:5]
phyloseq::otu_table(ps_rel_abund)[1:5, 1:5]

#Plot
phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Status, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#Agglomerate to phylum-level and rename
ps_phylum <- phyloseq::tax_glom(ps, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]

#Melt and plot
phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = Status, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  theme_bw()+
  facet_wrap(~ OTU, scales = "free")

#Statistical testing 
#test for a difference in the phylum-level abundance is to conduct a multivariate test for differences in the overall composition
#between groups of samples

#based on a paper:  Hypothesis Testing and Power Calculations for Taxonomic-Based Human Microbiome Data
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052078

#It is suggested that rare taxa be pooled into a single group to improve testing.

#Subset groups
controls <- phyloseq::subset_samples(ps_phylum, Status == "Control")
cf <- phyloseq::subset_samples(ps_phylum, Status == "Chronic Fatigue")
#Output OTU tables
control_otu <- data.frame(phyloseq::otu_table(controls))
cf_otu <- data.frame(phyloseq::otu_table(cf))
#Group rare phyla
control_otu <- control_otu %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(Other = Cyanobacteria + Euryarchaeota + Tenericutes + Verrucomicrobia + Fusobacteria) %>%
  dplyr::select(-Cyanobacteria, -Euryarchaeota, -Tenericutes, -Verrucomicrobia, -Fusobacteria)
cf_otu <- cf_otu %>%
  t(.) %>%
  as.data.frame(.) %>%
  mutate(Other = Cyanobacteria + Euryarchaeota + Tenericutes + Verrucomicrobia + Fusobacteria) %>%
  dplyr::select(-Cyanobacteria, -Euryarchaeota, -Tenericutes, -Verrucomicrobia, -Fusobacteria)
#HMP test
group_data <- list(control_otu, cf_otu)
(xdc <- HMP::Xdc.sevsample(group_data))  

# This results show the the HMP test fails to reject the null hypothesis that there is no difference in the dist. between the two groups. 


#Hierarchical clustering - here using Bray-Curtis distances
#Extract OTU table and compute BC
ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
head(ps_rel_otu)
ps_rel_otu <- t(ps_rel_otu) #transposing- Each row needs to be the sample and each column is the count per OTU/ASV
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
as.matrix(bc_dist)[1:5, 1:5]

#Save as dendrogram
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#Provide color codes
meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
colorCode <- c(Control = "coral", `Chronic Fatigue` = "purple")
labels_colors(ward) <- colorCode[meta$Status][order.dendrogram(ward)]
#Plot
plot(ward)

## now i need to try a heat map

#Alpha Diversity - the diversity in a single ecosystem or sample. 
#The simplest measure is richness, the number of species (or OTUs) observed in the sample. 

ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(ps),
                         "observed" = phyloseq::estimate_richness(ps, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")

#Subsample reads
(ps_rare <- phyloseq::rarefy_even_depth(ps, rngseed = 123, replace = FALSE)) #using this function to equaly subsample from all the rows (samples
head(phyloseq::sample_sums(ps_rare))

adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(ps_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps_rare, measures = "Shannon"),
  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps_rare)))), tree = phyloseq::phy_tree(ps_rare),include.root = F)[, 1],
  "Status" = phyloseq::sample_data(ps_rare)$Status
)

head(adiv)

#Plot adiv measures
adiv %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "PD")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "PD"))) %>%
  ggplot(aes(x = Status, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Status), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme_bw()+
  theme(legend.position="none")

#Summarize
adiv %>%
  group_by(Status) %>%
  dplyr::summarise(median_observed = median(Observed),
                   median_shannon = median(Shannon),
                   median_pd = median(PD))

#Wilcoxon test of location
wilcox.test(Observed ~ Status, data = adiv, exact = FALSE, conf.int = TRUE)

# we can do it on all the 3 parameters that were checked. But, the test are also visivle on the plot
# We find evidence to support the null hypothesis that there is no significant difference in location between the groups.


#Estimate the richness using breakaway - this is a method to reduce issue in depth
#https://github.com/adw96/breakaway
#Obtain breakaway estimates

ba_adiv <- breakaway(ps)
ba_adiv[1]
#Plot estimates
plot(ba_adiv, ps, color = "Status")     
#Examine models
summary(ba_adiv) %>%
  add_column("SampleNames" = ps %>% otu_table %>% sample_names) %>% 
  add_column("Observed" = phyloseq::estimate_richness(ps, measures = "Observed")$Observed) #I added the original value here
#Test for group differnce
bt <- breakaway::betta(summary(ba_adiv)$estimate,
                       summary(ba_adiv)$error,
                       make_design_matrix(ps, "Status"))
bt$table  


## Beta-diversity
#Beta-diversity provides a measure of similarity, or dissimilarity, 
#of one microbial composition to another. 

# There is a lot of Information about *compositional data analysis approach* in the tutorial. read it to understand more

#First we transform the data
#CLR transform
(ps_clr <- microbiome::transform(ps, "clr"))
phyloseq::otu_table(ps)[1:5, 1:5]
phyloseq::otu_table(ps_clr)[1:5, 1:5] # see the difference, the data is no longer counts but transfromed to the geometric mean of all taxa on the log scale.

#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig)) # calculates the proportions of each PC and returns a vector

#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)


phyloseq::plot_ordination(ps, ord_clr, type="samples", color="Status") + 
  geom_point(size = 2) +
  stat_ellipse(aes(group = Status), linetype = 2)

#statistical examination of diffrences between control and the fatigue using PERMANOVA & ADONIS

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$Status)

#Dispersion test and plot - type of assumption testing for PERMANOVA
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$Status)
dispr

plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr)

#Calculating UniFrac metric 
#Generate distances
ord_unifrac <- ordinate(ps_rare, method = "PCoA", distance = "wunifrac")
ord_unifrac_un <- ordinate(ps_rare, method = "PCoA", distance = "unifrac")   
#Plot ordinations
a <- plot_ordination(ps_rare, ord_unifrac, color = "Status") + geom_point(size = 2)+
  geom_point(size = 2) +
  theme_bw()
b <- plot_ordination(ps_rare, ord_unifrac_un, color = "Status") + geom_point(size = 2)+
  geom_point(size = 2) +
  theme_bw()
cowplot::plot_grid(a, b, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))

