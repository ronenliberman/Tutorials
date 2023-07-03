
#'This tutorial follows this link
#'https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
#'I will also implement examples from other sources
#'
#+ echo=FALSE

cran_packages <- c("tidyverse", "cowplot", "picante", "vegan", "HMP", "dendextend", "rms", "devtools")
.bioc_packages <- c("phyloseq", "DESeq2", "microbiome", "metagenomeSeq", "ALDEx2")
.bioc_packages2 <-c("ALDEx2","metagenomeSeq")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(.bioc_packages2)
 devtools::install_github("adw96/breakaway")
# Installing the files in the repository
devtools::install_github(repo = "malucalle/selbal")

install.packages("ALDEx2","metagenomeSeq","HMP","rms")

#' Loading the library
library("knitr")
#spin("~/R/Tutorials/Tutorials/Create-Giloteaux-2016-Phyloseq-Object-master/Statistical Analysis of Microbiome Data in R.R", knit = FALSE)
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

# There is a lot of Information about *compositional data analysis approach* in the
#tutorial. read it to understand more

#First we transform the data
#CLR transform
(ps_clr <- microbiome::transform(ps, "clr")) #centered log-ratio (CLR) transformation
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

##Differential abundance testing
#+The goal of differential abundance testing is 
#+to identify specific taxa associated with clinical metadata variables of interest.
#+
#+will also use *nested data* frames as advocated by Hadley Wickham 
#+to keep the data and test results together in a single data.frame

#Generate data.frame with OTUs and metadata
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps_clr))))
ps_wilcox$Status <- phyloseq::sample_data(ps_clr)$Status 
#Define functions to pass to map
wilcox_model <- function(df){
  wilcox.test(abund ~ Status, data = df) 
}
wilcox_pval <- function(df){
  wilcox.test(abund ~ Status, data = df)$p.value
}
#Create nested data frames by OTU and loop over each using map
wilcox_results <- ps_wilcox %>%
  gather(key = OTU, value = abund, -Status) %>%
  group_by(OTU) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))                       
#Show results
head(wilcox_results)
head(wilcox_results$data[[1]])
wilcox_results$wilcox_test[[1]]
wilcox_results$p_value[[1]]

#Un-Nesting
wilcox_results <- wilcox_results %>% 
  dplyr::select(OTU, p_value) %>%
  unnest(cols = c(p_value))
head(wilcox_results)

#Adding taxonomic labels
taxa_info <- data.frame(tax_table(ps_clr))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
view(taxa_info)

#Computing FDR corrected p-values
wilcox_results <- wilcox_results %>%
  full_join(taxa_info) %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(OTU, p_value, BH_FDR, everything())
## Joining, by = "OTU"
#Printing results
print.data.frame(wilcox_results) 

#Another approach - ANOVA-like differential expression (ALDEx2)  
#CoDA method for differential abundance testing. ALDEx2 can be run via a single command.
#But, the are actually many complicated steps in the backround, check the Tutorial for more explanations. 
#There are several steps that are occurring in the background. the steps include: check it: ?aldex 
#+Generate a large number (here n=128) of posterior probabilities for the observance of each taxon (i.e. output many data.frames where the counts have been converted to proportions). This is done by Monte-Carlo sampling from a Dirichlet distribution with a small non-zero prior to deal with zeros. The total read count therefore only contributes to the precision of the proportions.
#+Apply the centered log-ratio transformation to each instance.
#+Apply the Wilcoxon test to each taxon for each simulated instance.
#+Estimate the effect size as the difference between conditions divided by the maximum difference within conditions averaging over all instances. Scaling the between group difference by the maximum within group difference gives us a standardized effect size measure.
#+Obtain the expected p-values for each taxon by averaging over all instances.
#+Apply the BH-FDR correction to control the false positive rate.

#Run ALDEx2
aldex2_da <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps)),
                           phyloseq::sample_data(ps)$Status,
                           test="t", effect = TRUE, denom="iqlr")

ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)


#Clean up presentation
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  #filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)
#head(sig_aldex2)
sig_aldex2 <- left_join(sig_aldex2, taxa_info)
  
sig_aldex2
#We see that again that several Clostridiales organisms are identified as differentially abundant. Consistent with the results of running the Wilcoxon test outside of ALDEx2, we see that OTU48, OTU38, OTU44, and OTU8 are listed as differentially abundant. The others do not reach the FDR cut-off used here; although, they likely have “largish” effect sizes. Try and see if you can obtain these values. The reason for the discrepancy is hard to discern, 
#+but may be related to differences in the use of the CLR basis (geometric mean of all taxa versus the IQLR)
#+ and/or the use of the Bayesian resampling with a non-zero prior
#+ 
#+ 
##Prediction model for microbiome data
#+One problem we face when building a predictive model from metagenomic data is
#+ that we often have more features (taxon) than we have samples.

#First we will create a data.frame that contains the Status and the first 3 PCs from the centered-log ratio
#transformed abundance table we generated before. We will then plot the unconditional association for each PC with the outcome of CF versus control.

