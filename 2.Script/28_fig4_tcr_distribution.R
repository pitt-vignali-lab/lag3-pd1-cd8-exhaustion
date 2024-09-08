# ----------------------------------------------------------------------------
# Fig 4 | Distribution of clonotypes across clusters
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/4.1.3.TCR.R
# Input    : 1.Data/tcrCombined.csv
#            1.Data/TCR_defined_clusters.csv
# Output   : 3.Output/4.TCR/TCR_distribution.png
# Upstream : 27_fig4_tcr_clonotypes.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(scRepertoire)


dat <- readRDS("1.Data/datInt.RDS")
dat <- AddMetaData(dat, 
                   read.csv("1.Data/TCR_defined_clusters.csv", row.names = 1) %>% 
                     rename(newClst = dat.newClst))

DimPlot(dat, group.by = "newClst", label = T)

datNor <- subset(dat, # newClst %in% c(1, 3, 10) & 
                   sample == "Normal CD8")

dfMeta <- datNor@meta.data

dfMeta <- dfMeta %>% 
  mutate(genotype = gsub("E8ICre.GFP", "", genotype),
         genotype = gsub("-yfp", "", genotype))

contig_list <- list()

for (i in unique(dfMeta$genotype)) {
  contig_list[[i]] <- dfMeta[dfMeta$genotype == i, 
                             c("barcode","is_cell", "contig_id", 
                               "high_confidence", "length", "chain", "v_gene", 
                               "d_gene", "j_gene", "c_gene", "full_length", 
                               "productive", "cdr3", "cdr3_nt", "reads", "umis", 
                               "raw_clonotype_id", "raw_consensus_id")]
}

combined <- combineTCR(contig_list, 
                       samples = names(contig_list))


# Quantify Clonotypes, the total or relative numbers of unique clonotypes
quantContig(combined, cloneCall="gene+nt")
quantContig(combined, cloneCall="aa")

# Clonotype Abundance,  examine the relative distribution of clonotypes by abundance.
# a line graph with a total number of clonotypes by the number of instances within the sample
abundanceContig(combined, cloneCall = "gene", scale = F)

# Length of Clonotypes, the length distribution of the CDR3 sequences
lengthContig(combined, cloneCall="aa", chain = "both") 

# Clonal Space Homeostasis, By examining the clonal space, we are effectively 
# looking at the relative space occupied by clones at specific proportions. 
# Another way to think about this would be thinking of the total immune receptor 
# sequencing run as a measuring cup. In this cup, we will fill liquids of 
# different viscosity - or different number of clonal proportions. Clonal space 
# homeostasis is asking what percentage of the cup is filled by clones in distinct 
# proportions (or liquids of different viscosity, to extend the analogy). 
# The proportional cutpoints are set under the cloneType variable in the function 
# and can be adjusted, at baseline the bins are as follows:
clonalHomeostasis(combined, cloneCall = "gene", 
                  cloneTypes = c(Rare = 1e-04, 
                                 Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1, 
                                 Hyperexpanded = 1))

clonalHomeostasis(combined, cloneCall = "aa")

# Clonal Proportion,  clonal proportion acts to place clones into separate bins. 
# The key difference is instead of looking at the relative proportion of the clone 
# to the total, the clonalProportion() function will rank the clones by total 
# number and place them into bins.
clonalProportion(combined, cloneCall = "gene",
                 split = c(10, 100, 1000, 10000, 30000, 1e+05)) 

# Overlap Analysis, 1) overlap coefficient, 2) Morisita index, or 3) Jaccard index. 
# The former is looking at the overlap of clonotypes scaled to the length of unique 
# clonotypes in the smaller sample. The Morisita index is more complex, it is an 
# ecological measure of the dispersion of individuals within a population, 
# incorporating the size of the population. The Jaccard Similarity Index is
# very similar to the overlap coefficient - instead of using the length of 
# the smaller sample, the denominator for the Jaccard Index is the union of 
# the two comparisons, leading to a much smaller number.
clonalOverlap(combined, cloneCall = "gene+nt", 
              method = "overlap") # "overlap", "morisita", "jaccard", "raw"

clonesizeDistribution(combined, cloneCall = "gene+nt", 
                      method="ward.D2")

# Diversity Analysis, Diversity can also be measured for samples or by other variables. Diversity is calculated using five metrics: 1) Shannon, 2) inverse Simpson, 3) Chao1, 4) Abundance-based Coverage Estimator (ACE), and 5) inverse Pielou's measure of species evenness. With the former two generally used to estimate baseline diversity and Chao/ACE indices used to estimate the richness of the samples. New implementation of this function include downsampling with 100 boot straps (n.boots) using the minimum number of unique clonotypes, as a more robust diversity estimate.
clonalDiversity(combined, cloneCall = "gene", group.by = "sample",
                exportTable = T, n.boots = 500)
