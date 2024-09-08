# ----------------------------------------------------------------------------
# Analysis | Differential expression behind the similarity structure
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/5.3.similarityDeFunction.R
# Input    : 1.Data/datInt.RDS
# Output   : 3.Output/5.Similarity/OverallDeG*.csv
# Upstream : 09_fig2_similarity_pca.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(presto)

dat <- readRDS("1.Data/datInt.RDS")

datNor <- readRDS("1.Data/202104.RDS")

baseLine <- subset(datNor, tissues == "ndLN")

baseLine$tem <- "base"

# DE 
for (i in c(1, 3)) {
  # ployclonal
  setNor <- subset(dat, idents = i, sample == "Normal CD8")
  setNor$tem <- "target"
  setNor <- merge(setNor, baseLine)
  Idents(setNor) <- "tem"
  
  setNorDe <- wilcoxauc(setNor) %>% 
    filter(group == "target") %>% 
    write.csv(paste0("3.Output/5.Similarity/De_Polyclonal_Clst", 
                     i, "_tumor_vs_ndLN.csv"))
  
  # pMEL
  setPml <- subset(dat, idents = i, sample == "PML CD8")
  setPml$tem <- "target"
  setPml <- merge(setPml, baseLine)
  
  Idents(setPml) <- "tem"
  
  setPmlDe <- wilcoxauc(setPml) %>% 
    filter(group == "target")%>% 
    write.csv(paste0("3.Output/5.Similarity/De_PMel_Clst", 
                     i, "_tumor_vs_ndLN.csv"))
}
