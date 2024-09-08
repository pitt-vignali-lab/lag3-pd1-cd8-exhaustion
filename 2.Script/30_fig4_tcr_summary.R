# ----------------------------------------------------------------------------
# Fig 4 | Summary statistics of the TCR repertoire
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/4.3.TCRComboSummary.R
# Input    : 1.Data/tcrCombinedUpdated.csv
# Output   : 3.Output/4.TCR/
# Upstream : 29_fig4_tcr_combined.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(scRepertoire)

dat <- readRDS("1.Data/datInt.RDS")

datNor <- subset(dat, sample == "Normal CD8")

trustDf <- read.csv("1.Data/TcrSplited/TRUST_merged_r2_report_modified.csv")[-1] %>%
  mutate(barcode = gsub("\\_|[0-9]", "", barcode),
         barcode = paste0("nor_", barcode),
         v_gene = gsub("\\*0[0-9]", "", v_gene),
         d_gene = gsub("\\*0[0-9]|\\.", "", d_gene),
         j_gene = gsub("\\*0[0-9]|\\.", "", j_gene),
         j_gene = gsub("\\*0[0-9]|\\.", "", j_gene)
  ) %>% 
  filter(cdr3 != "out_of_frame")

trustDf$c_gene <- substr(trustDf$j_gene, 1, 6)

for (i in 1:length(trustDf$c_gene)) {
  if(stringi::stri_sub(trustDf$c_gene[i], -1, -1) != "-"){
    trustDf$c_gene[i] = paste0(substr(trustDf$c_gene[i], 1, 3), "C")
  } else{
    trustDf$c_gene[i] = gsub("-", "", gsub("J", "C", trustDf$c_gene[i]))
  }
}



head(trustDf)


tcrDf <- read.csv("1.Data/TcrSplited/filtered_contig_annotations.csv") %>%
  mutate(barcode = gsub("\\-1", "", barcode),
         barcode = paste0("nor_", barcode)) %>% 
  select(colnames(trustDf))


cmbdDf <- rbind(trustDf, tcrDf) %>% 
  filter(barcode %in% rownames(datNor@meta.data)) %>% 
  mutate(d_gene = "") %>% 
  unique

cmbdDf$productive <- 'true'
# cmbdDf$reads <- 0
contig_list <- list()

for (i in unique(datNor$genotype)) {
  contig_list[[i]] <- cmbdDf %>% 
    filter(barcode %in% rownames(subset(datNor, genotype == i)@meta.data))
}

combined <- combineTCR(contig_list,
                       cells = "T-AB",
                       samples = c("E8ICre(WT)", "DK", "Pdcd1 K", "Lag3 K", "WT"), 
                       ID = c("WT", "KD", "KD", "KD", "WT"))
# Quantify Clonotypes
quantContig(combined, cloneCall="aa")

abundanceContig(combined, cloneCall = "aa", scale = FALSE, exportTable = T)

lengthContig(combined, cloneCall="aa")

clonalHomeostasis(combined, cloneCall = "gene", 
                  cloneTypes = c(Rare = 1e-04, 
                                 Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1, 
                                 Hyperexpanded = 1))

clonalHomeostasis(combined, cloneCall = "aa")

clonalProportion(combined, cloneCall = "gene",
                 split = c(10, 100, 1000, 10000, 30000, 1e+05)) 

clonalOverlap(combined, cloneCall = "gene+nt", 
              method = "overlap")

clonesizeDistribution(combined, cloneCall = "gene+nt", 
                      method="ward.D2")

clonalDiversity(combined, cloneCall = "aa",
                n.boots = 100)
