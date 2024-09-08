# ----------------------------------------------------------------------------
# Analysis | Pseudobulk differential expression across genotypes
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/6.tonyDEGs.R
# Input    : 1.Data/datInt.RDS
# Output   : 3.Output/PseudoBulk/upGforORA.csv
#            3.Output/PseudoBulk/dnGforORA.csv
# Upstream : 02_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(readxl)
library(presto)

shN <- excel_sheets("1.Data/18071_CD8_TIL_DEGs_by_treatment_220425.xlsx")

for (i in shN) {
  t <- read_excel("1.Data/18071_CD8_TIL_DEGs_by_treatment_220425.xlsx",
                       sheet = i)
  cat(i, "number:", nrow(t), 
      "padj:", max(t$p_val_adj), 
      "logFC:", range(t$avg_log2FC),
      "PC1:", range(t$pct.1),
      "PC2:", range(t$pct.2), "\n")
}

# threshold: logFC > 0.25, pc > 0.1

dat <- readRDS("1.Data/datInt.RDS")

datSub <- subset(dat, sample == "PML CD8" & genotype != "Pdcd1L/L Lag3L/L-yfp (WT)")

Idents(datSub) <- "genotype"

# using same threshold in seurat ----
deG <- FindAllMarkers(datSub) %>% 
  filter(p_val_adj < 0.05)



sapply(shN, function(x) 
  nrow(read_excel("1.Data/18071_CD8_TIL_DEGs_by_treatment_220425.xlsx",
                  sheet = x)))[-3] %>%
  data.frame %>% 
  rownames_to_column(var = 'Genotype') %>% 
  mutate(Genotype = factor(Genotype, levels = c("LAG3_Tx", "PD1_Tx", "Combo_Tx"))) %>% 
  rename(Freq = 2) %>% 
  ggplot(aes(x = Genotype, y = Freq)) +
  geom_bar(stat = "identity") +
  ggtitle("Tony's gene set") +
  theme_minimal() +
  deG %>% 
  filter(p_val_adj < 0.05) %>% 
  select(cluster) %>% 
  table %>% 
  data.frame %>% 
  rename(Genotype = 1) %>% 
  mutate(Genotype = gsub("E8ICre.GFP", "", Genotype),
         Genotype = gsub("-yfp", "", Genotype)) %>% 
  ggplot(aes(x = Genotype, y = Freq)) +
  ggtitle("Lawrence's gene set") +
  geom_bar(stat = "identity") +
  theme_minimal()
    


p1 + p2

for (i in unique(deG$cluster)) {
  
  t <- deG[deG$cluster == i & deG$p_val_adj < 0.05,]
  cat(i, "number:", nrow(t), 
      "padj:", max(t$p_val_adj), 
      "logFC:", range(t$avg_log2FC),
      "PC1:", range(t$pct.1),
      "PC2:", range(t$pct.2), "\n")
}

# same threshold with wilcoauc ----
deG <- wilcoxauc(datSub)

for (i in unique(deG$group)) {
  deGSub <- deG[deG$group == i,]
  deGSub <- deGSub[deGSub$pct_in > 0.1 | deGSub$pct_out > 0.1,]
  deGSub <- deGSub[deGSub$logFC > 0.25,]
  
  cat(i, range(deGSub$logFC), range(deGSub$pct_in), range(deGSub$pct_out), nrow(deGSub), "\n")
}





thre <- sapply(shN, function(x) 
  max(read_excel("1.Data/18071_CD8_TIL_DEGs_by_treatment_220425.xlsx",
             sheet = x)$p_val_adj))

gNo <- sapply(shN, function(x) 
  nrow(read_excel("1.Data/18071_CD8_TIL_DEGs_by_treatment_220425.xlsx",
                 sheet = x)))

gNo

barplot(gNo[-3])

threSub <- thre[-3]

for (i in 1:length(unique(lDegs$group))) {
  cat("Lawrence:", unique(lDegs$group)[i],
      nrow(lDegs[lDegs$group == unique(lDegs$group)[i] & 
                   lDegs$padj <= threSub[i],]),
      max(lDegs[lDegs$group == unique(lDegs$group)[i] & 
                   lDegs$padj <= threSub[i],]$padj),
      "Tony:",
      names(threSub[i]),threSub[i], "\n")
}

gNo <- data.frame(gNo[-3])



