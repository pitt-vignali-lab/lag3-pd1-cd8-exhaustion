# ----------------------------------------------------------------------------
# Analysis | GSEA of the pseudobulk comparison
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/gsea.r
# Input    : 3.Output/PseudoBulk/
# Output   : 3.Output/PseudoBulk/GSEA_GOBP_*.csv
#            3.Output/6.GSEA/KEGGsig.csv
# Upstream : 12_pseudobulk_degs.R
# Note     : msigdbr 7.5.1 gene sets
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/6.GSEA")){
  dir.create("3.Output/6.GSEA")
}

library(Seurat)
library(presto)
library(fgsea)
library(msigdbr)
library(tidyverse)

# data read in ----
dat <- readRDS("1.Data/datInt.RDS")

datNor <- subset(dat, sample == "Normal CD8")

datPml <- subset(dat, sample == "PML CD8")

# GSEA(GOBP)----
m_df<- msigdbr(species = "Mus musculus", 
               category = "C5",
               subcategory = "BP")

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)


unique(datPml$genotype)

goBpDf <- data.frame()

for (genotype in c("Pdcd1L/L Lag3L/L-yfp E8ICre.GFP",
                   "Pdcd1L/L E8ICre.GFP",
                   "Lag3L/L-yfp E8ICre.GFP")) {
  
  rank <- wilcoxauc(datNor, group_by = "genotype", 
                     groups_use = c(genotype,
                                    "Pdcd1L/L Lag3L/L-yfp (WT)")) %>% 
    filter(group == genotype) %>% 
    arrange(auc) %>% 
    select(feature, auc) %>% 
    deframe()
  
  fgseaRes<- fgsea(fgsea_sets, stats = rank, nperm = 1000) %>%
    as_tibble() %>%
    arrange(desc(NES)) %>% 
    select(-leadingEdge, -ES, -nMoreExtreme) %>%
    filter(padj < 0.05) %>% 
    add_column(genotype = genotype)
  
  goBpDf <- rbind(goBpDf, fgseaRes)
}

write.csv(goBpDf, "3.Output/6.GSEA/KEGGsig.csv")

commonPtw <- goBpDf$pathway[table(goBpDf$pathway) == 1]

ggplot(goBpDf[goBpDf$pathway %in% commonPtw,] %>%  head, aes(genotype, pathway)) +
  geom_point()

# C8 catagory ----
m_df<- msigdbr(species = "Mus musculus", 
               category = "C8")

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

BpDf <- data.frame()

for (clst in unique(Idents(dat))) {
  
  rank <- wilcoxauc(dat) %>% 
    filter(group == clst) %>% 
    arrange(logFC) %>% 
    select(feature, logFC) %>% 
    deframe()
  
  fgseaRes<- fgsea(fgsea_sets, stats = rank, nperm = 1000) %>%
    as_tibble() %>%
    arrange(desc(NES)) %>% 
    select(-leadingEdge, -ES, -nMoreExtreme) %>%
    filter(padj < 0.05) %>% 
    add_column(cluster = clst)
  
  BpDf <- rbind(BpDf, fgseaRes)
}

t <- BpDf %>% filter(NES > 0) %>% 
  group_by(cluster) %>% 
  top_n(5, NES)
