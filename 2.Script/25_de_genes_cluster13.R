# ----------------------------------------------------------------------------
# Analysis | Differential expression in clusters 1 and 3
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/rev6_DeGin13.R
# Input    : 1.Data/datForRevision.RDS
# Output   : 3.Output/Revision/sigDeGensDKOvsWTClst*.csv
# Upstream : 04_cluster_annotation.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(presto)

dat <- readRDS('1.Data/datForRevision.RDS')

datSub <- subset(dat, 
                 idents = 1,
                 genotype %in% c("E8ICre.GFP (WT)", "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", "Pdcd1L/L Lag3L/L-yfp (WT)"))


datSub$newGeno <- 'DKO'


datSub$newGeno[grep('WT', datSub$genotype)] <- 'WT'

deG1 <- wilcoxauc(datSub, group_by = 'newGeno') %>% 
  filter(group == 'DKO',
         padj < 0.05)
write.csv(deG1, '3.Output/Revision/sigDeGensDKOvsWTClst1.csv')

datSub <- subset(dat, 
                 idents = 3,
                 genotype %in% c("E8ICre.GFP (WT)", "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", "Pdcd1L/L Lag3L/L-yfp (WT)"))


datSub$newGeno <- 'DKO'


datSub$newGeno[grep('WT', datSub$genotype)] <- 'WT'

deG3 <- wilcoxauc(datSub, group_by = 'newGeno') %>% 
  filter(group == 'DKO',
         padj < 0.05)
write.csv(deG3, '3.Output/Revision/sigDeGensDKOvsWTClst3.csv')





datSub$sampGeno <- paste0(datSub$sample, datSub$genotype) 

Idents(datSub) <- 'sampGeno'


deG <- wilcoxauc(datSub) %>% 
  filter(padj < 0.01)



deG64 <- wilcoxauc(datSub64) %>% 
  filter(padj < 0.05,
         group == 4)

datSub64$sampGeno <- paste(datSub64$sample, datSub64$genotype)
datSub64$sampGenoAnno <- paste(datSub64$sample, datSub64$genotype, datSub64$anno)

datSce <- as.SingleCellExperiment(datSub64)

milo <- Milo(datSce)

# reducedDim(milo, "UMAP") <- reducedDim(datSce, "UMAP")

milo

milo <- buildGraph(milo, k = 50, d = 30)

milo <- makeNhoods(milo, prop = 0.2, 
                   k = 50, 
                   d = 30, 
                   refined = TRUE)

plotNhoodSizeHist(milo)

milo <- countCells(milo, 
                   meta.data = data.frame(colData(milo)), 
                   sample="sampGenoAnno")

head(nhoodCounts(milo))

designDf <- data.frame(colData(milo))[,c('sample', 'genotype', 
                                         'anno', 'sampGeno', 'sampGenoAnno')] %>% 
  distinct() %>% 
  remove_rownames() %>% 
  # mutate(tissueType = paste(tissue, anno, sep = '_')) %>% 
  column_to_rownames(var = "sampGenoAnno") 

# designDf <- designDf[colnames(nhoodCounts(milo)), , drop=FALSE]

# designDf$sampGenoAnno <- rownames(designDf)

designDf

milo <- calcNhoodDistance(milo, d = 30)

# DA test 
da_results <- testNhoods(milo, 
                         design = ~ anno, 
                         design.df = designDf)
head(da_results)

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

milo <- buildNhoodGraph(milo)

plotReducedDim(milo, dimred = "UMAP", colour_by="anno", 
               text_size = 3, point_size=0.5) +
  guides(fill="none") +
  coord_fixed() +
  
  plotNhoodGraphDA(milo, 
                   da_results, 
                   layout="UMAP",
                   alpha=0.1) +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0) +
  labs(fill = 'logFC') +
  coord_fixed() 

da_results <- annotateNhoods(milo, 
                             da_results, 
                             coldata_col = "anno")


da_results$NhoodGroup <- da_results$anno

da_nhood_markers <- findNhoodGroupMarkers(milo, 
                                          da_results, 
                                          subset.row = rownames(milo))

mk <- da_nhood_markers %>% 
  filter(`adj.P.Val_Tumor primed-TCM` < 0.05)

inteG <- intersect(mk$GeneID, deG64$feature)

inteG[-grep('^Hist', inteG)]

write.csv(deG64[deG64$feature %in% inteG[-grep('^Hist', inteG)],],
          '3.Output/Revision/milo&DeOfClst64.csv')





dat$sampGeno <- paste(dat$sample, dat$genotype)
dat$sampGenoAnno <- paste(dat$sample, dat$genotype, dat$anno)

# dat <- subset(dat, sample == 'Normal CD8')

datSce <- as.SingleCellExperiment(dat)

milo <- Milo(datSce)

# reducedDim(milo, "UMAP") <- reducedDim(datSce, "UMAP")

milo

milo <- buildGraph(milo, k = 50, d = 30)

milo <- makeNhoods(milo, prop = 0.2, 
                   k = 50, 
                   d = 30, 
                   refined = TRUE)

plotNhoodSizeHist(milo)

milo <- countCells(milo, 
                   meta.data = data.frame(colData(milo)), 
                   sample="sampGenoAnno")

head(nhoodCounts(milo))

designDf <- data.frame(colData(milo))[,c('sample', 'genotype', 
                                         'anno', 'sampGeno', 'sampGenoAnno')] %>% 
  distinct() %>% 
  remove_rownames() %>% 
  # mutate(tissueType = paste(tissue, anno, sep = '_')) %>% 
  column_to_rownames(var = "sampGenoAnno") 

# designDf <- designDf[colnames(nhoodCounts(milo)), , drop=FALSE]

# designDf$sampGenoAnno <- rownames(designDf)

designDf

milo <- calcNhoodDistance(milo, d = 30)

# DA test 
da_results <- testNhoods(milo, 
                         design = ~ genotype, 
                         design.df = designDf)
head(da_results)

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

milo <- buildNhoodGraph(milo)

plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.15) +
  scale_fill_gradient2(high = 'red', mid = 'white', low = 'blue', midpoint = 0) +
  labs(fill = 'logFC') +
  coord_fixed() 


# milo <- logNormCounts(milo)

da_results <- annotateNhoods(milo, 
                             da_results, 
                             coldata_col = "sample")


da_results$NhoodGroup <- da_results$sample

da_nhood_markers <- findNhoodGroupMarkers(milo, 
                                          da_results, 
                                          subset.row = rownames(milo))

mk <- da_nhood_markers %>% 
  filter(`adj.P.Val_PML CD8` < 0.05)