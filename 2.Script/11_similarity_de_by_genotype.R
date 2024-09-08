# ----------------------------------------------------------------------------
# Analysis | Same comparison broken down by genotype
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/5.4.similarityDeFunctionByGenotype.R
# Input    : 1.Data/datInt.RDS
# Output   : 3.Output/5.Similarity/De_PMel_Clst
#            3.Output/5.Similarity/De_Polyclonal_Clst
# Upstream : 10_similarity_de_function.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(presto)
library(ggvenn)
library(factoextra)

if(!dir.exists("3.Output/5.Similarity")){
  dir.create("3.Output/5.Similarity")
}

dat <- readRDS("1.Data/datInt.RDS")

# baseline 
datBase <- readRDS("1.Data/202104.RDS")

datBase <- subset(datBase, tissues == "ndLN")

datBase$sample <- "ndLN"

# "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"
# "Pdcd1L/L E8ICre.GFP"
# "Lag3L/L-yfp E8ICre.GFP"
# "Pdcd1L/L Lag3L/L-yfp (WT)" 

datSub <- subset(dat, genotype == "Pdcd1L/L E8ICre.GFP")

datNew <- merge(datSub, datBase)

datTest <- subset(datNew, sample %in% c("Normal CD8", "ndLN"))

Idents(datTest) <- "sample"

deGnvn <- wilcoxauc(datTest) %>% 
  filter(group == "Normal CD8")

datTest <- subset(datNew, sample %in% c("PML CD8", "ndLN"))

Idents(datTest) <- "sample"

deGpvn <- wilcoxauc(datTest) %>% 
  filter(group == "PML CD8")

# # VENN ----
# x <- list(
#   `Bulk CD8 vs ndLN UP` = deGnvn[deGnvn$logFC > 0 & deGnvn$padj < 0.05,]$feature,
#   `Bulk CD8 vs ndLN DN` = deGnvn[deGnvn$logFC < 0 & deGnvn$padj < 0.05,]$feature,
#   `pMEL CD8 vs ndLN UP` = deGpvn[deGpvn$logFC > 0 & deGpvn$padj < 0.05,]$feature,
#   `pMEL CD8 vs ndLN DN` = deGpvn[deGpvn$logFC < 0 & deGpvn$padj < 0.05,]$feature)
# 
# ggvenn(
#   x, 
#   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#   stroke_size = 0.5, set_name_size = 4
# ) + xlim(-5, 5)+
#   theme(text = element_text(size = 500, face = 2))
# 
# ggsave("3.Output/5.Similarity/DK_DE_gene(padj0.05)_VENN.png",
#        width = 25, height = 8, units = "in", dpi = 300, bg = "white")
# 
# # synergistic regulating genes ----
# synUp <- intersect(deGnvn[deGnvn$logFC > 0 & deGnvn$padj < 0.05,]$feature,
#                    deGpvn[deGpvn$logFC > 0 & deGpvn$padj < 0.05,]$feature)
# 
# synDn <- intersect(deGnvn[deGnvn$logFC < 0 & deGnvn$padj < 0.05,]$feature,
#                    deGpvn[deGpvn$logFC < 0 & deGpvn$padj < 0.05,]$feature)

# matrix generation ----
deMat <- data.frame(feature = rownames(dat@assays$RNA))

for (i in unique(Idents(dat))) {
  tryCatch({
  # ployclonal
  setNor <- subset(datSub, idents = i, sample == "Normal CD8")
  setNor <- merge(setNor, datBase)
  
  Idents(setNor) <- "sample"
  
  setNorDe <- wilcoxauc(setNor) %>% 
    filter(group == "Normal CD8") %>% 
    select(feature, logFC) 
  names(setNorDe)[2] = paste("norm", i, sep = "_") 
  deMat <- merge(deMat, setNorDe, by = "feature")
  
  # pMEL
  setPml <- subset(datSub, idents = i, sample == "PML CD8")
  setPml <- merge(setPml, datBase)
  
  Idents(setPml) <- "sample"
  
  setPmlDe <- wilcoxauc(setPml) %>% 
    filter(group == "PML CD8") %>% 
    select(feature, logFC) 
  names(setPmlDe)[2] = paste("pmel", i, sep = "_") 
  deMat <- merge(deMat, setPmlDe, by = "feature")
  }, error=function(e){})
}

deMatNew <- na.omit(deMat) %>%
  # filter(feature %in% c(synUp, synDn)) %>% 
  column_to_rownames(var = "feature") %>% 
  t

deMatNew <- deMatNew[, colSums(deMatNew) != 0]

# PCA ----
xxx.pca1<-prcomp(deMatNew, center=FALSE, scale.=FALSE, rank. = 4)
results <- xxx.pca1$x

# # K-means ----
# fviz_nbclust(results, FUNcluster=kmeans, k.max = 8)
# fviz_nbclust(results, FUNcluster=kmeans, method="gap_stat", k.max = 8)+ 
#   theme_bw() +
#   theme(text = element_text(size = 15, face = 2))
# 
# km1<-eclust(results, "kmeans", hc_metric="eucliden",k=7)
# 
# fviz_silhouette(km1)

dm<-dist(results)
hc<-hclust(dm, method="complete")
plot(hc, hang=-1, xlab = "", sub="")
rect.hclust(hc, k=6)

# manualy clustering
cluster <- c(rep(1, 5),
             rep(2, 4),
             rep(3, 5),
             rep(4, 4))

names(cluster) <- hc$labels[hc$order]
cluster <- cluster[hc$labels]
rect.hclust(hc, k=5, cluster = cluster)
