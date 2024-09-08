# ----------------------------------------------------------------------------
# Fig 2 | PCA and k-means clustering of sample similarity
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/5.2.similarityTestPCA.R
# Input    : 1.Data/datInt.RDS
# Output   : 3.Output/5.Similarity/kmeansCluster*.png
# Upstream : 08_fig2_similarity_key_features.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(presto)
library(factoextra)

# synergistic regulating genes ----
synUp <- intersect(read.csv("3.Output/5.Similarity/OverallDeGNorvsNdLN.csv") %>% 
                     filter(logFC > 0) %>% 
                     select(feature) %>% 
                     deframe(),
                   read.csv("3.Output/5.Similarity/OverallDeGpMlvsNdLN.csv") %>% 
                     filter(logFC > 0) %>% 
                     select(feature) %>% 
                     deframe())

synDn <- intersect(read.csv("3.Output/5.Similarity/OverallDeGNorvsNdLN.csv") %>% 
                     filter(logFC < 0) %>% 
                     select(feature) %>% 
                     deframe(),
                   read.csv("3.Output/5.Similarity/OverallDeGpMlvsNdLN.csv") %>% 
                     filter(logFC < 0) %>% 
                     select(feature) %>% 
                     deframe())


# matrix generation ----
datNor <- readRDS("1.Data/202104.RDS")

baseLine <- subset(datNor, tissues == "ndLN")

baseLine$tem <- "base"

dat <- readRDS("1.Data/datInt.RDS")


deMat <- data.frame(feature = rownames(dat@assays$RNA))
for (i in unique(Idents(dat))) {
  
  # ployclonal
  setNor <- subset(dat, idents = i, sample == "Normal CD8")
  setNor$tem <- "target"
  setNor <- merge(setNor, baseLine)
  
  Idents(setNor) <- "tem"
  
  setNorDe <- wilcoxauc(setNor) %>% 
    filter(group == "target") %>% 
    select(feature, logFC) 
  names(setNorDe)[2] = paste("norm", i, sep = "_") 
  deMat <- merge(deMat, setNorDe, by = "feature")
  
  # pMEL
  setPml <- subset(dat, idents = i, sample == "PML CD8")
  setPml$tem <- "target"
  setPml <- merge(setPml, baseLine)
  
  Idents(setPml) <- "tem"
  
  setPmlDe <- wilcoxauc(setPml) %>% 
    filter(group == "target") %>% 
    select(feature, logFC) 
  names(setPmlDe)[2] = paste("pmel", i, sep = "_") 
  deMat <- merge(deMat, setPmlDe, by = "feature")
  
}

deMatNew <- na.omit(deMat) %>%
  filter(feature %in% c(synUp, synDn)) %>% 
  column_to_rownames(var = "feature") %>% 
  t


# PCA ----
xxx.pca1<-prcomp(deMatNew, center=FALSE, scale.=FALSE, rank. = 4)
results <- xxx.pca1$x

# K-means ----
fviz_nbclust(results, FUNcluster=kmeans, k.max = 8)
fviz_nbclust(results, FUNcluster=kmeans, method="gap_stat", k.max = 8)+ 
  theme_bw() +
  theme(text = element_text(size = 15, face = 2))
# ggsave("3.Output/5.Similarity/kmeansClusterNoGapStatis.png", width = 6,
#        height = 3, dpi = 300, bg = "white")

km1<-eclust(results, "kmeans", hc_metric="eucliden",k=5)
# ggsave("3.Output/5.Similarity/kmeansCluster.png", width = 6,
#        height = 5, dpi = 300, bg = "white")

fviz_silhouette(km1)
# ggsave("3.Output/5.Similarity/kmeansClusterSilhouette.png", width = 6,
#        height = 5, dpi = 300, bg = "white")

# Hierarchical Clustering ----
dm<-dist(results)
hc<-hclust(dm, method="complete")

# png("3.Output/5.Similarity/kmeansClusterHierarchical.png", width = 6,
#     height = 5, units = "in", res = 300)
plot(hc, hang=-1, xlab = "", sub="")
rect.hclust(hc, k=5, border=c("red", "pink", "red", "pink", "pink"))
# dev.off()

