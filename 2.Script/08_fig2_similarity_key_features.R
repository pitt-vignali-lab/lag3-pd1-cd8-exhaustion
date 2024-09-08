# ----------------------------------------------------------------------------
# Fig 2 | Key features distinguishing the genotypes
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/5.1.similarityKeyFeat.R
# Input    : 1.Data/datInt.RDS
# Output   : 3.Output/5.Similarity/
# Upstream : 02_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(presto)

if(!dir.exists("3.Output/5.Similarity")){
  dir.create("3.Output/5.Similarity")
}

dat <- readRDS("1.Data/datInt.RDS")

table(dat$sample)

newClust <- read.csv("1.Data/TCR_defined_clusters.csv", row.names = 1)

names(newClust) <- "newClst"

dat <- AddMetaData(dat, newClust)

datSub <- subset(dat, newClst %in% c(1, 3, 10))

DimPlot(datSub, group.by = "newClst", split.by = "genotype") 

datBase <- readRDS("1.Data/202104.RDS")

datBase <- subset(datBase, tissues == "ndLN")

datBase$sample <- "ndLN"

# DE 

datNew <- merge(datSub, datBase)

datTest <- subset(datNew, sample %in% c("Normal CD8", "ndLN"))

Idents(datTest) <- "sample"

deGnvn <- wilcoxauc(datTest) %>% 
  filter(group == "Normal CD8")

datTest <- subset(datNew, sample %in% c("PML CD8", "ndLN"))

Idents(datTest) <- "sample"

deGpvn <- wilcoxauc(datTest) %>% 
  filter(group == "PML CD8")

# VENN ----

library(ggvenn)

x <- list(
  `Bulk CD8 vs ndLN UP` = deGnvn[deGnvn$logFC > 0 & deGnvn$padj < 0.05,]$feature,
  `Bulk CD8 vs ndLN DN` = deGnvn[deGnvn$logFC < 0 & deGnvn$padj < 0.05,]$feature,
  `pMEL CD8 vs ndLN UP` = deGpvn[deGpvn$logFC > 0 & deGpvn$padj < 0.05,]$feature,
  `pMEL CD8 vs ndLN DN` = deGpvn[deGpvn$logFC < 0 & deGpvn$padj < 0.05,]$feature)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + xlim(-5, 5)+
  theme(text = element_text(size = 500, face = 2))

ggsave("3.Output/5.Similarity/DE_gene(padj0.05)_VENN.png",
       width = 25, height = 8, units = "in", dpi = 300, bg = "white")

# plot distributions 

plotDf <- merge(Embeddings(datSub, "umap"),
                datSub@meta.data[c("sample", "genotype", "newClst")],
                by = 0)

plotDf$genotype <- factor(plotDf$genotype, 
                          levels = c("E8ICre.GFP (WT)",
                                     "Pdcd1L/L Lag3L/L-yfp (WT)",
                                     "Pdcd1L/L E8ICre.GFP",
                                     "Lag3L/L-yfp E8ICre.GFP",
                                     "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))

labDf <- plotDf %>% select("UMAP_1", "UMAP_2", "newClst") %>% 
  group_by(newClst) %>% 
  summarise(UMAP_1 = mean(UMAP_1), 
            UMAP_2 = mean(UMAP_2))

ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2, col = factor(newClst, levels = c(1, 3, 10)))) +
  geom_point(alpha = 0.6) +
  facet_grid(~ sample) +
  coord_fixed(ratio = 1) +
  geom_text(data = labDf, 
            aes(UMAP_1, UMAP_2, label = newClst),
            size = 8, fontface = 2, col = "black") +
  theme_minimal() +
  labs(col = "Clusters") +
  theme(text = element_text(size = 15, face = 2),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))

ggsave("3.Output/5.Similarity/cell_distribution_in_either_group.png",
       width = 8, height = 5, units = "in", dpi = 300, bg = "white")

ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2, col = factor(newClst, levels = c(1, 3, 10)))) +
  geom_point(alpha = 0.6) +
  facet_grid(sample ~ genotype) +
  coord_fixed(ratio = 1) +
  geom_text(data = labDf, 
            aes(UMAP_1, UMAP_2, label = newClst),
            size = 8, fontface = 2, col = "black") +
  theme_minimal() +
  labs(col = "Clusters") +
  theme(text = element_text(size = 15, face = 2),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))

ggsave("3.Output/5.Similarity/cell_distribution_in_either_group_genotyp.png",
       width = 15, height = 8, units = "in", dpi = 300, bg = "white")
