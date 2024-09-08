# ----------------------------------------------------------------------------
# Fig 4 | Load and merge CellRanger vdj and TRUST4 TCR calls
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/4.1.1.TCR.R
# Input    : 1.Data/filtered_contig_annotations.csv
#            1.Data/Trust4Output/TRUST_merged_r2_report.tsv
# Output   : 1.Data/tcrCombined.csv
# Upstream : 01_data_integration.R
# Note     : CellRanger vdj plus TRUST4 1.0.12
#            duplicate sequences trimmed
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

if(!dir.exists("3.Output/4.TCR")){
  dir.create("3.Output/4.TCR")
}

dat <- readRDS("1.Data/datInt.RDS")

datNor <- subset(dat, sample == "Normal CD8")

# TCR frequency summary ----
barplot(table(datNor$cdr3))

cdr3Df <- data.frame(table(datNor$cdr3))

barp <- ggplot(cdr3Df[cdr3Df$Freq > 20,], aes(y = reorder(Var1, Freq), x = Freq)) +
  geom_bar(stat = "identity") +
  labs(y = "Clonaltype", x = "Freqence" ) +
  theme_minimal()+
  theme(text = element_text(size = 10, face = 2))

ggsave("3.Output/4.TCR/TCR_clonaltypes_freqence.png", width = 6, height = 3,
       units = "in", dpi = 300, bg = "white")


targetTcr <- cdr3Df[order(cdr3Df$Freq, decreasing = T),]$Var1[1:3]

# Tcr distribution ----

plotDf <- Reduce(function(x, y) merge(x, y , by = "cells"),
                 list(data.frame(Embeddings(datNor, reduction = "umap")) %>% 
                        rownames_to_column(var = "cells"),
                      data.frame(Idents(datNor))%>% 
                        rownames_to_column(var = "cells"),
                      datNor@meta.data["cdr3"]%>% 
                        rownames_to_column(var = "cells")))

labDf <- plotDf %>% select("UMAP_1", "UMAP_2", "Idents.datNor.") %>% 
  group_by(Idents.datNor.) %>% 
  summarise(UMAP_1 = mean(UMAP_1), 
            UMAP_2 = mean(UMAP_2))

ggplot(data= plotDf[c(1:4)], aes(UMAP_1, UMAP_2, fill = Idents.datNor.)) +
  geom_point(shape = 21, alpha = 0.1, size = 5, col = "gray50") +
  geom_point(data = plotDf[plotDf$cdr3 %in% targetTcr,] %>% 
               mutate(cdr3 = factor(cdr3, levels = c("CASSPPGGERLFF", "CSSRDGGSTEVFF", "CASSLELGGREQYF"))),
             aes(UMAP_1, UMAP_2), col = "#231f20", size = 2, alpha = 0.5)+
  geom_text(data = labDf, 
            aes(UMAP_1, UMAP_2, label = Idents.datNor.),
            size = 8, fontface = 2) +
  facet_grid(~cdr3) +
  coord_fixed(ratio = 1) +
  # scale_color_brewer(palette="Set2") +
  guides(fill = "none") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 15, face = 2))

ggsave("3.Output/4.TCR/TCR_distribution.png", width = 12, height = 6, units = "in",
       dpi = 300, bg = "white")

table(Idents(subset(dat, cdr3 %in% targetTcr)))

# updating clusters ----

newClst <- dat@meta.data["integrated_snn_res.0.35"]

names(newClst) <- "newClst"
 
newClst$newClst <- ifelse(rownames(newClst) %in% Cells(subset(datInt, idents = c(1, 2))),
                          11,
                          newClst$newClst)

newClst$newClst <- newClst$newClst-1

dat <- AddMetaData(dat, newClst)

dat <- AddMetaData(dat, read.csv("1.Data/TCR_defined_clusters.csv", row.names = 1) %>% rename(newClst = dat.newClst))

datNor <- subset(dat, sample == "Normal CD8")

plotDf <- Reduce(function(x, y) merge(x, y , by = "cells"),
                 list(data.frame(Embeddings(datNor, reduction = "umap")) %>% 
                        rownames_to_column(var = "cells"),
                      data.frame(datNor$newClst)%>% 
                        rownames_to_column(var = "cells"),
                      data.frame(datNor$genotype)%>% 
                        rownames_to_column(var = "cells"),
                      datNor@meta.data["cdr3"]%>% 
                        rownames_to_column(var = "cells")))

labDf <- plotDf %>% select("UMAP_1", "UMAP_2", "datNor.newClst") %>% 
  group_by(datNor.newClst) %>% 
  summarise(UMAP_1 = mean(UMAP_1), 
            UMAP_2 = mean(UMAP_2))

plotDf$datNor.genotype <- factor(plotDf$datNor.genotype, 
                          levels = c("E8ICre.GFP (WT)",
                                     "Pdcd1L/L Lag3L/L-yfp (WT)",
                                     "Pdcd1L/L E8ICre.GFP",
                                     "Lag3L/L-yfp E8ICre.GFP",
                                     "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))

ggplot() +
  geom_point(data= plotDf[-6] %>% mutate(datNor.newClst = as.character(datNor.newClst)), 
             aes(UMAP_1, UMAP_2, fill = datNor.newClst), 
             shape = 21, alpha = 0.1, size = 5, col = "gray50") +
  geom_point(data = plotDf[plotDf$cdr3 %in% targetTcr,] %>%
               mutate(cdr3 = factor(cdr3, levels = c("CASSPPGGERLFF", "CSSRDGGSTEVFF", "CASSLELGGREQYF"))),
             aes(UMAP_1, UMAP_2), col = "#231f20", size = 2, alpha = 0.5)+
  geom_text(data = labDf,
            aes(UMAP_1, UMAP_2, label = datNor.newClst),
            size = 5, fontface = 2) +
  facet_grid( cdr3 ~ datNor.genotype) +
  coord_fixed(ratio = 1) +
  # scale_color_brewer(palette="Set2") +
  guides(fill = "none") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 15, face = 2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2))

ggsave("3.Output/4.TCR/TCR_distribution_updated_genotype.png", width = 15, height = 12, units = "in",
       dpi = 300, bg = "white")

table(subset(dat, cdr3 %in% targetTcr)$newClst)

write.csv(data.frame(dat$newClst), "1.Data/TCR_defined_clusters.csv")

DimPlot(dat, group.by = "newClst")

# dat <- FindClusters(dat, resolution = 0.4)

datNor <- subset(dat, sample == "Normal CD8")

DimPlot(datNor, label = T, label.size = 5)

Idents(dat) <- 'newClust'

table(Idents(subset(dat, cdr3 %in% c("CASSPPGGERLFF", "CSSRDGGSTEVFF"))))

t <- data.frame(cdr3 = table(datNor$cdr3_nt))

t[order(t$cdr3.Freq, decreasing = T),] %>% head

DimPlot(datNor, 
        cells.highlight = rownames(datNor@meta.data[datNor@meta.data$cdr3 %in% c("CASSPPGGERLFF", "CSSRDGGSTEVFF", "CASSLELGGREQYF"),]),
        label = T)



ggplot(plotDf, aes(x = UMAP_1, y = UMAP_2, fill = Idents.datNor.)) + 
  stat_density_2d(geom = "polygon")


ggplot(plotDf[plotDf$Idents.datNor. %in% c(2, 3, 4, 6),], 
       aes(x = UMAP_1, y = UMAP_2, fill = Idents.datNor.)) +
  stat_density_2d(geom="polygon", bins=10, alpha=0.2)+
  geom_point(data = plotDf[plotDf$cdr3 %in% c("CASSPPGGERLFF", "CSSRDGGSTEVFF"),],
             aes(UMAP_1, UMAP_2), col = "black")

table(plotDf[plotDf$cdr3 %in% c("CASSPPGGERLFF", "CSSRDGGSTEVFF"),]$Idents.datNor.)
