# ----------------------------------------------------------------------------
# Fig 2 | Render the published volcano panels
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/PictureFroManual/volcano.R
# Input    : 1.Data/datInt.RDS
# Output   : 3.Output/FormalFigures/volcano_*.png
# Upstream : 02_cluster_labeling.R
# ----------------------------------------------------------------------------

library(Seurat)
library(presto)
library(tidyverse)
library(ggrepel)
library(extrafont)

windowsFonts("Arial" = windowsFont("Arial"))

dat <- readRDS("1.Data/datInt.RDS")

# cluster 7 ----
datSub <- subset(dat, genotype == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP")

de <- wilcoxauc(datSub)

de <- de %>% filter(group == 7, logFC != 0, padj < 0.05)

tarFea <- strsplit("Ccl5, Gzmk, Ifit3, Cxcr3, Ifngr1", split = ", ")[[1]]
topFea <- c(de %>% arrange(desc(logFC)) %>% head(5) %>% dplyr::select(feature) %>% deframe(),
            de %>% arrange(desc(logFC)) %>% tail(5) %>% dplyr::select(feature) %>% deframe(),
            de %>% arrange(desc(padj)) %>% tail(5) %>% dplyr::select(feature) %>% deframe())

selecFea <- unique(c(tarFea, topFea))

ggplot(de, aes(logFC, -log10(padj))) +
  geom_point(aes(col = ifelse(logFC > 0, "up", "dn")),
             alpha = 0.3, size = 1.5)+
  geom_label_repel(aes(label = ifelse(feature %in% selecFea, feature, ""),
                       col = ifelse(logFC > 0, "up", "dn"),
                       fontface = ifelse(feature %in% tarFea, 2, 1)),
                   alpha = ifelse(de$feature %in% tarFea, 1, 0.5),
                   max.overlaps = 1000, show.legend = F) +
  labs(x = bquote(bold(log^{fold~change})),
       y = bquote(bold(Significance~(-log^{adjuste~italic(p)~value})))) +
  scale_y_continuous(trans = scales::log10_trans())+
  scale_color_manual(name = "Gene staus",
                     breaks = c("up", "dn", "none-sig"),
                     labels = c("Significantly up-regulated genes",
                                "Significantly down-regulated genes", 
                                "none-sig"),
                     limits = c("up","dn"),
                     values=c("red", "blue", "gray50"))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "Arial", size = 10, face = 2))
ggsave("3.Output/FormalFigures/volcano_Clustr7_2.png",  
       width = 6, height = 3, units = "in", dpi = 300, bg = "white")

# cluster 0 ----
datSub <- subset(dat, genotype == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP")

de <- wilcoxauc(datSub)

de <- de %>% filter(group == 0, logFC != 0, padj < 0.05)

tarFea <- strsplit("Ccl5, Gzmb, Ifit3, Cxcr6, Ifngr1, Nkg7", split = ", ")[[1]]
topFea <- c(de %>% arrange(desc(logFC)) %>% head(5) %>% select(feature) %>% deframe(),
            de %>% arrange(desc(logFC)) %>% tail(5) %>% select(feature) %>% deframe(),
            de %>% arrange(desc(padj)) %>% tail(5) %>% select(feature) %>% deframe())

selecFea <- unique(c(tarFea, topFea))

ggplot(de, aes(logFC, -log10(padj))) +
  geom_point(aes(col = ifelse(logFC > 0, "up", "dn")),
             alpha = 0.3, size = 1.5)+
  geom_label_repel(aes(label = ifelse(feature %in% selecFea, feature, ""),
                       col = ifelse(logFC > 0, "up", "dn"),
                       fontface = ifelse(feature %in% tarFea, 2, 1)),
                   alpha = ifelse(de$feature %in% tarFea, 1, 0.5),
                   max.overlaps = 1000, show.legend = F) +
  labs(x = bquote(bold(log^{fold~change})),
       y = bquote(bold(Significance~(-log^{adjuste~italic(p)~value})))) +
  scale_y_continuous(trans = scales::log10_trans())+
  scale_color_manual(name = "Gene staus",
                     breaks = c("up", "dn", "none-sig"),
                     labels = c("Significantly up-regulated genes",
                                "Significantly down-regulated genes", 
                                "none-sig"),
                     limits = c("up","dn"),
                     values=c("red", "blue", "gray50"))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "Arial", size = 10, face = 2))
ggsave("3.Output/FormalFigures/volcano_Clustr0_2.png",  
       width = 6, height = 4, units = "in", dpi = 300, bg = "white")

# comparison between WT vs DN in normal and pml----
tarFea <- c("Klrc1", "Klrd1")
selecFea <- c("Klrc1", "Lag3", "Klrd1", "Pdcd1","Cd94", "Klrc2", "Klre1", "Klrc3")

# single wt
datPoly <- subset(dat, 
                  genotype %in% c(# "E8ICre.GFP (WT)",  
                                  "Pdcd1L/L Lag3L/L-yfp (WT)",
                                  "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP") &
                    sample == "Normal CD8")

# combo wt
datPoly <- subset(dat, 
                  genotype %in% c("E8ICre.GFP (WT)", 
                                  "Pdcd1L/L Lag3L/L-yfp (WT)",
                                  "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP") &
                    sample == "Normal CD8")
datPoly$genotype[grep("WT", datPoly$genotype)] <- "WT"

# plotting 
table(datPoly$genotype)

Idents(datPoly) <- "genotype"

dePoly <- wilcoxauc(datPoly) %>% 
  filter(group == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP",
         padj < 0.05)


ggplot(dePoly, aes(x = logFC, y = -log10(padj))) +
  geom_point(aes(col = ifelse(logFC > 0, "up", "dn")),
             alpha = 0.5, size = 1.5)+
  geom_label_repel(aes(label = ifelse(feature %in% selecFea, feature, ""),
                       col = ifelse(logFC > 0, "up", "dn"),
                       fontface = ifelse(feature %in% tarFea, 2, 1)),
                   alpha = ifelse(dePoly$feature %in% tarFea, 1, 0.5),
                   max.overlaps = Inf, 
                   show.legend = F) +
  labs(x = bquote(bold(log^{fold~change})),
       y = bquote(bold(Significance~(-log^{adjuste~italic(p)~value})))) +
  scale_y_continuous(trans = scales::log10_trans())+
  scale_color_manual(name = "Gene staus",
                     breaks = c("up", "dn", "none-sig"),
                     labels = c("Significantly up-regulated genes",
                                "Significantly down-regulated genes", 
                                "none-sig"),
                     limits = c("up","dn"),
                     values=c("red", "blue", "gray50"))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "Arial", size = 10, face = 2))

ggsave("3.Output/FormalFigures/volcano_PolyClonal_DNvsComboWT.png",  
       width = 6, height = 4, units = "in", dpi = 300, bg = "white")
#--
datPmel <- subset(dat, 
                  genotype %in% c("Pdcd1L/L Lag3L/L-yfp (WT)", 
                                  "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP") &
                    sample == "PML CD8")

Idents(datPmel) <- "genotype"

dePmel <- wilcoxauc(datPmel) %>% 
  filter(group == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP",
         padj < 0.05)


ggplot(dePmel, aes(x = logFC, y = -log10(padj))) +
  geom_point(aes(col = ifelse(logFC > 0, "up", "dn")),
             alpha = 0.5, size = 1.5)+
  geom_label_repel(aes(label = ifelse(feature %in% selecFea, feature, ""),
                       col = ifelse(logFC > 0, "up", "dn"),
                       fontface = ifelse(feature %in% tarFea, 2, 1)),
                   alpha = ifelse(dePmel$feature %in% tarFea, 1, 0.5),
                   max.overlaps = Inf, show.legend = F) +
  labs(x = bquote(bold(log^{fold~change})),
       y = bquote(bold(Significance~(-log^{adjuste~italic(p)~value})))) +
  scale_y_continuous(trans = scales::log10_trans())+
  scale_color_manual(name = "Gene staus",
                     breaks = c("up", "dn", "none-sig"),
                     labels = c("Significantly up-regulated genes",
                                "Significantly down-regulated genes", 
                                "none-sig"),
                     limits = c("up","dn"),
                     values=c("red", "blue", "gray50"))+
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(family = "Arial", size = 10, face = 2))

ggsave("3.Output/FormalFigures/volcano_PMEL_DNvsWT.png",  
       width = 6, height = 4, units = "in", dpi = 300, bg = "white")

# sub clusters of 0 and 7 ----
# subClst <- read.csv("3.Output/3.SubClst/clstrsforsubClst_0_7.csv")
# 
# datSub <- subset(dat, cells = subClst$X)
# 
# datSub <- AddMetaData(datSub, 
#                       subClst %>% 
#                         column_to_rownames(var = "X") %>% 
#                         rename(subClst = 1))
# 
# DimPlot(datSub, group.by = "subClst")
# 
# datSub$newClst <- ifelse(datSub$subClst %in% c(0, 3), "poly", "pmel")
# 
# datSubSub <- subset(datSub, genotype == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP")
# 
# de <- wilcoxauc(datSubSub, group_by = "newClst")
# 
# de <- de %>% filter(group == "poly", logFC != 0, padj != 1)
# 
# # tarFea <- strsplit("Ccl5, Gzmk, Ifit3, Cxcr3, Ifngr1", split = ", ")[[1]]
# topFea <- c(de %>% arrange(desc(logFC)) %>% head(5) %>% dplyr::select(feature) %>% deframe(),
#             de %>% arrange(desc(logFC)) %>% tail(5) %>% dplyr::select(feature) %>% deframe(),
#             de %>% arrange(desc(padj)) %>% tail(5) %>% dplyr::select(feature) %>% deframe())
# 
# selecFea <- unique(c(# tarFea, 
#                      topFea))
# 
# ggplot(de, aes(logFC, -log10(padj))) +
#   geom_point(aes(col = ifelse(padj > 0.05, "none-sig", 
#                               ifelse(logFC > 0, "up", "dn"))),
#              alpha = 0.3, size = 1.5)+
#   geom_label_repel(aes(label = ifelse(feature %in% selecFea, feature, ""),
#                        col = ifelse(padj > 0.05, "none-sig", 
#                                     ifelse(logFC > 0, "up", "dn")),
#                        fontface = 2 #ifelse(feature %in% tarFea, 2, 1)
#                        ),
#                    # alpha = ifelse(de$feature %in% tarFea, 1, 0.5),
#                    max.overlaps = 1000, show.legend = F) +
#   labs(x = bquote(bold(log^{fold~change})),
#        y = bquote(bold(Significance~(-log^{adjuste~italic(p)~value})))) +
#   scale_y_continuous(trans = scales::log10_trans())+
#   scale_color_manual(name = "Gene staus",
#                      breaks = c("up", "dn", "none-sig"),
#                      labels = c("Significantly up-regulated genes in polyclonal cells",
#                                 "Significantly up-regulated genes in pMEL cells", 
#                                 "none-sig"),
#                      limits = c("up","dn"),
#                      values=c("red", "blue", "gray50"))+
#   guides(color = guide_legend(override.aes = list(size=3), nrow = 2))+
#   theme_bw() +
#   theme(legend.position = "bottom",
#         text = element_text(family = "Arial", size = 10, face = 2))
# ggsave("3.Output/FormalFigures/volcano_SubClustr0&7_1.png",  
#        width = 6, height = 6, units = "in", dpi = 300, bg = "white")


