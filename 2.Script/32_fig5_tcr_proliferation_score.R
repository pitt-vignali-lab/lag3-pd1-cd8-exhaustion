# ----------------------------------------------------------------------------
# Fig 5 | Proliferation and cytotoxicity module scores in clonal cells
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/4.5.TCRcellProliferationEvaluation.R
# Input    : 1.Data/datForRevision.RDS
#            1.Data/tcrCombinedUpdated.csv
# Output   : 3.Output/4.TCR/cellProliferationToxcicityEvaluation*.png
# Upstream : 29_fig4_tcr_combined.R
# Note     : Seurat AddModuleScore on the GO terms for positive regulation of CD8 alpha-beta T cell proliferation and of T cell mediated cytotoxicity against tumour targets
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/4.TCR")){
  dir.create("3.Output/4.TCR")
}


library(Seurat)
library(tidyverse)
library(scales)
library(MASS)


dat <- readRDS("1.Data/datForRevision.RDS")

newClst <- read.csv("1.Data/TCR_defined_clusters.csv") %>%
  column_to_rownames(var = "X") %>%
  rename(newClst = dat.newClst)

dat <- AddMetaData(dat, newClst)

Idents(dat) <- "newClst"

datNor <- subset(dat, sample == "Normal CD8")

tcrDf <- read.csv("1.Data/tcrCombinedUpdated.csv")[-1] %>%
  # mutate(barcode = paste0("nor_", barcode)) %>% 
  filter(barcode %in% rownames(datNor@meta.data)) 

tcrTabDf <- data.frame(table(tcrDf$cdr3)) %>% 
  arrange(desc(Freq)) %>% 
  rename(cdr3 = Var1)

# positive regulation of CD8-positive, alpha-beta T cell proliferation: Cd244a, Ptpn22
proFeature <- list(c("Cd244a", "Ptpn22", "Mr1", "Slc22a13"))

datTcr <- AddModuleScore(subset(datNor, cells = Cells(datNor)[Cells(datNor) %in% tcrDf$barcode]),
                         assay = "RNA",
                         features = proFeature,
                         ctrl = 100,
                         name = "CD8proliferating")

proScor <- data.frame(datTcr$CD8proliferating1) %>% 
  rename(score = 1)

proScorSub <- proScor %>% filter(score > 0) 

datNor$genotype <- factor(datNor$genotype, 
                       levels = c("E8ICre.GFP (WT)",
                                  "Pdcd1L/L Lag3L/L-yfp (WT)",
                                  "Lag3L/L-yfp E8ICre.GFP",
                                  "Pdcd1L/L E8ICre.GFP",
                                  "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))

# figure generation ----
DimPlot(datNor, cells.highlight = rownames(proScorSub)) + 
  coord_fixed()
ggsave("3.Output/4.TCR/cellProliferationToxcicityEvaluationAll.png", 
       height = 5, width = 5, units = "in", dpi = 300, bg = "white")

DimPlot(datNor, cells.highlight = rownames(proScorSub), split.by = "genotype") + 
  coord_fixed()
ggsave("3.Output/4.TCR/cellProliferationToxcicityEvaluation.png", 
       height = 4, width = 15, units = "in", dpi = 300, bg = "white")

# for manual
plotDf <- data.frame(Embeddings(datNor, reduction = "umap"),
                     ident = paste0('C', datNor$seurat_clusters),
                     type = datNor$types)

# plotDf$ident <- paste0("C", plotDf$ident)

plotDf$type <- factor(plotDf$type, 
                          levels = c("E8ICre.GFP (WT)",
                                     "Pdcd1L/L Lag3L/L-yfp (WT)",
                                     "Lag3L/L-yfp E8ICre.GFP",
                                     "Pdcd1L/L E8ICre.GFP",
                                     "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))

ggplot(plotDf[!rownames(plotDf) %in% rownames(proScorSub),], 
       aes(UMAP_1, UMAP_2)) +
  geom_point(col = "gray", size = 2, alpha = 0.6) +
  geom_point(data = plotDf[rownames(plotDf) %in% rownames(proScorSub),],
             aes(fill = ident), shape = 21, col = "gray20", size = 2.1,
             alpha = 0.8) +
  coord_fixed() +
  theme_void()+
  guides(fill = guide_legend(title = "Activated cells by clusters", 
                             ncol=2,
                             byrow=TRUE,
                             override.aes = list(size=3))) +
  theme(# legend.position = c(0.6, 0.1),
        # legend.direction = "vertical",
        legend.title.align = 0.5,
        legend.key=element_rect(fill = "gray90", color = NA),
        # legend.box.background = element_rect(colour = "black"),
        legend.text = element_text(size = 11.5, 
                                   face = 2),
        legend.title = element_text(size = 12.5, 
                                    face = 2))
ggsave("3.Output/CellFigure/UMAPcellProliferationToxcicityEvaluationAll.png", 
       height = 5, width = 6, units = "in", dpi = 300, bg = "white")

ggplot(plotDf[!rownames(plotDf) %in% rownames(proScorSub),], 
       aes(UMAP_1, UMAP_2)) +
  geom_point(col = "gray", size = 2, alpha = 1) +
  geom_point(data = plotDf[rownames(plotDf) %in% rownames(proScorSub),],
             aes(fill = ident), shape = 21, col = "gray20", size = 2.1,
             alpha = 0.8) +
  facet_grid(~ type) +
  coord_fixed() +
  theme_void()+
  theme(legend.position = "none",
        strip.placement = "outside",
        # text = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "cm")))
ggsave('3.Output/CellFigure/UMAPcellProliferationToxcicityEvaluationSplit.png',
       height = 5, width = 16, units = "in")

ggsave("3.Output/4.TCR/cellProliferationToxcicityEvaluationSplit(Update).png", 
       height = 5, width = 16, units = "in", dpi = 300, bg = "white")



scoreDf <- data.frame()
for (i in 1:200) {
  t <- data.frame(score = proScor[tcrDf[tcrDf$cmbdCdr3 == tcrTabDf$Var1[i],]$barcode, ])
  t$seq = tcrTabDf$Var1[i]
  scoreDf <- rbind(scoreDf, t)
  
}

scoreDf %>% 
ggplot(aes( y = seq, x = score))+
  geom_point() +
  geom_boxplot()

# log TCR abundance vs CD8 activity ----
actOverall <- data.frame()

corDf <- data.frame()

for (samp in unique(datNor$samples)) {
  datSub <- subset(datNor, samples %in% c(samp))
  
  tcrDfModi <- merge(tcrDf, 
                     tcrDf %>% 
                       filter(barcode %in% Cells(datSub)) %>% 
                       dplyr::select(cdr3) %>% 
                       deframe() %>% table() %>% 
                       data.frame() %>% rename("Var" = 1),
                     by.x = "cdr3", by.y = "Var")
  
  
  
  tcrDfModi$FreqN <- ""
  
  for (i in sort(unique(tcrDfModi$Freq))) {
    bars <- tcrDfModi %>% 
      filter(Freq == i) %>% 
      dplyr::select(barcode) %>% 
      deframe()
    
    tcrDfModi[tcrDfModi$barcode %in% bars,]$FreqN <- i
  }
  
  datNew <- AddMetaData(datSub, 
                        unique(tcrDfModi[c("barcode", "FreqN")]) %>%
                          remove_rownames() %>% 
                          column_to_rownames(var = "barcode"))
  
  datNew@meta.data[is.na(datNew@meta.data$FreqN),] <- 0
  
  datNew <- AddModuleScore(datNew,
                           assay = "RNA",
                           nbin = 50,
                           features = proFeature,
                           ctrl = 100,
                           name = "CD8Act")
  
  act <- datNew@meta.data[c('FreqN', 'CD8Act1')] %>% 
    mutate(FreqN = as.numeric(FreqN)) %>% 
    filter(CD8Act1 > 0, FreqN > 0 )
  
  act$sample <- samp
 
  corMod <- cor.test(log10(act$FreqN), act$CD8Act1)
  
  corDfSub <- data.frame(sample = samp,
                         cor = corMod[["estimate"]][["cor"]],
                         corP = corMod[["p.value"]])
  corDf <- rbind(corDf, corDfSub)
  
  actOverall <- rbind(actOverall, act)
}

# as long as cells are with tcr expression, and positive cd8 acitivity, the result is interesting

corDf

actOverall %>% 
  mutate(sample = factor(sample, 
                         levels = c("Pdcd1L/L Lag3L/L-yfp (WT) - Tumor",
                                    "E8ICre.GFP (WT) - Tumor",
                                    "Lag3L/L-yfp E8ICre.GFP - Tumor",
                                    "Pdcd1L/L E8ICre.GFP - Tumor",
                                    "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP - Tumor"))) %>% 
  ggplot(aes(x = CD8Act1, y= log10(FreqN), 
             col = sample)) +
  geom_point(size = 2, alpha = 0.15) +
  geom_smooth(aes(group = sample, fill = sample), 
              alpha = 0.1,
              method = lm,
              show.legend = F,
              se = T) +
  scale_color_manual(values=c("gray50", "black", "blue", "red", "purple")) +
  scale_fill_manual(values=c("gray50", "black", "blue", "red", "purple")) +
  guides(color = guide_legend(nrow=2,byrow=TRUE)) +
  labs(y = "Logarithm abundance of TCR",
       x = "CD8 activity") +
  coord_cartesian(ylim=c(0, 2)) +
  theme(
    text = element_text(size = 16, face = 2),
    legend.position = "none",
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, color = "grey75"),
    axis.ticks = element_line(color = "grey85"),
    panel.grid.major = element_line(color = "grey95", size = 0.2),
    panel.grid.minor = element_line(color = "grey95", size = 0.2),
    legend.key = element_blank())

ggsave("3.Output/4.TCR/correlate_CD8act_combined.png",
       width = 5, height = 4, units = "in")

# UCell test ----
# proFeature <- list(CD8Act1 = c("Cd244a", "Ptpn22", "Mr1", "Slc22a13"))
# 
# actOverall <- data.frame()
# 
# corDf <- data.frame()
# 
# for (samp in unique(datNor$samples)) {
#   datSub <- subset(datNor, samples %in% c(samp))
#   
#   tcrDfModi <- merge(tcrDf, 
#                      tcrDf %>% 
#                        filter(barcode %in% Cells(datSub)) %>% 
#                        dplyr::select(cdr3) %>% 
#                        deframe() %>% table() %>% 
#                        data.frame() %>% rename("Var" = 1),
#                      by.x = "cdr3", by.y = "Var")
#   
#   
#   
#   tcrDfModi$FreqN <- ""
#   
#   for (i in sort(unique(tcrDfModi$Freq))) {
#     bars <- tcrDfModi %>% 
#       filter(Freq == i) %>% 
#       dplyr::select(barcode) %>% 
#       deframe()
#     
#     tcrDfModi[tcrDfModi$barcode %in% bars,]$FreqN <- i
#   }
#   
#   datNew <- AddMetaData(datSub,
#                         unique(tcrDfModi[c("barcode", "FreqN")]) %>%
#                           remove_rownames() %>%
#                           column_to_rownames(var = "barcode"))
#   
#   datNew@meta.data[is.na(datNew@meta.data$FreqN),] <- 0
#   
#   datNew <- AddModuleScore_UCell(datNew, 
#                                  assay = "RNA",
#                                  features = proFeature)
#   
#   act <- datNew@meta.data[c('FreqN', 'CD8Act1_UCell')] %>% 
#     mutate(FreqN = as.numeric(FreqN)) %>% 
#     filter(FreqN > 0)
#   
#   act$sample <- samp
#   
#   corMod <- cor.test(log10(act$FreqN+1), act$CD8Act1_UCell)
#   
#   corDfSub <- data.frame(sample = samp,
#                          cor = corMod[["estimate"]][["cor"]],
#                          corP = corMod[["p.value"]])
#   corDf <- rbind(corDf, corDfSub)
#   
#   actOverall <- rbind(actOverall, act)
# }
# 
# # as long as cells are with tcr expression, and positive cd8 acitivity, the result is interesting
# 
# corDf
# 
# actOverall %>% 
#   mutate(sample = factor(sample, 
#                          levels = c("Pdcd1L/L Lag3L/L-yfp (WT) - Tumor",
#                                     "E8ICre.GFP (WT) - Tumor",
#                                     "Lag3L/L-yfp E8ICre.GFP - Tumor",
#                                     "Pdcd1L/L E8ICre.GFP - Tumor",
#                                     "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP - Tumor"))) %>% 
#   ggplot(aes(x = CD8Act1_UCell, y= log10(FreqN+1), 
#              col = sample)) +
#   geom_point(size = 2, alpha = 0.15) +
#   geom_smooth(aes(group = sample, fill = sample), 
#               alpha = 0.1,
#               method = lm,
#               show.legend = F,
#               se = F) +
#   scale_color_manual(values=c("gray50", "black", "blue", "red", "purple")) +
#   scale_fill_manual(values=c("gray50", "black", "blue", "red", "purple")) +
#   guides(color = guide_legend(nrow=2,byrow=TRUE)) +
#   labs(y = "Logarithm abundance of TCR",
#        x = "CD8 activity") +
#   theme(
#     text = element_text(size = 16, face = 2),
#     legend.position = "none",
#     panel.background = element_rect(fill = NA),
#     panel.border = element_rect(fill = NA, color = "grey75"),
#     axis.ticks = element_line(color = "grey85"),
#     panel.grid.major = element_line(color = "grey95", size = 0.2),
#     panel.grid.minor = element_line(color = "grey95", size = 0.2),
#     legend.key = element_blank())
# 
# ggsave("3.Output/4.TCR/correlate_CD8act_combined.png",
#        width = 5, height = 5, units = "in")
