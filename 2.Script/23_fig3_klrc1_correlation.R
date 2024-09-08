# ----------------------------------------------------------------------------
# Fig 3 | Correlation of clonal features with Klrc1 expression
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/4.6.corAnalysisVsKlrc1Exp.R
# Input    : 1.Data/datForRevision.RDS
#            1.Data/tcrCombinedUpdated.csv
# Output   : 3.Output/4.TCR/correlate_*.png
#            3.Output/CellFigure/correlate_hiVslow_DE_*.png
# Upstream : 29_fig4_tcr_combined.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(scales)
library(MASS)
library(presto)
library(scMerge)

if(!dir.exists("3.Output/4.TCR")){
  dir.create("3.Output/4.TCR")
}

dat <- readRDS("1.Data/datForRevision.RDS")

table(dat$sample)

DimPlot(subset(dat, sample == "Normal CD8"))

datNor <- subset(dat, sample == "Normal CD8")

tcrDf <- read.csv("1.Data/tcrCombinedUpdated.csv")[-1] %>%
  # mutate(barcode = paste0("nor_", barcode)) %>% 
  filter(barcode %in% rownames(datNor@meta.data)) 

colPanel = data.frame(sample = unique(datNor$samples),
                      col = c("black", "purple", "red", "blue", "gray50"))

deOverall <- data.frame()
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
  
  datNew <- subset(datNew, FreqN != 0)
  
  de <- wilcoxauc(datNew, group_by = "FreqN") %>% 
    filter(feature == "Klrc1") 
  
  de$group <- as.numeric(de$group)
  
  # linMod <- lm(group~logFC, data = de)  
  # linSum <- summary(linMod)
  
  corMod <- cor.test(de$logFC, log10(de$group))
  
  de$cor <- corMod[["estimate"]][["cor"]]
  de$corP <- corMod[["p.value"]]
  de$sample <- samp
  
  deOverall <- rbind(deOverall, de)
}

corDf <- deOverall %>% 
  mutate(sample = gsub("- Tumor", "", sample)) %>% 
  dplyr::select('sample', 'cor', 'corP') %>% 
  unique()
corDf$title <- paste0(c("Cor[E8I^{Cre.GFP}]:~",
                        "Cor[Pdcd1^{L/L}~Lag3^{L/L-YFP}~E8I^{Cre.GFP}]:~",
                        "Cor[Pdcd1^{L/L-YFP}~E8I^Cre.GFP]:~",
                        "Cor[Lag3^{L/L-YFP}~E8I^Cre.GFP]:~",
                        "Cor[Pdcd1^{L/L}~Lag3^{L/L-YFP}]:~"),
                      format(corDf$cor, digits = 3),
                      rep("~italic((p-value:~", 5),
                      format(corDf$corP, digits = 4, scientific = T),
                      rep("))", 5))

deOverall %>% 
  mutate(sample = factor(sample, 
                        levels = c("Pdcd1L/L Lag3L/L-yfp (WT) - Tumor",
                                   "E8ICre.GFP (WT) - Tumor",
                                   "Lag3L/L-yfp E8ICre.GFP - Tumor",
                                   "Pdcd1L/L E8ICre.GFP - Tumor",
                                   "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP - Tumor"))) %>% 
  ggplot(aes(x= logFC, y = log10(group), col = sample)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(aes(group = sample, fill = sample), 
              alpha = 0.1,
              method = lm,
              show.legend = F,
              se = T) +
  scale_color_manual(values=c("gray50", "black", "blue", "red", "purple")) +
  scale_fill_manual(values=c("gray50", "black", "blue", "red", "purple")) +
  coord_fixed(ratio = 0.5) +
  # annotate(geom = "text",
  #          label = corDf$title[c(1, 5, 4, 3, 2)],
  #          col = c("gray50", "black", "blue", "red", "purple"),
  #          vjust = 3:7, hjust = 0,
  #          x = -Inf + 1, y = Inf + 1, parse = T) +
  guides(color = guide_legend(nrow=2,byrow=TRUE)) +
  labs(x = "Logarithm fold change of NKG2A",
       y = "Logarithm abundance of TCR") +
  theme(
    text = element_text(size = 16, face = 2),
    legend.position = "none",
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, color = "grey75"),
    axis.ticks = element_line(color = "grey85"),
    panel.grid.major = element_line(color = "grey95", size = 0.2),
    panel.grid.minor = element_line(color = "grey95", size = 0.2),
    legend.key = element_blank())

ggsave("3.Output/4.TCR/correlate_DE_combinedN.png",
       width = 6, height = 6, units = "in")

# method 3 (set threshold to evaluate DE level of NKG2A) ----
plotDf <- data.frame()

for (samp in c("all", unique(datNor$samples))) {
  
  if(samp == "all"){
    datSub <- datNor
  } else {
    datSub <- subset(datNor, samples == samp )
  }
  
  tcrSub <- tcrDf %>% filter(barcode %in% Cells(datSub))
  
  tcrSum <- data.frame(table(tcrSub$cdr3)) %>% 
    arrange(-Freq)
  
  # set threshold with percentage ----
  for(i in 1:length(tcrSum$Freq)){
    if(sum(tcrSum$Freq[1:i]) > sum(tcrSum$Freq) * 0.10) {
      break
      print(i)}
  }
  # sum(tcrSum$Freq[1:i])

  tcrSum$Freq[i]

  plot(log10(tcrSum$Freq))
  abline(h = log10(tcrSum$Freq[i]), col = "red")
  
  # hiAb <- tcrSum %>% 
  #   filter(Freq >= tcrSum$Freq[i]) %>% 
  #   dplyr::select(Var1) %>% 
  #   deframe() %>% 
  #   unique()
  # 
  # hiAbBar <- tcrDf %>% 
  #   filter(cdr3 %in% hiAb) %>% 
  #   dplyr::select("barcode") %>% 
  #   deframe
  
  # set threshold with solid number just change tcrSum$Freq[i] to a number 
  hiAb <- tcrSum %>%
    filter(Freq >= tcrSum$Freq[i]) %>%
    dplyr::select(Var1) %>%
    deframe() %>%
    unique()

  hiAbBar <- tcrDf %>%
    filter(cdr3 %in% hiAb) %>%
    dplyr::select("barcode") %>%
    deframe

  tcrDf$cat <- "Low"
  
  tcrDf[tcrDf$barcode %in% hiAbBar,]$cat<- "High"
  
  datNew <- AddMetaData(datNor, 
                        tcrDf %>% 
                          dplyr::select("barcode", "cat") %>% 
                          unique() %>% 
                          remove_rownames() %>% 
                          column_to_rownames(var = "barcode"))
  
  datNew <- subset(datNew, 
                   cells = tcrDf %>% 
                     dplyr::select("barcode") %>% 
                     deframe() %>% 
                     unique())
  deDf = wilcoxauc(datNew, group_by = "cat") %>% 
    filter(feature == "Klrc1")
  
  plotDfSub <- VlnPlot(datNew, features = "Klrc1", group.by = "cat")$data
  plotDfSub$samples <- paste0(samp, "\npdaj = ", 
                              format(deDf$pval[1], digits = 3, scientific = T))
  
  plotDf <- rbind(plotDf, plotDfSub)
}

plotDf %>% 
  mutate(samples = gsub(" - Tumor", "", samples)) %>% 
  mutate(samples = factor(samples,
                          levels = unique(samples)[c(1, 2, 6, 5, 4, 3)])) %>%  
ggplot(aes(x = ident, y = Klrc1, fill = ident))+
  geom_violin(alpha = 0.6) +
  geom_jitter(alpha = 0.5, size = 0.1)+
  geom_boxplot(width=0.1, fill = "white") +
  facet_wrap(~samples, nrow = 1)+
  labs(x = "",
       y = "Gene expression level of NKG2A") +
  coord_fixed(ratio = 0.8) +
  theme(
    text = element_text(size = 15, face = 2),
    legend.position = "none",
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, color = "grey50"),
    axis.ticks = element_line(color = "grey85"),
    strip.background = element_rect(fill=" white"),
    strip.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey95", size = 0.2),
    panel.grid.minor = element_line(color = "grey95", size = 0.2),
    legend.key = element_blank(),
    strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "cm")))

ggsave("3.Output/CellFigure/correlate_hiVslow_DE_0.10_Labeled.png",
       width = 15, height = 5, units = "in")

plotDf %>% 
  mutate(samples = gsub(" - Tumor", "", samples)) %>% 
  mutate(samples = gsub("\n.*", "", samples)) %>% 
  mutate(samples = factor(samples,
                          levels = unique(samples)[c(1, 2, 6, 5, 4, 3)])) %>%  
  ggplot(aes(x = ident, y = Klrc1, fill = ident))+
  geom_violin(alpha = 0.6) +
  geom_jitter(alpha = 0.5, size = 0.1)+
  geom_boxplot(width=0.1, fill = "white") +
  facet_wrap(~samples, nrow = 1)+
  labs(x = "",
       y = "") +
  coord_fixed(ratio = 0.8) +
  theme(
    axis.text = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, color = "grey50"),
    axis.ticks = element_line(color = "grey85"),
    strip.background = element_rect(fill=" white"),
    strip.text = element_text(size = 10, face = 2),
    panel.grid.major = element_line(color = "grey95", size = 0.2),
    panel.grid.minor = element_line(color = "grey95", size = 0.2),
    legend.key = element_blank(),
    strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "cm")))

ggsave("3.Output/CellFigure/correlate_hiVslow_DE_0.10_NoLabel.png",
       width = 15, height = 5, units = "in")
#################################################################################
sce <- as.SingleCellExperiment(datNor)

exprs_mat <- SummarizedExperiment::assay(sce, 'logcounts')

scSEG = scSEGIndex(exprs_mat = exprs_mat)


plotDfSeg <- data.frame()

for (samp in c("all", unique(datNor$samples))) {
  
  if(samp == "all"){
    datSub <- datNor
  } else {
    datSub <- subset(datNor, samples == samp )
  }
  
  tcrSub <- tcrDf %>% filter(barcode %in% Cells(datSub))
  
  tcrSum <- data.frame(table(tcrSub$cdr3)) %>% 
    arrange(-Freq)
  
  # set threshold with percentage ----
  for(i in 1:length(tcrSum$Freq)){
    if(sum(tcrSum$Freq[1:i]) > sum(tcrSum$Freq) * 0.25) {
      break
      print(i)}
  }
  # sum(tcrSum$Freq[1:i])
  
  tcrSum$Freq[i]
  
  plot(log10(tcrSum$Freq))
  abline(h = log10(tcrSum$Freq[i]), col = "red")
  
  # hiAb <- tcrSum %>% 
  #   filter(Freq >= tcrSum$Freq[i]) %>% 
  #   dplyr::select(Var1) %>% 
  #   deframe() %>% 
  #   unique()
  # 
  # hiAbBar <- tcrDf %>% 
  #   filter(cdr3 %in% hiAb) %>% 
  #   dplyr::select("barcode") %>% 
  #   deframe
  
  # set threshold with solid number just change tcrSum$Freq[i] to a number 
  hiAb <- tcrSum %>%
    filter(Freq >= tcrSum$Freq[i]) %>%
    dplyr::select(Var1) %>%
    deframe() %>%
    unique()
  
  hiAbBar <- tcrDf %>%
    filter(cdr3 %in% hiAb) %>%
    dplyr::select("barcode") %>%
    deframe
  
  tcrDf$cat <- "Low"
  
  tcrDf[tcrDf$barcode %in% hiAbBar,]$cat<- "High"
  
  datNew <- AddMetaData(datNor, 
                        tcrDf %>% 
                          dplyr::select("barcode", "cat") %>% 
                          unique() %>% 
                          remove_rownames() %>% 
                          column_to_rownames(var = "barcode"))
  
  datNew <- subset(datNew, 
                   cells = tcrDf %>% 
                     dplyr::select("barcode") %>% 
                     deframe() %>% 
                     unique())
  deDf = wilcoxauc(datNew, group_by = "cat") %>% 
    filter(feature == "Mlana")
  
  plotDfSub <- VlnPlot(datNew, features = "Mlana", group.by = "cat")$data
  plotDfSub$samples <- paste0(samp, "\npdaj = ", 
                              format(deDf$pval[1], digits = 3, scientific = T))
  
  plotDfSeg <- rbind(plotDfSeg, plotDfSub)
}

plotDfSeg %>% 
  mutate(samples = gsub(" - Tumor", "", samples)) %>%
  mutate(samples = factor(samples,
                          levels = c("all\npdaj = 4.07e-01",
                                     "E8ICre.GFP (WT)\npdaj = 8.32e-01",
                                     "Pdcd1L/L Lag3L/L-yfp (WT)\npdaj = 2.57e-01",
                                     "Lag3L/L-yfp E8ICre.GFP\npdaj = 7.14e-01",
                                     "Pdcd1L/L E8ICre.GFP\npdaj = 3.78e-01",
                                     "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP\npdaj = 3.71e-01"))) %>% 
  ggplot(aes(x = ident, y = Mlana, fill = ident))+
  geom_violin(alpha = 0.6) +
  # geom_jitter(alpha = 0.5, size = 0.1)+
  geom_boxplot(width=0.1, fill = "white") +
  facet_wrap(~samples, nrow = 1)+
  labs(x = "",
       y = "Gene expression level of Mlana") +
  coord_fixed(ratio = 0.8) +
  theme(
    text = element_text(size = 15, face = 2),
    legend.position = "none",
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, color = "grey50"),
    axis.ticks = element_line(color = "grey85"),
    strip.background = element_rect(fill=" white"),
    strip.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey95", size = 0.2),
    panel.grid.minor = element_line(color = "grey95", size = 0.2),
    legend.key = element_blank(),
    strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "cm")))
ggsave("3.Output/4.TCR/correlate_hiVslow_DE_hskp.png",
       width = 15, height = 5, units = "in")

print(samp)
print(t)

# method 2 (TCR abundance vs actual NKG2A expression level) ----
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
                          mutate(FreqN = factor(FreqN, 
                                                levels = sort(unique(tcrDfModi$Freq)))) %>% 
                          remove_rownames() %>% 
                          column_to_rownames(var = "barcode"))
  
  
  datNew@meta.data[is.na(datNew@meta.data$FreqN),] <- 0
  
  datNew <- subset(datNew, FreqN != 0)
  
  t <- VlnPlot(datNew, features = "Klrc1", group.by = "FreqN")
  
  linMod <- t$data %>% 
    drop_na() %>% 
    mutate(ident = as.numeric(as.character(ident))) %>%
    lm(ident~Klrc1, data = .)  
  
  corMod <- cor.test(t$data %>% 
                       drop_na() %>% 
                       mutate(ident = as.numeric(as.character(ident))) %>% 
                       dplyr::select(Klrc1) %>% 
                       deframe(), 
                     t$data %>% 
                       drop_na() %>% 
                       mutate(ident = as.numeric(as.character(ident))) %>% 
                       dplyr::select(ident) %>% 
                       deframe())
  
  t$data %>% 
    drop_na() %>%
    filter(ident != 0) %>% 
    mutate(ident = log10(as.numeric(as.character(ident)))) %>% 
    ggplot(aes(Klrc1, ident)) +
    geom_point(size = 3, alpha = 0.5, col = colPanel[colPanel$sample == samp,]$col) +
    geom_smooth(method = lm, col = colPanel[colPanel$sample == samp,]$col) +
    annotate(geom = "text",
             label = paste("Correlation = ",
                           format(corMod[["estimate"]][["cor"]],
                                  digits = 3),
                           "\n",
                           "     p-value = ",
                           format(corMod[["p.value"]],
                                  digits = 4,
                                  scientific = T),
                           sep = ""),
             vjust = 3, hjust = -0.2, 
             x = -Inf, y = Inf) +
    labs(x = "Expression value of NKG2A",
         y = "Logarithm of expression abundance of TCR",
         title = samp) +
    theme(
      text = element_text(size = 10, face = 2),
      legend.position = "bottom",
      panel.background = element_rect(fill = NA),
      panel.border = element_rect(fill = NA, color = "grey75"),
      axis.ticks = element_line(color = "grey85"),
      panel.grid.major = element_line(color = "grey95", size = 0.2),
      panel.grid.minor = element_line(color = "grey95", size = 0.2),
      legend.key = element_blank())
  
  ggsave(paste0("3.Output/4.TCR/correlate_realValue", gsub("/", "_", samp), ".png"),
         width = 5, height = 6, units = "in")
  
  corMod <- cor.test(t$data %>% 
                       drop_na() %>%
                       filter(ident != 0) %>% 
                       mutate(ident = log10(as.numeric(as.character(ident)))) %>% 
                       group_by(ident) %>% 
                       summarise(mean = median(Klrc1)) %>% 
                       dplyr::select(mean) %>% 
                       deframe(), 
                     t$data %>% 
                       drop_na() %>%
                       filter(ident != 0) %>% 
                       mutate(ident = log10(as.numeric(as.character(ident)))) %>% 
                       group_by(ident) %>% 
                       summarise(mean = median(Klrc1)) %>% 
                       dplyr::select(ident) %>% 
                       deframe())
  
  t$data %>% 
    drop_na() %>%
    filter(ident != 0) %>% 
    mutate(ident = log10(as.numeric(as.character(ident)))) %>% 
    group_by(ident) %>% 
    summarise(mean = median(Klrc1)) %>% 
    ggplot(aes(mean, ident)) +
    geom_point(size = 3, alpha = 0.5, col = colPanel[colPanel$sample == samp,]$col) +
    geom_smooth(method = lm, col = colPanel[colPanel$sample == samp,]$col) +
    annotate(geom = "text",
             label = paste("Correlation = ",
                           format(corMod[["estimate"]][["cor"]],
                                  digits = 3),
                           "\n",
                           "     p-value = ",
                           format(corMod[["p.value"]],
                                  digits = 4,
                                  scientific = T),
                           sep = ""),
             vjust = 3, hjust = -0.2, 
             x = -Inf, y = Inf) +
    labs(x = "Expression value of NKG2A",
         y = "Logarithm of expression abundance of TCR",
         title = samp) +
    theme(
      text = element_text(size = 10, face = 2),
      legend.position = "bottom",
      panel.background = element_rect(fill = NA),
      panel.border = element_rect(fill = NA, color = "grey75"),
      axis.ticks = element_line(color = "grey85"),
      panel.grid.major = element_line(color = "grey95", size = 0.2),
      panel.grid.minor = element_line(color = "grey95", size = 0.2),
      legend.key = element_blank())
  
  ggsave(paste0("3.Output/4.TCR/correlate_realValue", gsub("/", "_", samp), "(median).png"),
         width = 5, height = 6, units = "in")
  
}
