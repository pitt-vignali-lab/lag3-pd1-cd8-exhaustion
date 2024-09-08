# ----------------------------------------------------------------------------
# Fig 4 | Combined TCR table across genotypes and samples
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/4.2.TCRComboUpdated.R
# Input    : 1.Data/TcrSplited/filtered_contig_annotations.csv
#            1.Data/TcrSplited/TRUST_merged_r2_report_modified.csv
# Output   : 1.Data/tcrCombinedUpdated.csv
#            3.Output/CellFigure/cmbdTCR_distribution_org(updated).png
# Upstream : 26_fig4_tcr_init.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)

if(!dir.exists("3.Output/4.TCR")){
  dir.create("3.Output/4.TCR")
}

dat <- readRDS("1.Data/datForRevision.RDS")

newClst <- read.csv("1.Data/TCR_defined_clusters.csv") %>%
  column_to_rownames(var = "X") %>%
  rename(newClst = dat.newClst)

dat <- AddMetaData(dat, newClst)

datNor <- subset(dat, sample == "Normal CD8")

# TCR information ognization ----
tcrSeq <- read.csv("1.Data/TcrSplited/filtered_contig_annotations.csv") %>% 
  filter(chain == "TRB",
         cdr3 != "None",
         productive == 'true',
         full_length == "true") %>% 
  group_by(barcode) %>% 
  top_n(1, reads) %>% 
  top_n(1, umis) %>% 
  select('barcode', 'cdr3', 'reads') %>% 
  mutate(barcode = gsub('-[0-9]+', '', barcode)) %>% 
  unique()

# test duplicated cells
tcrSeqSum <- tcrSeq %>% group_by(barcode, cdr3) %>% summarise(n = n())
table(tcrSeqSum$n)

tcrTrs <- read.table("1.Data/Trust4Output/TRUST_merged_r2_report.tsv", 
                     sep = "\t",
                     header = T) 

tcrTrs$chain <- apply(tcrTrs[c('V', 'D', 'J', 'C')],
               1,
               function(x) sort(unique(substr(x, 1, 3)),
                                decreasing = T)[1])
tcrTrs <- tcrTrs %>% 
  filter(CDR3aa != 'out_of_frame',
         chain == 'TRB') %>%
  mutate(barcode = gsub('_[0-9]+', '', cid)) %>% 
  group_by(barcode) %>% 
  top_n(1, count) %>% 
  select(barcode, CDR3aa, count) %>% 
  rename(cdr3 = 2, reads = 3)


# combine two datasets
tcrDf <- rbind(tcrSeq, tcrTrs) %>% 
  group_by(barcode) %>% 
  top_n(1, reads)

# check redundant cells with multiple TCR-beta reads and remove them 
tcrSum = tcrDf %>% group_by(barcode) %>% summarise(n = n())

tcrReddt <- tcrSum %>% filter(n > 1) %>% 
  pull(barcode)

tcrDf <- tcrDf %>% filter(!barcode %in% tcrReddt ) %>% 
  select(barcode, cdr3)

# checking if there is duplicated barcode
tcrDf %>% group_by(barcode) %>% summarise(n = n()) %>% pull(n) %>% table()

# combine sc with tcr ----
head(Cells(datNor))

tcrDf <- tcrDf %>% 
  mutate(barcode = paste('nor', barcode, sep = '_')) %>% 
  filter(barcode %in% Cells(datNor))

write.csv(tcrDf, "1.Data/tcrCombinedUpdated.csv")

# this combined TCR is from tcr seq(without modification) and tcr TRUST directly 
tcrDf <- read.csv("1.Data/tcrCombinedUpdated.csv")[-1] %>%
  # mutate(barcode = paste0("nor_", barcode)) %>% 
  filter(barcode %in% rownames(datNor@meta.data)) %>% 
  rename(cmbdCdr3 = cdr3) %>% 
  unique

# TCR frequency summary of combined file ----
barplot(table(tcrDf$cmbdCdr3))
cdr3Df <- data.frame(table(tcrDf$cmbdCdr3))

ggplot(cdr3Df %>% arrange(desc(Freq)) %>% head(15), 
       aes(y = reorder(Var1, Freq), x = Freq)) +
  geom_bar(stat = "identity") +
  labs(y = "Clonaltype", x = "Freqence" ) +
  theme_minimal()+
  labs(x = '', y = ' ') +
  theme(text = element_text(size = 10, face = 2),
        axis.text.x = element_blank())

ggsave("3.Output/CellFigure/cmbdTCR_clonaltypes_freqence(Updated).png", width = 6, height = 3,
       units = "in", dpi = 300, bg = "white")

ggsave("3.Output/4.TCR/cmbdTCR_clonaltypes_freqence(Updated).png", width = 6, height = 3,
       units = "in", dpi = 300, bg = "white")

# Tcr distribution ----
targetTcr <- cdr3Df[order(cdr3Df$Freq, decreasing = T),] %>% 
  arrange(desc(Freq)) %>% head(5)

plotDf <- Reduce(function(x, y) merge(x, y , by = "cells"),
                 list(data.frame(Embeddings(datNor, reduction = "umap")) %>% 
                        rownames_to_column(var = "cells"),
                      data.frame(Idents(datNor))%>% 
                        rownames_to_column(var = "cells"),
                      data.frame(datNor$newClst)%>% 
                        rownames_to_column(var = "cells"))) %>% 
  mutate(Idents.datNor. = factor(as.character(Idents.datNor.), levels = 0:9),
         datNor.newClst = factor(as.character(datNor.newClst), levels = 0:10))

plotDf <- unique(merge(plotDf, tcrDf, by.x = "cells", by.y = "barcode", 
                       all.x = T)) %>% 
  unique

labDf <- plotDf %>% select("UMAP_1", "UMAP_2", "Idents.datNor.") %>% 
  group_by(Idents.datNor.) %>% 
  summarise(UMAP_1 = mean(UMAP_1), 
            UMAP_2 = mean(UMAP_2))

# labDf <- plotDf %>% select("UMAP_1", "UMAP_2", "Idents.datNor.") %>% 
#   group_by(Idents.datNor.) %>% 
#   summarise(UMAP_1 = mean(UMAP_1), 
#             UMAP_2 = mean(UMAP_2))

ggplot(data= plotDf[c(1:5)], aes(UMAP_1, UMAP_2, fill = Idents.datNor.)) +
  geom_point(shape = 21, alpha = 0.1, size = 5, col = "gray50") +
  geom_point(data = plotDf[plotDf$cmbdCdr3 %in% targetTcr$Var1,] %>% 
               mutate(cmbdCdr3 = factor(cmbdCdr3, levels = targetTcr$Var1)),
             aes(UMAP_1, UMAP_2), col = "#231f20", size = 2, alpha = 0.5) +
  geom_text(data = labDf,
            aes(UMAP_1, UMAP_2, label = Idents.datNor.),
            size = 8, fontface = 2) +
  facet_wrap(~cmbdCdr3, ncol = 5) +
  coord_fixed(ratio = 1) +
  # scale_color_brewer(palette="Set2") +
  guides(fill = "none") +
  theme_minimal() +
  labs(x = '', y = '') +
  theme(legend.position = "none",
        text = element_text(size = 20, face = 2),
        axis.text = element_blank())
ggsave("3.Output/CellFigure/cmbdTCR_distribution_org(updated).png", width = 15, height = 5, units = "in",
       dpi = 300, bg = "white")

ggsave("3.Output/4.TCR/cmbdTCR_distribution_org(updated).png", width = 15, height = 5, units = "in",
       dpi = 300, bg = "white")

# Tcr genotype wise ----
plotDf <- Reduce(function(x, y) merge(x, y , by = "cells"),
                 list(data.frame(Embeddings(datNor, reduction = "umap")) %>%
                        rownames_to_column(var = "cells"),
                      data.frame(Idents(datNor))%>%
                        rownames_to_column(var = "cells"),
                      data.frame(datNor$genotype)%>%
                        rownames_to_column(var = "cells")))


plotDf <- unique(merge(plotDf, tcrDf, by.x = "cells", by.y = "barcode", all = T))


tcr <- merge(tcrDf, plotDf %>% select(cells, datNor.genotype), 
             by.x = "barcode", by.y = "cells", all = T) %>% 
  unique %>% 
  rename(genotype = datNor.genotype)

cells <- "Pdcd1L/L E8ICre.GFP"

overall <- data.frame(Embeddings(dat, reduction = "umap"))
for (cells in unique(datNor$genotype)) {
  
  targetTcr <- tcr %>% filter(genotype == cells) %>%
    na.omit() %>% 
    group_by(cmbdCdr3) %>% 
    summarise(Freq = n()) %>% 
    arrange(desc(Freq)) %>%
    # select(cmbdCdr3) %>%
    head(10) # %>% tibble::deframe()
  
  targetTcr$cmbdCdr3 <- factor(targetTcr$cmbdCdr3, 
                               levels = targetTcr[order(targetTcr$Freq, decreasing = T),]$cmbdCdr3)
  
  ggplot(targetTcr,
         aes(x = Freq, y = cmbdCdr3, fill = cmbdCdr3)) +
    geom_bar(stat = "identity") +
    geom_label(aes(label = Freq), fill = "white") +
    scale_y_discrete(limits=rev) +
    xlim(0, 300) +
    theme_minimal() +
    labs(x = '') +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          text = element_text(size = 15, face = 2))
  
  ggsave(paste0("3.Output/CellFigure/TCRSpcific_", gsub("/", "_", cells), "_cmbdTCR_bar(Updtd).png"), 
         width = 8, height = length(targetTcr$cmbdCdr3)/3, units = "in", 
         dpi = 300, bg = "white")
  
  
  ggplot() +
    geom_point(data = overall, aes(x = UMAP_1, y = UMAP_2), 
               col = "gray90", alpha = 0) +
    geom_point(data = plotDf[plotDf$datNor.genotype == cells,][-6] %>% 
                 mutate(datNor.newClst = as.character(Idents.datNor.)), 
               aes(UMAP_1, UMAP_2), fill = "gray10", 
               shape = 21, alpha = 0.3, size = 5) +
    geom_point(data = plotDf[plotDf$datNor.genotype == cells,] %>% 
                 filter(cmbdCdr3 %in% targetTcr$cmbdCdr3[1]) %>%  # remove [1], it will show top 10 tcr
                 mutate(cmbdCdr3 = factor(cmbdCdr3, levels = targetTcr$cmbdCdr3)),
               aes(UMAP_1, UMAP_2, color = cmbdCdr3), size = 3, alpha = 0.8)+
    facet_wrap( ~ cmbdCdr3, ncol = 5) +
    coord_fixed(ratio = 1) +
    # scale_color_brewer(breaks = targetTcr) +
    guides(fill = "none") +
    theme_minimal() +
    theme(legend.position = "none",
          title = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = 20, face = 2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.2))
  
  ggsave(paste0("3.Output/CellFigure/TCRSpcific_", gsub("/", "_", cells), "_cmbdTCR_distr(Updtd).png"), width = 3, 
         height = 4, units = "in", 
         bg = "white")
  
  # ggsave(paste0("3.Output/4.TCR/", gsub("/", "_", cells), "_cmbdTCR_distr.png"), width = 9, 
  #        height = ifelse(length(targetTcr$cmbdCdr3) >5, 5, 2.5), units = "in", 
  #        bg = "white")
  
}




