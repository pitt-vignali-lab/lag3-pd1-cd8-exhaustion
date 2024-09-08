# ----------------------------------------------------------------------------
# Fig 4 | Clonotype calling from TCRbeta CDR3 sequences
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/4.1.2.TCR.R
# Input    : 1.Data/tcrCombined.csv
# Output   : 3.Output/4.TCR/TCR_clonaltypes_freqence.png
# Upstream : 26_fig4_tcr_init.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(scRepertoire)

if(!dir.exists("3.Output/4.TCR")){
  dir.create("3.Output/4.TCR")
}

dat <- readRDS("1.Data/datInt.RDS")
dat <- AddMetaData(dat, 
                   read.csv("1.Data/TCR_defined_clusters.csv", row.names = 1) %>% 
                     rename(newClst = dat.newClst))

DimPlot(dat, group.by = "newClst", label = T)

datNor <- subset(dat, newClst %in% c(1, 3, 10) & 
                   sample == "Normal CD8")

DimPlot(datNor, group.by = "newClst")

# cell plot ----

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

# TCR plot ----

tcr <- data.frame(table(datNor$genotype, datNor$cdr3))

tcr <- tcr[tcr$Freq != 0,]

names(tcr)[c(1, 2)] <- c("genotype", "cdr3")

tcr$genotype <- factor(tcr$genotype, 
                     levels = c("E8ICre.GFP (WT)",
                                "Pdcd1L/L Lag3L/L-yfp (WT)",
                                "Pdcd1L/L E8ICre.GFP",
                                "Lag3L/L-yfp E8ICre.GFP",
                                "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))
cells <- "Pdcd1L/L E8ICre.GFP"

for (cells in c("Pdcd1L/L E8ICre.GFP",
                "Lag3L/L-yfp E8ICre.GFP",
                "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP")) {
  
  targetTcr <- tcr %>% filter(genotype == cells) %>%
    arrange(desc(Freq)) %>%
    select(cdr3) %>%
    head(10) %>% tibble::deframe()

  ggplot(tcr[tcr$genotype == cells & tcr$cdr3 %in% targetTcr,]%>%
                 mutate(cdr3 = factor(cdr3, levels = targetTcr)),
               aes(x = Freq, y = reorder(cdr3, Freq), fill = cdr3)) +
    geom_bar(stat = "identity") +
    geom_label(aes(label = Freq), fill = "white") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          text = element_text(size = 15, face = 2))
  
  ggsave(paste0("3.Output/4.TCR/", gsub("/", "_", cells), "_TCR_bar.png"), 
         width = 8, height = length(targetTcr)/3, units = "in", 
         dpi = 300, bg = "white")
  # p2 <- ggplot() +
  #   geom_point(data= plotDf[plotDf$datNor.genotype == cells,][-6] %>% 
  #                mutate(datNor.newClst = as.character(datNor.newClst)), 
  #              aes(UMAP_1, UMAP_2, fill = datNor.newClst), 
  #              shape = 21, alpha = 0.2, size = 5) +
  #   geom_point(data = plotDf %>% filter(cdr3 %in% targetTcr) %>% 
  #                mutate(cdr3 = factor(cdr3, levels = targetTcr)),
  #              aes(UMAP_1, UMAP_2), col = "black", size = 2, alpha = 0.5)+
  #   # facet_wrap( ~ cdr3, ncol = 5) +
  #   coord_fixed(ratio = 1) +
  #   # scale_color_brewer(breaks = targetTcr) +
  #   guides(fill = "none") +
  #   theme_minimal() +
  #   theme(# legend.position = "none",
  #     text = element_blank(),
  #     title = element_blank(),
  #     panel.border = element_rect(colour = "black", fill=NA, size=0.2))
  
  
 ggplot() +
    geom_point(data= plotDf[plotDf$datNor.genotype == cells,][-6] %>% 
                 mutate(datNor.newClst = as.character(datNor.newClst)), 
               aes(UMAP_1, UMAP_2), fill = "gray10", 
               shape = 21, alpha = 0.3, size = 5) +
    geom_point(data = plotDf[plotDf$datNor.genotype == cells,] %>% 
                 filter(cdr3 %in% targetTcr) %>% 
                 mutate(cdr3 = factor(cdr3, levels = targetTcr)),
               aes(UMAP_1, UMAP_2, color = cdr3), size = 3, alpha = 0.8)+
    facet_wrap( ~ cdr3, ncol = 5) +
    coord_fixed(ratio = 1) +
    # scale_color_brewer(breaks = targetTcr) +
    guides(fill = "none") +
    theme_minimal() +
    theme(legend.position = "none",
          title = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = 15, face = 2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.2))
  
  ggsave(paste0("3.Output/4.TCR/", gsub("/", "_", cells), "_TCR_distr.png"), width = 12, 
         height = ifelse(length(targetTcr) >5, 5, 2.5), units = "in", 
         dpi = 300, bg = "white")
  
}

ggplot(t, aes(x = Freq, y = reorder(cdr3, Freq))) +
  geom_bar(stat = "identity") +
  facet_wrap(~genotype)

tSub <- t %>% filter(genotype == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP")

ggplot(tSub, aes(x = Freq, y = reorder(cdr3, Freq))) +
  geom_bar(stat = "identity")

targetTcr <- tSub$cdr3



# diversity evaluation ----

dfMeta <- datNor@meta.data

dfMeta$genotype <- gsub("/", "_", dfMeta$genotype)

contig_list <- list()

for (i in unique(dfMeta$genotype)[2:4]) {
  contig_list[[i]] <- dfMeta[dfMeta$genotype == i, 
                             c("barcode","is_cell", "contig_id", 
                               "high_confidence", "length", "chain", "v_gene", 
                               "d_gene", "j_gene", "c_gene", "full_length", 
                               "productive", "cdr3", "cdr3_nt", "reads", "umis", 
                               "raw_clonotype_id", "raw_consensus_id")]
}

combined <- combineTCR(contig_list, 
                       samples = names(contig_list))

# the number of unique clonotypes in each cell type_cluster
quantContig(combined, cloneCall="gene+nt", scale = T)

clonalDiversity(combined, cloneCall = "aa", n.boots = 100, scale = T)
