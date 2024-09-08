# ----------------------------------------------------------------------------
# Fig 4 | Clonotype abundance and diversity by genotype
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/4.4TCRComboAbundance.R
# Input    : 1.Data/tcrCombinedUpdated.csv
# Output   : 3.Output/4.TCR/TCR_clonaltypes_abundance*.png
# Upstream : 29_fig4_tcr_combined.R
# ----------------------------------------------------------------------------

if(!dir.exists("3.Output/4.TCR")){
  dir.create("3.Output/4.TCR")
}

library(Seurat)
library(tidyverse)
library(scales)
library(MASS)

dat <- readRDS("1.Data/datInt.RDS")

newClst <- read.csv("1.Data/TCR_defined_clusters.csv") %>%
  column_to_rownames(var = "X") %>%
  rename(newClst = dat.newClst)

dat <- AddMetaData(dat, newClst)

Idents(dat) <- "newClst"

datNor <- subset(dat, sample == "Normal CD8")

tcrDf <- read.csv("1.Data/tcrCombinedUpdated.csv")[-1] %>%
  # mutate(barcode = paste0("nor_", barcode)) %>% 
  filter(barcode %in% rownames(datNor@meta.data)) %>% 
  rename(cmbdCdr3 = cdr3)

unique(datNor$genotype)

abdDf <- data.frame()
for (i in unique(datNor$genotype)) {
  tcrDfSub <- tcrDf %>% 
    filter(barcode %in% Cells(subset(datNor, genotype == i)))
  
  tabDf <- data.frame(table(tcrDfSub$cmbdCdr3)) %>% 
    arrange(Freq = desc(Freq))
  tabDf$index <- 1:nrow(tabDf)
  tabDf$sample <- i
  
  abdDf <- rbind(abdDf, tabDf)
}

abdDf$index <- abs(abdDf$index - max(abdDf$index) + 1)

abdDf$sample <- factor(abdDf$sample, 
                       levels = c("E8ICre.GFP (WT)",
                                  "Pdcd1L/L Lag3L/L-yfp (WT)",
                                  "Lag3L/L-yfp E8ICre.GFP",
                                  "Pdcd1L/L E8ICre.GFP",
                                  "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))



abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp (WT)", ]$index <- abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp (WT)", ]$index - min(abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp (WT)", ]$index)
abdDf[abdDf$sample == "E8ICre.GFP (WT)", ]$index <- abdDf[abdDf$sample == "E8ICre.GFP (WT)", ]$index - min(abdDf[abdDf$sample == "E8ICre.GFP (WT)", ]$index)
abdDf[abdDf$sample == "Lag3L/L-yfp E8ICre.GFP", ]$index <- abdDf[abdDf$sample == "Lag3L/L-yfp E8ICre.GFP", ]$index - min(abdDf[abdDf$sample == "Lag3L/L-yfp E8ICre.GFP", ]$index)
abdDf[abdDf$sample == "Pdcd1L/L E8ICre.GFP", ]$index <- abdDf[abdDf$sample == "Pdcd1L/L E8ICre.GFP", ]$index - min(abdDf[abdDf$sample == "Pdcd1L/L E8ICre.GFP", ]$index) 
abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", ]$index <- abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", ]$index - min(abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", ]$index)

abdDf %>% group_by(sample) %>% top_n(1, Freq)

ggplot(abdDf, aes(Freq, index, col = sample,
                  alpha = ifelse(Freq < 3, 0.0001, 0.5))) +
  geom_point(size = 2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values=c( "black", "gray40","blue", "red", "purple")) +
  geom_vline(xintercept = abdDf %>% 
               group_by(sample) %>% 
               top_n(1, Freq) %>% 
               filter(sample %in% c('E8ICre.GFP (WT)', 
                                    'Pdcd1L/L Lag3L/L-yfp (WT)')) %>% 
               pull(Freq),
             col = "#c21807",
             lty = 2, linewidth = 0.8) +
  geom_hline(yintercept = abdDf %>% 
               group_by(sample) %>% 
               top_n(1, index) %>% 
               filter(sample %in% c('E8ICre.GFP (WT)', 
                                    'Pdcd1L/L Lag3L/L-yfp (WT)')) %>% 
               pull(index),
             col = "#c21807",
             lty = 2, linewidth = 0.8) +
  facet_grid(~sample, ) +
  coord_fixed(1/200) +
  labs(x = "Number of clonotype", y = "Abundance") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.spacing = unit(0.5, "cm"),
        text = element_text(size = 15, face = 2),
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "cm")))

ggsave("3.Output/4.TCR/TCR_clonaltypes_abundance_with_stdLine(Updated).png", 
       height = 4, width = 15, units = "in", dpi = 300, bg = "white")


# abdDf$newN <- NA
# abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp (WT)", ]$newN <- abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp (WT)", ]$Freq + 0
# abdDf[abdDf$sample == "E8ICre.GFP (WT)", ]$newN <- abdDf[abdDf$sample == "E8ICre.GFP (WT)", ]$Freq + 10^0.1 - 1
# abdDf[abdDf$sample == "Lag3L/L-yfp E8ICre.GFP", ]$newN <- abdDf[abdDf$sample == "Lag3L/L-yfp E8ICre.GFP", ]$Freq + 10^0.2 - 1
# abdDf[abdDf$sample == "Pdcd1L/L E8ICre.GFP", ]$newN <- abdDf[abdDf$sample == "Pdcd1L/L E8ICre.GFP", ]$Freq + 10^0.3 -1 
# abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", ]$newN <- abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", ]$Freq + 10^0.4 -1 
# 
# 
# abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp (WT)", ]$index <- abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp (WT)", ]$index - min(abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp (WT)", ]$index)
# abdDf[abdDf$sample == "E8ICre.GFP (WT)", ]$index <- abdDf[abdDf$sample == "E8ICre.GFP (WT)", ]$index - min(abdDf[abdDf$sample == "E8ICre.GFP (WT)", ]$index)
# abdDf[abdDf$sample == "Lag3L/L-yfp E8ICre.GFP", ]$index <- abdDf[abdDf$sample == "Lag3L/L-yfp E8ICre.GFP", ]$index - min(abdDf[abdDf$sample == "Lag3L/L-yfp E8ICre.GFP", ]$index)
# abdDf[abdDf$sample == "Pdcd1L/L E8ICre.GFP", ]$index <- abdDf[abdDf$sample == "Pdcd1L/L E8ICre.GFP", ]$index - min(abdDf[abdDf$sample == "Pdcd1L/L E8ICre.GFP", ]$index) 
# abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", ]$index <- abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", ]$index - min(abdDf[abdDf$sample == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", ]$index)
# 
# abdDf |> group_by(sample) |> summarise(max(Freq))
# 
# abdDf |> group_by(sample) |> summarise(max(index))
# 
# ggplot(abdDf, aes(newN, index, col = sample)) +
#   geom_point(size = 1, alpha = 0.8) +
#   scale_x_log10(breaks = c(1, 10^0.1, 10^0.2, 10^0.3, 10^0.4, 283, 296+10^0.1-1, 87+10^0.2-1, 65+10^0.3-1, 34+10^0.4-1),
#                 labels = c(1, 1, 1, 1, 1, 283, 296, 87, 65, 341)) +
#   geom_vline(xintercept = c(283, 296+10^0.1-1, 87+10^0.2-1, 65+10^0.3-1, 34+10^0.4-1),
#              col = c("gray60", "black", "blue", "red", "purple"),
#              lty = 2) +
#   geom_hline(yintercept = c(531, 610, 889, 846, 1195),
#              col = c("gray60", "black", "blue", "red", "purple"),
#              lty = 2) +
#   # scale_y_reverse() +
#   # facet_wrap(~sample, nrow = 1) +
#   geom_smooth() +
#   labs(x = "Number of clonotype", y = "Abundance") +
#   scale_color_manual(values=c("gray60", "black", "blue", "red", "purple")) +
#   theme_bw() +
#   theme(legend.position = "none",
#         strip.background = element_blank(),
#         text = element_text(size = 15, face = 2),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, 
#                                    colour = rep(c("gray60", "black", "blue", "red", "purple"), 2)))
# ggsave("3.Output/4.TCR/TCR_clonaltypes_abundance.png", 
#        height = 4, width = 15, units = "in", dpi = 300, bg = "white")
# 
# # ----
# df <- data.frame(dat.anno = datNor$genotype,
#            cdr3 = datNor$cdr3) %>% 
#   filter(cdr3 != "None",
#          cdr3 != ".",
#          !is.na(cdr3))  %>% 
#   group_by(cdr3, dat.anno) %>% 
#   summarise(n = n()) %>% 
#   mutate(dat.anno =  factor(dat.anno, 
#                             levels = c("Pdcd1L/L Lag3L/L-yfp (WT)",
#                                        "E8ICre.GFP (WT)",
#                                        "Lag3L/L-yfp E8ICre.GFP",
#                                        "Pdcd1L/L E8ICre.GFP",
#                                        "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))) |> 
#   mutate(n = as.numeric(n))
# 
# df$newN <- NA
# df[df$dat.anno == "Pdcd1L/L Lag3L/L-yfp (WT)", ]$newN <- df[df$dat.anno == "Pdcd1L/L Lag3L/L-yfp (WT)", ]$n + 0
# df[df$dat.anno == "E8ICre.GFP (WT)", ]$newN <- df[df$dat.anno == "E8ICre.GFP (WT)", ]$n + 5
# df[df$dat.anno == "Lag3L/L-yfp E8ICre.GFP", ]$newN <- df[df$dat.anno == "Lag3L/L-yfp E8ICre.GFP", ]$n + 10
# df[df$dat.anno == "Pdcd1L/L E8ICre.GFP", ]$newN <- df[df$dat.anno == "Pdcd1L/L E8ICre.GFP", ]$n + 15
# df[df$dat.anno == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", ]$newN <- df[df$dat.anno == "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP", ]$n + 20
#   
# ggplot(df, aes(x  = n, y = reorder(n, n), group = dat.anno, col = dat.anno)) +
#   geom_line(size = 1)+
#   # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#   #               labels = trans_format("log10", math_format(10^.x))) +
#   facet_wrap(~dat.anno, nrow = 1) +
#   labs(x = "Genotypes", y = "Abundance") +
#   scale_color_manual(values=c("gray60", "black", "blue", "red", "purple")) +
#   theme_bw() +
#   #coord_fixed(1/10) +
#   theme(# legend.position = "bottom",
#         strip.background = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         text = element_text(size = 15, face = 2))
# 
# ggsave(paste0("3.Output/4.TCR/", Sys.Date(), "TCR_clonaltypes_abundance.png"), 
#        height = 6, width = 9, units = "in", dpi = 300, bg = "white")
# # ----
# data.frame(dat.anno = datNor$genotype,
#            cdr3 = datNor$cdr3) %>%
#   filter(cdr3 != "None",
#          cdr3 != ".",
#          !is.na(cdr3))  %>% 
#   group_by(cdr3, dat.anno) %>% 
#   summarise(n = n()) %>% 
#   # filter(dat.anno == "Pdcd1L/L Lag3L/L-yfp (WT)") |> 
#   ggplot(aes(x  = n, y = reorder(n, n), group = dat.anno, col = dat.anno)) +
#   geom_line()
# 
