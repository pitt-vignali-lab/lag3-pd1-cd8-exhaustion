# ----------------------------------------------------------------------------
# Fig 2 | Render the published sub-clustering panels
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/PictureFroManual/subclustering.R
# Input    : 1.Data/datIntNoSmallClst.RDS
# Output   : 3.Output/FormalFigures/subClusteringSplit.png
#            3.Output/CellFigure/UMAPSubcluster*.png
# Upstream : 03_subclustering.R
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(extrafont)
library(scales)
library(grid)
library(ggrepel)

dat <- readRDS("1.Data/datInt.RDS")

subClst <- read.csv("3.Output/3.SubClst/clstrsforsubClst_0_7.csv")

datSub <- subset(dat, cells = subClst$X)

datSub <- AddMetaData(datSub, 
                      subClst %>% 
                        column_to_rownames(var = "X") %>% 
                        rename(newClst = 1))

DimPlot(datSub, group.by = "newClst") +
  coord_fixed()

df <- data.frame(x = datSub[["umap"]]@cell.embeddings[,1],
                 y = datSub[["umap"]]@cell.embeddings[,2],
                 sample = datSub$sample,
                 tissue = datSub$tissue,
                 genotype = datSub$genotype,
                 newLab = datSub$newClst) |> 
  filter(y > -3)

df$newLab <- factor(paste0("SubC", df$newLab) , levels = paste0("SubC", 0:3))

df$genotype <- factor(df$genotype, levels = c("E8ICre.GFP (WT)",
                                              "Pdcd1L/L Lag3L/L-yfp (WT)",
                                              "Lag3L/L-yfp E8ICre.GFP",
                                              "Pdcd1L/L E8ICre.GFP",
                                              "Pdcd1L/L Lag3L/L-yfp E8ICre.GFP"))

dfLab <- df %>% group_by(newLab) %>% 
  summarise(x = mean(x), y = mean(y))

ggplot() +
  geom_point(data = df, aes(x , y, col = newLab), alpha = 0.5) +
  coord_fixed(ratio = 1) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_bw() +
  guides(color = guide_legend(title = "Clusters", 
                              nrow=2,
                              byrow=TRUE,
                              override.aes = list(size=3))) +
  scale_color_brewer(palette="Set2") +
  theme(legend.position = "none", # c(0.6, 0.1),
        legend.direction = "vertical",
        legend.title.align = 0.5,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key=element_rect(fill = "gray90", color = NA),
        legend.text = element_text(size = 11.5, 
                                   face = 2),
        legend.title = element_text(size = 12.5, 
                                    face = 2))

ggsave(paste0("3.Output/FormalFigures/", Sys.Date(), "_subClusteringAll.png"), 
       width = 4, height = 4, units = "in", 
       dpi = 300, bg = "white")



p <- ggplot()

for (i in unique(df$newLab)) {
  p = p +
    geom_point(data = df %>% filter(newLab == i), 
               aes(x, y, fill = newLab), shape = 21, stroke = 1) +
    geom_point(data = df %>% filter(newLab == i), 
               aes(x, y, col = newLab), shape = 16)
}

# wwith label 
p +
  # geom_label_repel(data = labDf,
  #                  aes(UMAP_1, UMAP_2, label = clst, col = clst),
  #                  size = 5, fontface = 2) +
  scale_color_brewer(palette="Set2") +
  geom_label(data = dfLab %>% mutate(labs = LETTERS[1:4]), 
                   aes(x, y, label = labs, col = newLab), 
             fontface = 2, size = 8) +
  coord_fixed() +
  theme_bw() +
  ggforce::theme_no_axes() +
  NoLegend()

ggsave('3.Output/CellFigure/UMAPSubclusterOverall.png', width = 5, height = 5, units = 'in')

splitp <- p +
  scale_color_brewer(palette="Set2") +
  # geom_label(data = dfLab %>% mutate(labs = LETTERS[1:4]), 
  #            aes(x, y, label = labs, col = newLab), 
  #            fontface = 2, size = 8) +
  facet_grid(sample ~ genotype, switch = "y") +
  coord_fixed() +
  theme_bw() +
  ggforce::theme_no_axes() +
  NoLegend() +
  theme(text = element_blank())

g <- ggplotGrob(splitp)
rm_grobs <- g$layout$name %in% c("panel-2-1")
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]
# g$layout[g$layout$name == "axis-b-3-3", c("t", "b")] = c(14.5, 14.5)
grid.newpage()

png(paste0("3.Output/CellFigure/UMAPSubclusterSplit.png"), 
    width = 8, height = 4, units = "in", res = 300)
grid.draw(g)
dev.off()




p <- ggplot() +
  geom_point(data = df, aes(x , y, col = newLab), stroke = 0) +
  facet_grid(sample ~ genotype, switch = "y")+
  scale_color_brewer(palette="Set2") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, 
                                  face = 2))


p

g <- ggplotGrob(p)
rm_grobs <- g$layout$name %in% c("panel-2-1")
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]
# g$layout[g$layout$name == "axis-b-3-3", c("t", "b")] = c(14.5, 14.5)
grid.newpage()

png(paste0("3.Output/FormalFigures/", Sys.Date(), "_subClusteringSplit.png"), 
    width = 8, height = 4, units = "in", res = 300)
grid.draw(g)
dev.off()


ggsave("3.Output/FormalFigures/subClusteringSplit.png", 
       width = 8, height = 5, units = "in", dpi = 300, bg = "white")
