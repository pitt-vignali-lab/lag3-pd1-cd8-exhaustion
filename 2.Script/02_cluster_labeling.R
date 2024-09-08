# ----------------------------------------------------------------------------
# Setup | Assign cluster identities and rank marker genes
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/2.clstLab.R
# Input    : 1.Data/datInt.RDS
# Output   : 3.Output/2.ClstLab/makerGs.csv
#            3.Output/2.ClstLab/topRankedMarkerGs.csv
# Upstream : 01_data_integration.R
# Note     : Markers from presto wilcoxauc, adjusted P < 0.05
# ----------------------------------------------------------------------------

library(Seurat)
library(tidyverse)
library(presto)

if(!dir.exists("3.Output/2.ClstLab")){
  dir.create("3.Output/2.ClstLab")
}

dat <- readRDS("1.Data/datInt.RDS")

dat <- subset(dat, idents = 9, invert = T)

DefaultAssay(dat) <- "RNA"

freq1 = table(Idents(dat), dat$sample)
freq1[,2] = freq1[,2] * (table(dat$sample)[1]/table(dat$sample)[2])

prop.table(freq1, 1) %>% 
  data.frame() %>% 
  ggplot(aes(x = Freq, y = Var1, fill = Var2)) +
  geom_bar(stat = 'identity') +
  scale_y_discrete(limits=rev)

DimPlot(dat, label = T) +
  coord_fixed() +
  NoLegend()

# Marker genes ----
DimPlot(dat, 
        label = T) +
  coord_fixed()

deGs <- wilcoxauc(dat)

write.csv(deGs, "3.Output/2.ClstLab/makerGs.csv")

deGsSelec <- deGs %>% 
  filter(padj < 0.05,
         logFC > 0) %>% 
  group_by(group) %>%  
  top_n(n = 5, wt = logFC)

write.csv(deGsSelec, "3.Output/2.ClstLab/topRankedMarkerGs.csv")

gDf <- DotPlot(dat, 
        features = unique(deGsSelec$feature))$data %>% 
  rename('feature' = 'features.plot')

gDf <- merge(gDf, deGsSelec %>% select('feature', 'group'),
             by = 'feature')
 
gDf %>% 
  ggplot(aes(x = feature,
             y = id)) +
  geom_point(shape = 21,
             aes(size = pct.exp,
                 fill = avg.exp.scaled)) +
  facet_grid(~group, scales = 'free_x') +
  scale_fill_gradient2(high = 'red', 
                       mid = 'white', 
                       low = 'blue',
                       midpoint = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) 

ggsave("3.Output/2.ClstLab/top5MkGPerClst.png",
       width = 15, height = 5, units = 'in')


# Dot plot from Pauken et al ----

# naive and/or central memory-like 
DotPlot(dat, features = c(strsplit("Sell, Tcf7, Ccr7, Il7r, S1pr1, Klf2", ", ")[[1]], # high
                          strsplit("Klrg1, Cx3cr1, S1pr5, Tnf, Ifng, Gzmb, Prf1, Mki67, Pdcd1, Havcr2, Ctla4", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue")) +
  geom_rect(data = data.frame(start = c(0, 6.5), end = c(6.5, 17.5), 
                              Status = c("Shold be high", "should be low")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# effector and/or effector memory-like 3
DotPlot(dat, features = c(strsplit("Klrg1, Cx3cr1, S1pr5, Tnf, Ifng, Gzmb, Prf1", ", ")[[1]], # high
                          strsplit("Sell, Tcf7, Lef1, Ccr7, Il7r", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue"))+
  geom_rect(data = data.frame(start = c(0, 7.5), end = c(7.5, 12.5), 
                              Status = c("Shold be high", "Should be low")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Exhausted
DotPlot(dat, features = c(strsplit("Pdcd1, Havcr2, Lag3, Tigit", ", ")[[1]], # high
                          strsplit("Klrg1, Cx3cr1", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue")) +
  geom_rect(data = data.frame(start = c(0, 4.5), end = c(4.5, 6.5), 
                              Status = c("Shold be high", "Should be low")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Exhausted progenitor-like
DotPlot(dat, features = c(strsplit("Pdcd1, Lag3, Tigit, Tcf7", ", ")[[1]], # high
                          strsplit("Klrg1, Cx3cr1, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue")) +
  geom_rect(data = data.frame(start = c(0, 4.5), end = c(4.5, 7.5), 
                              Status = c("Shold be high", "Should be low")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Exhausted intermediate-like
DotPlot(dat, features = c(strsplit("Pdcd1, Ctla4, Lag3, Cd160", ", ")[[1]], # high
                          strsplit("Tcf7, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue")) +
  geom_rect(data = data.frame(start = c(0, 4.5), end = c(4.5, 6.5), 
                              Status = c("Shold be high", "Should be low")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Exhausted terminal-like
DotPlot(dat, features = c(strsplit("Pdcd1, Havcr2, Ctla4, Lag3, Cd160", ", ")[[1]]#, # high
                          # strsplit("Tcf7, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue")) +
  geom_rect(data = data.frame(start = c(0), end = c(5.5), 
                              Status = c("Shold be high")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# IFN-stimulated
DotPlot(dat, features = c(strsplit("Bst2, Irf1, Irf2, Irf7, Stat1, Stat2, Mx1", ", ")[[1]]#, # high
                          # strsplit("Tcf7, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue")) +
  geom_rect(data = data.frame(start = c(0), end = c(7.5), 
                              Status = c("Shold be high")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# cell cycle
DotPlot(dat, features = c(strsplit("Mki67, Cdk1, Cdk2", ", ")[[1]]#, # high
                          # strsplit("Tcf7, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue")) +
  geom_rect(data = data.frame(start = c(0), end = c(3.5), 
                              Status = c("Shold be high")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# resident memory-like
DotPlot(dat, features = c(strsplit("Itgae, Itga1, Cxcr6", ", ")[[1]]#, # high
                          # strsplit("Tcf7, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue")) +
  geom_rect(data = data.frame(start = c(0), end = c(3.5), 
                              Status = c("Shold be high")),
            aes(ymin=0, ymax=10.5, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# dotplot from Carmona et al ----
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7153840/

# Tcm markers (CD62Lhi, CD27hi, CD127hi, KLRG-1lo)

# exhaustion-related
DotPlot(dat, features = c("Pdcd1", "Havcr2", "Ctla4", "Tigit", "Lag3", "Batf", "Tox", "Gzmb", "Prf1", "Fasl"), col.min = 0, cols = c("white", "blue"))

# naïve/memory
DotPlot(dat, features = strsplit("Tcf7, Sell, Ccr7, Lef1, Il7r, Cd44", ", ")[[1]], col.min = 0, cols = c("white", "blue"))

DotPlot(dat, features = strsplit("Tcf7, Lef1, Il7r, Gzmk, Gzmb", ", ")[[1]], col.min = 0, cols = c("white", "blue"))

# exhausted
DotPlot(dat, features = c("Entpd1", "Havcr2", "Tigit", "Lag3", "Pdcd1", "Tnfsf13b", "Tox", "Id2"), col.min = 0, cols = c("white", "blue"))

DotPlot(dat, features = c("Ccl3", "Prf1", "Fasl", "Gzmb", "Ccl5"), col.min = 0, cols = c("white", "blue"))

DotPlot(dat, features = c("Entpd1", "Havcr2", "Tigit", "Lag3", "Pdcd1", "Tnfsf13b", "Tox", "Id2"), col.min = 0, cols = c("white", "blue"))