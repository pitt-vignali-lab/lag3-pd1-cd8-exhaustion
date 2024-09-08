# ----------------------------------------------------------------------------
# Setup | Marker gene expression across clusters
#
# Paper    : LAG-3 and PD-1 synergize on CD8+ T cells to drive T cell
#            exhaustion and hinder autocrine IFN-gamma-dependent anti-tumor
#            immunity
# Original : 2.Script/markers.R
# Input    : 1.Data/datInt.RDS
# Output   : 3.Output/2.ClstLab/
# Upstream : 02_cluster_labeling.R
# ----------------------------------------------------------------------------

# naive and/or central memory–like 
DotPlot(datInt, features = c(strsplit("Sell, Tcf7, Ccr7, Il7r, S1pr1, Klf2", ", ")[[1]], # high
                          strsplit("Klrg1, Cx3cr1, S1pr5, Tnf, Ifng, Gzmb, Prf1, Mki67, Pdcd1, Havcr2, Ctla4", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue"), assay = "RNA") +
  geom_rect(data = data.frame(start = c(0, 6.5), end = c(6.5, 17.5), 
                              Status = c("Shold be high", "should be low")),
            aes(ymin=0, ymax=12, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# effector and/or effector memory–like 3
DotPlot(datInt, features = c(strsplit("Klrg1, Cx3cr1, S1pr5, Tnf, Ifng, Gzmb, Prf1", ", ")[[1]], # high
                          strsplit("Sell, Tcf7, Lef1, Ccr7, Il7r", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue"), assay = "RNA")+
  geom_rect(data = data.frame(start = c(0, 7.5), end = c(7.5, 12.5), 
                              Status = c("Shold be high", "Should be low")),
            aes(ymin=0, ymax=12, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Exhausted
DotPlot(datInt, features = c(strsplit("Pdcd1, Havcr2, Lag3, Tigit", ", ")[[1]], # high
                          strsplit("Klrg1, Cx3cr1", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue"), assay = "RNA") +
  geom_rect(data = data.frame(start = c(0, 4.5), end = c(4.5, 6.5), 
                              Status = c("Shold be high", "Should be low")),
            aes(ymin=0, ymax=12, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Exhausted progenitor-like
DotPlot(datInt, features = c(strsplit("Pdcd1, Lag3, Tigit, Tcf7", ", ")[[1]], # high
                          strsplit("Klrg1, Cx3cr1, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue"), assay = "RNA") +
  geom_rect(data = data.frame(start = c(0, 4.5), end = c(4.5, 7.5), 
                              Status = c("Shold be high", "Should be low")),
            aes(ymin=0, ymax=12, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Exhausted intermediate-like
DotPlot(datInt, features = c(strsplit("Pdcd1, Ctla4, Lag3, Cd160", ", ")[[1]], # high
                          strsplit("Tcf7, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue"), assay = "RNA") +
  geom_rect(data = data.frame(start = c(0, 4.5), end = c(4.5, 6.5), 
                              Status = c("Shold be high", "Should be low")),
            aes(ymin=0, ymax=12, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Exhausted terminal-like
DotPlot(datInt, features = c(strsplit("Pdcd1, Havcr2, Ctla4, Lag3, Cd160", ", ")[[1]]#, # high
                          # strsplit("Tcf7, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue"), assay = "RNA") +
  geom_rect(data = data.frame(start = c(0), end = c(5.5), 
                              Status = c("Shold be high")),
            aes(ymin=0, ymax=12, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# IFN-stimulated
DotPlot(datInt, features = c(strsplit("Bst2, Irf1, Irf2, Irf7, Stat1, Stat2, Mx1", ", ")[[1]]#, # high
                          # strsplit("Tcf7, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue"), assay = "RNA") +
  geom_rect(data = data.frame(start = c(0), end = c(7.5), 
                              Status = c("Shold be high")),
            aes(ymin=0, ymax=12, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# cell cycle
DotPlot(datInt, features = c(strsplit("Mki67, Cdk1, Cdk2", ", ")[[1]]#, # high
                          # strsplit("Tcf7, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue"), assay = "RNA") +
  geom_rect(data = data.frame(start = c(0), end = c(3.5), 
                              Status = c("Shold be high")),
            aes(ymin=0, ymax=12, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# resident memory–like
DotPlot(datInt, features = c(strsplit("Itgae, Itga1, Cxcr6", ", ")[[1]]#, # high
                          # strsplit("Tcf7, Havcr2", ", ")[[1]] # low
), col.min = 0, cols = c("white", "blue"), assay = "RNA") +
  geom_rect(data = data.frame(start = c(0), end = c(3.5), 
                              Status = c("Shold be high")),
            aes(ymin=0, ymax=12, xmin=start,
                xmax=end, fill = Status), alpha =0.2, inherit.aes = FALSE ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
