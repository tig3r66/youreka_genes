library(tidyr)
library(dplyr)
library(mclust)
library(factoextra)  # PCA
library(cowplot)     # for consistently-sized plots
library(plotly)      # for 3d PCA plots

# READING CELL LINE DATA AND APPENDING COSMIC IDENTIFIER ================
## GDSC data
dfGdsc <- read.csv("temp1GDSC.csv")
dfGdsc <- dfGdsc[, which(colMeans(!is.na(dfGdsc)) > 0.8)]
dfGdsc <- na.omit(dfGdsc)
dfGdsc$Cell_name <- gsub("\\.", "-", dfGdsc$Cell_name)
## cell lines with cosmic identifiers
meta_data <- read.csv("mut.csv")
meta_data <- (meta_data %>% rename(Cell_name=cellid))
meta_data <- meta_data[c("Cell_name", "TP53", "KRAS")]
## merging on "Cell_name"
dfGdsc <- merge(dfGdsc, meta_data, "Cell_name")
## normalizing data
scaledGdsc <- as.data.frame(t(scale(t(scale(dfGdsc[, 2:117])))))
dfGdsc[2:117] <- scaledGdsc
dfGdsc <- as.data.frame(na.omit(dfGdsc))
dfGdsc$KRAS <- as.factor(ifelse(dfGdsc$KRAS=="0", "No mutation", "Mutation"))
dfGdsc$TP53 <- as.factor(ifelse(dfGdsc$TP53=="0", "No mutation", "Mutation"))


# PCA Plot ============
## PCA
pc <- princomp(dfGdsc[2:117])

# plotly PCA
main_comps <- data.frame("PC1"=pc$scores[,1], "PC2"=pc$scores[,2], "PC3"=pc$scores[,3])

# TP53 PCA
tp53fig <- plot_ly(main_comps, x=~PC1, y=~PC2, z=~PC3,
                   type="scatter3d", mode="markers", marker=list(size=3),
                   color=dfGdsc$TP53, colors=c("#FF0000", "#2565AE"))
tp53fig <- tp53fig %>%
    layout(
        legend=list(font=list(size=20), itemsizing="constant", orientation="h", x=0.3, y=0),
        font=list(size=14),
        scene=list(
            xaxis=list(title="PC1 (12.8%)",
                       zerolinecolor="white",
                       backgroundcolor="rgb(200, 200, 230)",
                       showbackground=T,
                       gridcolor="white"),
            yaxis=list(title="PC2 (12.6%)",
                       zerolinecolor="white",
                       backgroundcolor="rgb(230, 200, 230)",
                       showbackground=T,
                       gridcolor="white"),
            zaxis=list(title="PC3 (10.5%)",
                       zerolinecolor="white",
                       backgroundcolor="rgb(230, 230, 200)",
                       showbackground=T,
                       gridcolor="white"),
            camera=list(eye=list(x=-0.75,y=2,z=1))
        ))
tp53fig

# KRAS PCA
krasfig <- plot_ly(main_comps, x=~PC1, y=~PC2, z=~PC3,
                   type="scatter3d", mode="markers", marker=list(size=3),
                   color=dfGdsc$KRAS, colors=c("#FF0000", "#2565AE"))
krasfig <- krasfig %>%
    layout(
        legend=list(font=list(size=20), itemsizing="constant", orientation="h", x=0.3, y=0),
        font=list(size=14),
        scene=list(
            xaxis=list(title="PC1 (12.8%)",
                       zerolinecolor="white",
                       backgroundcolor="rgb(200, 200, 230)",
                       showbackground=T,
                       gridcolor="white"),
            yaxis=list(title="PC2 (12.6%)",
                       zerolinecolor="white",
                       backgroundcolor="rgb(230, 200, 230)",
                       showbackground=T,
                       gridcolor="white"),
            zaxis=list(title="PC3 (10.5%)",
                       zerolinecolor="white",
                       backgroundcolor="rgb(230, 230, 200)",
                       showbackground=T,
                       gridcolor="white"),
            camera=list(eye=list(x=-0.75,y=2,z=1))
        ))
krasfig


# fviz PCA (2d) ============
# ## KRAS
# KRASPlot <- fviz_pca_ind(pc, geom="point", habillage=dfGdsc$KRAS,
#                          alpha.ind=1, pointsize=2, legend.title="KRAS\nMutation\nStatus",
#                          palette=c("#FF0000", "#2565AE")) +
#     theme_bw() +
#     theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
#     theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
#     theme(legend.text=element_text(size=12), legend.title=element_text(size=15)) +
#     ggtitle("")
# KRASPlot
# 
# ## TP53 plot
# TP53Plot <- fviz_pca_ind(pc, geom="point", habillage=dfGdsc$TP53,
#                          alpha.ind=1, pointsize=2, legend.title="TP53\nMutation\nStatus",
#                          palette=c("#FF0000", "#2565AE")) +
#     theme_bw() +
#     theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
#     theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
#     theme(legend.text=element_text(size=12), legend.title=element_text(size=15)) +
#     ggtitle("")
# TP53Plot
# 
# # MERGING PLOTS ACROSS tissueid.r, mutations.r (need to run all files)
# plot_grid(tissueIdPlot,
#           tumourTypePlot + theme(legend.justification="left"),
#           TP53Plot + theme(legend.justification="left"),
#           KRASPlot + theme(legend.justification="left"),
#           nrow=2, ncol=2, labels=c("A", "B", "C", "D"), align="v")
# 
# # 3D PCA plots
# library(pca3d)
# # tp53
# pca3d(pc, group=dfGdsc$TP53, show.ellipses=F, ellipse.ci=0.95,
#       legend="bottomleft", show.plane=F,
#       radius=0.75, axes.color="black", palette=c("#FF0000", "#2565AE"))
# snapshotPCA3d(file="tp53.png")
