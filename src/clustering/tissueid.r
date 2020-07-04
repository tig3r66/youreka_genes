library(tidyr)
library(dplyr)
library(stringr)
library(mclust)
library(viridis)
library(factoextra)  # PCA
library(cowplot)     # for consistently-sized plots
library(plotly)      # for 3d pca plot

# READING CELL LINE DATA AND APPENDING COSMIC IDENTIFIER ================
## GDSC data
dfGdsc <- read.csv("temp1GDSC.csv")
dfGdsc <- dfGdsc[, which(colMeans(!is.na(dfGdsc)) > 0.8)]
dfGdsc <- na.omit(dfGdsc)
dfGdsc$Cell_name <- gsub("\\.", "-", dfGdsc$Cell_name)
## cell lines with cosmic identifiers
meta_data <- read.csv("cells_new.csv")
meta_data <- (meta_data %>% rename(Cell_name=cellid))
## merging on "Cell_name"
dfGdsc <- merge(dfGdsc, meta_data, "Cell_name")
## normalizing data
scaledGdsc <- as.data.frame(t(scale(t(scale(dfGdsc[, 2:117])))))
dfGdsc[2:117] <- scaledGdsc


# PCA Plot ============
## PCA
pc <- princomp(dfGdsc[2:117])
comp <- data.frame(pc$scores[, 1:2])
## K-means
fit <- kmeans(dfGdsc[2:117], centers=2, iter.max=10, nstart=25)
class <- as.data.frame(fit$cluster)
class$CELL_LINE_NAME <- rownames(class)
dfGdsc$CELL_LINE_NAME <- rownames(dfGdsc)
dfGdsc <- merge(dfGdsc, class, "CELL_LINE_NAME")
rownames(dfGdsc) <- dfGdsc$CELL_LINE_NAME
dfGdsc$pred <- dfGdsc$`fit$cluster`

# # PCA PLOTS ============
# # tumour type plot
# tumourTypePlot <- fviz_pca_ind(pc, geom="point", habillage=dfGdsc$solidity,
#                              alpha.ind=1, pointsize=2, legend.title="Tumour\nType",
#                              palette=c("#FF0000", "#2565AE")) +
#     theme_bw() +
#     theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
#     theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
#     theme(legend.text=element_text(size=12), legend.title=element_text(size=15)) +
#     ggtitle("")
# tumourTypePlot
# 
# # tissue id plot
# tissueIdPlot <- fviz_pca_ind(pc, geom="point", habillage=str_wrap(dfGdsc$body.system, 20),
#                              alpha.ind=1, pointsize=1.5, palette=magma(12)) +
#     theme_bw() +
#     theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
#     theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
#     theme(legend.text=element_text(size=12), legend.title=element_text(size=15)) +
#     ggtitle("")
# tissueIdPlot

over_comps <- data.frame("PC1"=pc$scores[,1], "PC2"=pc$scores[,2], "PC3"=pc$scores[,3])

# body system plot
bsplot <- plot_ly(over_comps, x=~PC1, y=~PC2, z=~PC3,
                  type="scatter3d", mode="markers", marker=list(size=3),
                  color=str_wrap(dfGdsc$body.system,15), colors=magma(12))
bsplot <- bsplot %>%
    layout(
        legend=list(font=list(size=20), itemsizing="constant", orientation="h", x=0, y=0),
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
bsplot


# tumour type plot
tumourplot <- plot_ly(over_comps, x=~PC1, y=~PC2, z=~PC3,
                   type="scatter3d", mode="markers", marker=list(size=3),
                   color=dfGdsc$solidity, colors=c("#2565AE", "#FF0000", "#2565AE", "#FF0000"))
tumourplot <- tumourplot %>%
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
tumourplot
