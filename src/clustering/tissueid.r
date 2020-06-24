library(tidyr)
library(dplyr)
library(stringr)
library(mclust)
library(viridis)
library(factoextra)  # PCA
library(cowplot)     # for consistently-sized plots

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

# PCA PLOTS ============
# tumour type plot
tumourTypePlot <- fviz_pca_ind(pc, geom="point", habillage=dfGdsc$solidity,
                             alpha.ind=1, pointsize=2, legend.title="Tumour\nType",
                             palette=c("#FF0000", "#2565AE")) +
    theme_bw() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(legend.text=element_text(size=12), legend.title=element_text(size=15)) +
    ggtitle("")
tumourTypePlot

# tissue id plot
tissueIdPlot <- fviz_pca_ind(pc, geom="point", habillage=str_wrap(dfGdsc$body.system, 20),
                             alpha.ind=1, pointsize=1.5, palette=magma(12)) +
    theme_bw() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(legend.text=element_text(size=12), legend.title=element_text(size=15)) +
    ggtitle("")
tissueIdPlot
