library(tidyr)
library(dplyr)
library(mclust)
library(ggplotify)
library(factoextra)  # PCA
library(cowplot)     # for consistently-sized plots

# READING CELL LINE DATA AND APPENDING COSMIC IDENTIFIER ================
## GDSC data
dfGdsc <- read.csv("temp1GDSC.csv")
dfGdsc <- dfGdsc[, which(colMeans(!is.na(dfGdsc)) > 0.8)]
dfGdsc <- na.omit(dfGdsc)
dfGdsc$Cell_name <- gsub("\\.", "-", dfGdsc$Cell_name)
## cell lines with cosmic identifiers
meta_data <- read.csv("cells.csv")
meta_data <- (meta_data %>% rename(Cell_name=Sample.Name))
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

# re-running PCA
dfGdsc$pred <- as.factor(ifelse(class$`fit$cluster`==1, "A", "B"))
dfGdsc$type <- as.factor(ifelse(dfGdsc$Growth.Properties=="Suspension", "solid", "non-solid"))

# PCA PLOTS ============
# therapeutic response clusters
responseGroupPlot <- fviz_pca_ind(pc, geom="point", habillage=dfGdsc$pred,
                                  addEllipses=T, ellipse.level=0.95,
                                  alpha.ind=0.75, pointsize=1,
                                  legend.title="Response\ncluster",
                                  palette=c("#FF0000", "#2565AE")) +
    theme_bw() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(legend.text=element_text(size=12), legend.title=element_text(size=15)) +
    ggtitle("")
responseGroupPlot


# GRID (NEED TO RUN GDSC2.R AND PCA.R TO GET PLOT OBJECTS) ============
cowResp <- plot_grid(responseGroupPlot, basic,
                     ncol=2, nrow=1, align="h",
                     labels=c("A", "B"))
cowHeatmap <- plot_grid(as.grob(heatmapPlot), nrow=1, labels=c("C"))
cowCombined <- plot_grid(cowResp, cowHeatmap, nrow=2, ncol=1)
cowCombined
