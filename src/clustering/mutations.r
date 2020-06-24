library(tidyr)
library(dplyr)
library(mclust)
library(factoextra)  # PCA
library(cowplot)     # for consistently-sized plots

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
dfGdsc$KRAS <- as.factor(ifelse(dfGdsc$KRAS=="0", "no mutation", "mutation"))
dfGdsc$TP53 <- as.factor(ifelse(dfGdsc$TP53=="0", "no mutation", "mutation"))


# PCA Plot ============
## PCA
pc <- princomp(dfGdsc[2:117])
## KRAS
KRASPlot <- fviz_pca_ind(pc, geom="point", habillage=dfGdsc$KRAS,
                         alpha.ind=1, pointsize=2, legend.title="KRAS Mutation\nStatus",
                         palette=c("#FF0000", "#2565AE")) +
    theme_bw() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(legend.text=element_text(size=12), legend.title=element_text(size=15)) +
    ggtitle("")
KRASPlot

## TP53 plot
TP53Plot <- fviz_pca_ind(pc, geom="point", habillage=dfGdsc$TP53,
                         alpha.ind=1, pointsize=2, legend.title="TP53 Mutation\nStatus",
                         palette=c("#FF0000", "#2565AE")) +
    theme_bw() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    theme(legend.text=element_text(size=12), legend.title=element_text(size=15)) +
    ggtitle("")
TP53Plot

# MERGING PLOTS ACROSS tissueid.r, mutations.r (need to run all files)
plot_grid(tissueIdPlot,
          tumourTypePlot + theme(legend.justification="left"),
          TP53Plot + theme(legend.justification="left"),
          KRASPlot + theme(legend.justification="left"),
          nrow=2, ncol=2, labels=c("A", "B", "C", "D"), align="v")

