library(tidyr)
library(dplyr)
library(mclust)
library(ggplotify)
library(factoextra)  # PCA
library(cowplot)     # for consistently-sized plots
library(plotly)      # for 3D PCA plots

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
zoom <- 1.15
main_comps <- data.frame("PC1"=pc$scores[,1], "PC2"=pc$scores[,2], "PC3"=pc$scores[,3])
respPlot <- plot_ly(main_comps, x=~PC1, y=~PC2, z=~PC3,
                   type="scatter3d", mode="markers", marker=list(size=3),
                   color=dfGdsc$pred, colors=c("#FF0000", "#2565AE"))
respPlot <- respPlot %>%
    layout(
        legend=list(font=list(size=20), itemsizing="constant", orientation="h", x=0.4, y=0),
        font=list(size=14),
        scene=list(
            xaxis=list(title="PC1 (12.8%)",
                       zerolinecolor="white",
                       backgroundcolor="rgb(200, 200, 230)",
                       showbackground=T,
                       gridcolor="white"),
            yaxis=list(title="PC2 (13.4%)",
                       zerolinecolor="white",
                       backgroundcolor="rgb(230, 200, 230)",
                       showbackground=T,
                       gridcolor="white"),
            zaxis=list(title="PC3 (10.6%)",
                       zerolinecolor="white",
                       backgroundcolor="rgb(230, 230, 200)",
                       showbackground=T,
                       gridcolor="white"),
            camera=list(eye=list(x=2*zoom,y=0.1*zoom,z=0.75))
        ))
respPlot

# for gif
# cam_coords <- list(x=2, y=0.1, z=0.75)
# gifnum <- 0
# for (i in seq(0, 2*pi, by=pi/100)) {
#     # 3d rotation around z-axis (rotation matrix)
#     x <- zoom*(cos(i)*as.numeric(cam_coords['x']) - sin(i)*as.numeric(cam_coords['y']))
#     y <- zoom*(sin(i)*as.numeric(cam_coords['x']) + cos(i)*as.numeric(cam_coords['y']))
#     z <- as.numeric(cam_coords['z'])
# 
#     # creating plot
#     respPlot <- plot_ly(main_comps, x=~PC1, y=~PC2, z=~PC3,
#                         type="scatter3d", mode="markers", marker=list(size=4),
#                         color=dfGdsc$pred, colors=c("#FF0000", "#2565AE"))
#     respPlot <- respPlot %>%
#         layout(
#             legend=list(font=list(size=20), itemsizing="constant", orientation="h", x=0.4, y=0),
#             font=list(size=14),
#             scene=list(
#                 xaxis=list(title="PC1 (12.8%)",
#                            zerolinecolor="white",
#                            backgroundcolor="rgb(200, 200, 230)",
#                            showbackground=T,
#                            gridcolor="white"),
#                 yaxis=list(title="PC2 (13.4%)",
#                            zerolinecolor="white",
#                            backgroundcolor="rgb(230, 200, 230)",
#                            showbackground=T,
#                            gridcolor="white"),
#                 zaxis=list(title="PC3 (10.6%)",
#                            zerolinecolor="white",
#                            backgroundcolor="rgb(230, 230, 200)",
#                            showbackground=T,
#                            gridcolor="white"),
#                 camera=list(
#                            eye=list(x=x, y=y, z=z),
#                            center = list(x=0, y=0, z=0)
#                     )
#             ))
#     # saving
#     orca(respPlot, paste("gif/resp", gifnum, ".png", sep=""))
#     gifnum <- gifnum + 1
# }



# fviz plots
# therapeutic response clusters
# responseGroupPlot <- fviz_pca_ind(pc, geom="point", habillage=dfGdsc$pred,
#                                   addEllipses=T, ellipse.level=0.95,
#                                   alpha.ind=0.75, pointsize=1,
#                                   legend.title="Response\ncluster",
#                                   palette=c("#FF0000", "#2565AE")) +
#     theme_bw() +
#     theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
#     theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
#     theme(legend.text=element_text(size=12), legend.title=element_text(size=15)) +
#     ggtitle("")
# responseGroupPlot


# GRID (NEED TO RUN GDSC2.R AND PCA.R TO GET PLOT OBJECTS) ============
# cowResp <- plot_grid(responseGroupPlot, basic,
#                      ncol=2, nrow=1, align="h",
#                      labels=c("A", "B"))
# cowHeatmap <- plot_grid(as.grob(heatmapPlot), nrow=1, labels=c("C"))
# cowCombined <- plot_grid(cowResp, cowHeatmap, nrow=2, ncol=1)
# cowCombined
