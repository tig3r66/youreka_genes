library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)    # heatmap plot
library(viridis)     # for heatmap colours


# READING CELL LINE DATA AND APPENDING COSMIC IDENTIFIER ================
## GDSC data
dfGdsc <- read.csv("temp1GDSC.csv")
dfGdsc <- dfGdsc[, which(colMeans(!is.na(dfGdsc)) > 0.8)]
dfGdsc <- na.omit(dfGdsc)
dfGdsc$Cell_name <- gsub("\\.", "-", dfGdsc$Cell_name)
## cell lines with cosmic identifiers
meta_data <- read.csv("Cell_Lines_Details.csv")
meta_data <- (meta_data %>% rename(Cell_name=Sample.Name))
## merging on "Cell_name"
dfGdsc <- merge(dfGdsc, meta_data, "Cell_name")
rownames(dfGdsc) <- dfGdsc$COSMIC.identifier
dfGdsc$Cell_name <- NULL
dfGdsc$COSMIC.identifier <- NULL
## normalizing data
dfGdsc <- as.data.frame(t(scale(t(scale(dfGdsc)))))


# HEATMAP ANALYSIS ================
## reading in data
dfDrugFunc <- read.csv("drug_func.csv")
dfDrugFunc <- dfDrugFunc[!duplicated(dfDrugFunc$DRUG_NAME),]
dfDrugFunc$DRUG_NAME <- gsub(" ", "\\.", dfDrugFunc$DRUG_NAME)
dfDrugFunc <- filter(dfDrugFunc, DRUG_NAME %in% colnames(dfGdsc))
## drug annotations
colanno <- subset(dfDrugFunc, select=c(DRUG_NAME, PATHWAY_NAME))
colanno$PATHWAY_NAME <- as.factor(colanno$PATHWAY_NAME)
rownames(colanno) <- colanno$DRUG_NAME
colanno$DRUG_NAME <- NULL
## generating Heatmap
pheatmap(dfGdsc, color=inferno(10),
    cluster_rows = F, cluster_cols = T,
    kmeans_k = 2, show_rownames=F)

# CLUSTERING ============
kmeansGdsc <- kmeans(dfGdsc, center=2, iter.max = 10, nstart = 25)
class <- as.data.frame(kmeansGdsc$cluster)
class$CELL_LINE_NAME <- rownames(class)
dfGdsc$CELL_LINE_NAME <- rownames(dfGdsc)
data <- merge(dfGdsc, class, "CELL_LINE_NAME")
rownames(dfGdsc) <- dfGdsc$CELL_LINE_NAME
dfGdsc$pred <- dfGdsc$`kmeansGdsc$cluster`
dfGdsc$`kmeansGdsc$cluster` <- NULL
dfGdsc$CELL_LINE_NAME <- NULL

## Merging with expression data
exprs <- read.csv("expressiondata_GDSC.csv")
exprs$CELL_LINE_NAME <- gsub("DATA.", "", exprs$CELL_LINE_NAM)
combined <- merge (class, exprs, "CELL_LINE_NAME")
write.csv(combined, "combined_expression.csv")


# Differential Expression Analysis ============
## Volcano DESEQ2
expression <- read.csv("combined_expression.csv")
identification <- subset(expression, select=c("CELL_LINE_NAME", "kmeansGdsc.cluster"))
identification$COSMIC.identifier <- identification$CELL_LINE_NAME
identification$CELL_LINE_NAME <- NULL

## storing cosmic identifiers
dfGdsc$COSMIC.identifier <- rownames(dfGdsc)
for_use <- merge(dfGdsc, identification, "COSMIC.identifier")
rownames(for_use) <- for_use$COSMIC.identifier
for_use$COSMIC.identifier <- NULL

storage <- data.frame()
foldchange <- data.frame()
for (i in 1:(length(for_use))) {
  x <- filter(for_use, for_use$kmeansGdsc.cluster==1)
  a <- x[,i]
  y <- filter(for_use, for_use$kmeansGdsc.cluster==2)
  b <- y[,i]

  temp1<-broom::tidy(wilcox.test(a,b))
  temp2<-log2(median(a)/median(b))
  storage <- rbind(storage, temp1)
  foldchange <- rbind(foldchange, temp2)
}

results_pval <- as.data.frame(storage)
results_pval$p.adj <- p.adjust(results_pval$p.value, method="fdr")
results_fc <- as.data.frame(foldchange)
tempor <- colnames(for_use)

results <- as.data.frame(cbind(results_pval,results_fc))
results <- as.data.frame(cbind(results,tempor))
write.csv(results, "volcano_data.csv")

results <- read.csv("volcano_data.csv")
results$id <- as.vector(colnames(data[-1]))
results$classifier <- ifelse(results$p.adj <= 0.01, 1, 0)
results$padj <- -log10(results$p.adj)

## creating volcano plot
basic <- ggplot(results, aes(x=statistic, y=padj)) +
  theme_bw() +
  geom_hline(yintercept=2, linetype="dashed", color = "black", size=0.3) +
  geom_vline(xintercept=33301, linetype="dashed", color = "black", size=0.3) +
  geom_point(aes(colour=as.factor(classifier))) +
  theme(legend.position="none") +
  scale_colour_manual(values=c("#2565ae", "#ff0000")) +
  theme(axis.title=element_text(size=15), axis.text=element_text(size=12)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  labs(x="Test statistic", y=bquote('-'*~Log[10]*' P'))
basic


# FUNCTIONAL ANNOTATION ============
temporary <- read.csv("drugs.csv")
temporary$class<- ifelse(temporary$p.adj > 0.05,1,0)
seta <- subset(temporary, temporary$class==0)
setb <- subset(temporary, temporary$class==1)
pathA <- as.data.frame(table(seta$Pathway))
pathB <- as.data.frame(table(setb$Pathway))
pathA$proportion <- pathA$Freq / nrow(seta)
pathB$proportion <- pathB$Freq / nrow(setb)
comb <- merge(pathA, pathB, "Var1")
write.csv(comb, "drug_data.csv")
## gprofiler plot
gprofiler <- read.csv("gProfiler.csv")
gprofiler <- subset(gprofiler, select=c("adjusted_p_value", "term_name"))
gprofiler$logp <- -log10(gprofiler$adjusted_p_value)
ggplot(gprofiler, aes(x = reorder(term_name, logp), y = logp, fill=scale(logp))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_hline(yintercept=1.30103, linetype="dashed", color = "black", size=0.5) +
  scale_fill_gradient2(low="#f1f7fc", mid = "#83CEFF", high ="#2565AE", name="Z-score") +
  guides(fill=guide_colourbar(ticks.colour="black", ticks.linewidth=1.5,
                              frame.colour="black", frame.linewidth=1.5),
         border=element_line(color="black")) +
  labs(x="", y="-Log10(FDR)") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(color="black", size=.1)) +
  theme(legend.text=element_text(size=13),
        legend.title=element_text(size=13, vjust=0.85),
        legend.position="bottom") +
  theme(axis.text=element_text(size=14, color="black"),
        axis.title=element_text(size=14, color="black"))

