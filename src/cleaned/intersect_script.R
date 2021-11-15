library(tidyverse)
library(dplyr)
install.packages("sets")

genes <- c('BLM', 'MET', 'LMNA', 'RPL22', 'SMAD3', 'WWTR1', 'UBR5', 'ITGAV', 'IGF2BP2', 'ARID1A', 'MYB', 'VHL', 'NONO', 'EPAS1', 'MED12', 'EXT1', 'CCND1', 'LPP', 'EZH2', 'BCL2', 'RAD21', 'PRKCB', 'EGFR', 'BCL9L', 'RBM10', 'RHOH', 'NACA', 'LATS2', 'CLIP1', 'SDC4', 'CXCR4')
census_data <- read_csv("cancer_gene_census.csv")


dataf1 <- subset(census_data, census_data$`Gene Symbol` %in% genes)
