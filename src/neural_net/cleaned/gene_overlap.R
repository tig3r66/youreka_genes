library(dplyr)

blist <- read_csv("boruta.csv", col_names=F)
reflist <- read_csv("cancer_gene_census.csv", col_names = F)

intersect(blist$X1, reflist$X1)
