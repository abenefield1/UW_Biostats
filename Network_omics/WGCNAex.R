##########
## WGCNA
##########

## Acknowledgement: 
## This code is based on an example from Mike Inouye (UQ), which was primarily 
## written by Scott Ritchie (UQ). 

## Install the package
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("impute", "preprocessCore", "GO.db", "AnnotationDbi"))
# install.packages("WGCNA")

rm(list=ls())

## read the data
read.matrix <- function(file) {
  df <- read.table(file, header=TRUE, row.names=1, sep="\t", 
                   check.names=FALSE)
  mat <- as.matrix(df)
  # Avoid having numeric IDs
  rownames(mat) <- paste0("Probe_", rownames(mat))
  colnames(mat) <- paste0("Sample_", colnames(mat))
  return(mat)
}

## NOTE: Update this path and make sure that you have saved the data 
setwd("/Users/Amy\ 1/Desktop/UW_Biostats/Network_omics")
liver_ge <- read.matrix("liver_expression_head2500.txt")


## Removes probes with >5% of observations missing
filter_probes <- function(x) {
  # How many samples are missing for each probe?
  nMissing <- apply(x, 2, function(probe) {
    sum(is.na(probe))
  })
  missingness <- nMissing/nrow(x)
  return(x[, missingness <= 0.05])
}


## Removes samples that failed > 5% of their assays
filter_samples <- function(x) {
  # How many samples are missing for each probe?
  nMissing <- apply(x, 1, function(sample) {
    sum(is.na(sample))
  })
  missingness <- nMissing/ncol(x)
  return(x[missingness <= 0.05,])
}

liver_ge_fil <- filter_samples(filter_probes(liver_ge))


## Impute the remaining missing values
library(impute) 

if (any(is.na(liver_ge_fil)))
  liver_ge_fil_imp <- impute.knn(liver_ge_fil)$data


anyNA(liver_ge_fil)
anyNA(liver_ge_fil_imp)
dim(liver_ge_fil_imp_vary)

## Chose the top varying probes
most_varying <- function(ge, topN=1000) {
  standard_deviation <- apply(ge, 1, sd)
  sorted <- sort(standard_deviation, decreasing=TRUE)
  sorted_names <- names(sorted)
  topN_names <- sorted_names[seq_len(topN)]
  return(topN_names)
}

top_varying <- most_varying(liver_ge_fil_imp)
liver_ge_fil_imp_vary <- liver_ge_fil_imp[top_varying,]


## Network inference
library(WGCNA)
calculate_coexpression <- function(ge) {
  coexpression <- cor(t(ge), method="pearson")
}

infer_network <- function(coexpression) {
  # First pick the soft threshold to use to define the interaction network
  sft <- WGCNA::pickSoftThreshold(abs(coexpression), dataIsExpr=FALSE)
  if (is.na(sft$powerEstimate)) {
    sft$powerEstimate <- 1
    warning("Could not satisfy the scale-free topology criterion")
  }
  network <- abs(coexpression)^sft$powerEstimate
}

liver_coexpression <- calculate_coexpression(liver_ge_fil_imp_vary)
liver_network <- infer_network(liver_coexpression)
dim(liver_coexpression)
## Power of about 7 is best here

## Identify tightly coexpressed modules in the network 
detect_modules <- function(ge, network) {
  # Calculate the distance between probes based on their topological similarity:
  # i.e. the strength of their coexpression as well as the similarity of their
  # patterns of coexpression to all other probes
  probe_dist <- WGCNA::TOMdist(network)
  dimnames(probe_dist) <- dimnames(network)
  
  # Hierarchically cluster based on this distance metric
  dendro <- hclust(as.dist(probe_dist), method="average")
  
  # Detect modules. `cutreeDynamic` is a function in the `dynamicTreeCut`
  # package, which is loaded in by the `WGCNA` package.
  module_labels <- cutreeDynamic(dendro, distM = probe_dist)
  names(module_labels) <- colnames(network)
  
  # Merge similar modules 
  merged <- mergeCloseModules(t(ge), module_labels, relabel=TRUE)
  module_labels <- merged$colors
  
  return(module_labels)
}

liver_modules <- detect_modules(liver_ge_fil_imp_vary, liver_network)
## Each probe is now assigned a (numeric) module label:
table(liver_modules)


## Save the data for analysis using NetRep
#dir.create("temp_data")

write.matrix <- function(x, file) {
  write.csv(x, file, quote=FALSE)
}

write.vector <- function(x, file, col.names) {
  column_matrix <- t(t(x))
  colnames(column_matrix) <- col.names
  write.csv(column_matrix, file, quote=FALSE)
}


## NetRep requires the probes to be columns, so we will transpose when saving
write.matrix(t(liver_ge_fil_imp_vary), "temp_data/liver_expression.csv")
write.matrix(liver_coexpression, "temp_data/liver_coexpression.csv")
write.matrix(liver_network, "temp_data/liver_network.csv")
write.vector(liver_modules, "temp_data/liver_modules.csv", col.names="Module")


## Phenotype association analysis (using NetRep, based on eigengenes)
MEs <- moduleEigengenes(t(liver_ge_fil_imp_vary), liver_modules)

names(MEs$varExplained) <- colnames(MEs$eigengenes)
MEs$varExplained

## Add sample names to the eigengenes 
rownames(MEs$eigengenes) <- colnames(liver_ge_fil_imp_vary)

pheno<-read.table("pheno.txt", header=TRUE)
data<-cbind(pheno,MEs$eigengenes)

plot(data$ME6,data$pheno,xlab="module summary expression",ylab="phenotype levels",pch=20)
lmfit <- lm(pheno ~ ME6, data=data)
summary(lmfit)
abline(a=lmfit$coefficients[[1]], b=lmfit$coefficients[[2]], col="red", lwd=2)

