setwd(dirname(parent.frame(2)$ofile))
library(HTSet)
library(readxl)
library(tidyverse)
library(KEGGREST)

enzyme <- read.csv("../data-raw/kotable-enzyme-filtered-absolute.tsv", 
                     header = T, sep = "\t", comment.char = "") 
orthology <- read.csv("../data-raw/kotable-filtered-absolute.tsv", 
                      header = T, sep = "\t", comment.char = "") 
# module <- read.csv("../data-raw/kotable-module-filtered-absolute.tsv", 
#                      header = T, sep = "\t", comment.char = "") 
# pathwayL2 <- read.csv("../data-raw/kotable-pathway-L2-filtered-absolute.tsv", 
#                         header = T, sep = "\t", comment.char = "") 
# pathwayL3 <- read.csv("../data-raw/kotable-pathway-L3-filtered-absolute.tsv", 
#                         header = T, sep = "\t", comment.char = "") 
pdata <- readRDS("taxa.rds")$pdata 

f_edata <- function(x) {
    edata <- column_to_rownames(x, colnames(x)[1])
    colnames(edata) <- sub("X", "", colnames(edata))
    edata <- as.matrix(edata)[, rownames(pdata)]
    return(edata)
}
f_fdata <- function(x) {
    fdata <- data.frame(category = rownames(x))
    rownames(fdata) <- rownames(x)
    return(fdata)
}

#get KO fdata
orthology_info <- keggList("ko")
names(orthology_info) <- sub("ko:", "", names(orthology_info))
# remove ids that don't have data available in db
both <- intersect(orthology$X.ko, names(orthology_info))
# get detailed information of KO
orthology_info <- orthology_info[both]
# generate fdata table
orthology_fdata <- data.frame(
    category = orthology_info
    # module = keggLink("module", "ko")
)
# add rownames for fdata
rownames(orthology_fdata) <- names(orthology_info)
# subset edata, remove ids that don't have data available in db
orthology_edata <- f_edata(orthology)[both, ]

enzymeSet <- HTSet(edata = f_edata(enzyme), fdata = f_fdata(f_edata(enzyme)), pdata = pdata)
orthologySet <- HTSet(edata = orthology_edata, fdata = orthology_fdata, pdata = pdata)
# moduleSet <- HTSet(edata = f_edata(module), fdata = f_fdata(f_edata(module)), pdata = pdata)
# pathwayL2Set <- HTSet(edata = f_edata(pathwayL2), fdata = f_fdata(f_edata(pathwayL2)), pdata = pdata)
# pathwayL3Set <- HTSet(edata = f_edata(pathwayL3), fdata = f_fdata(f_edata(pathwayL3)), pdata = pdata)

functions <- list(enzyme = enzymeSet, ko = orthologySet)
saveRDS(functions, "functions.rds")