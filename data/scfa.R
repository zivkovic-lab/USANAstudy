setwd(dirname(parent.frame(2)$ofile))
library(HTSet)
library(readxl)
library(tidyverse)

scfa <- read.csv("../raw-data/scfa_concentration_normalized.csv", header = T) 
pdata <- readRDS('taxa.rds')$pdata %>%
    cbind(scfa[, 1:3])
edata <- scfa %>%
    select(ends_with("kg")) %>%
    as.matrix() %>%
    t()
colnames(edata) <- rownames(pdata)
fdata <- data.frame(
    rowname = sub(".mmol.per.kg", '', rownames(edata)),
    unit = sub(".+.mmol.per.kg", 'mmol.per.kg', rownames(edata))
) %>%
    column_to_rownames()
rownames(edata) <- rownames(fdata)
scfa <- HTSet(edata, fdata, pdata)
saveRDS(scfa, 'scfa.rds')
