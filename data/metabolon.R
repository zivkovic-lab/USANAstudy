setwd(dirname(parent.frame(2)$ofile))
library(HTSet)
library(readxl)
library(tidyverse)

file <- "../raw-data/UCDA-03-20MD DATA TABLES.XLSX"
# row = sampleID; col = compounfID
edata <- read_xlsx(file, sheet = 6) %>%
    column_to_rownames("PARENT_SAMPLE_NAME") %>%
    as.matrix() %>%
    t()
    
pdata <- read_xlsx(file, sheet = 3) %>%
    mutate(
        GENDER = factor(GENDER, levels = c("Male", "Female")),
        GROUP_NAME = factor(GROUP_NAME, levels = c("Baseline", "Post-Tx A", "Post Wash", "Post-Tx B")),
        SUBJECT_ID = factor(SUBJECT_ID),
        TP = ifelse(grepl("Baseline", GROUP_NAME, ignore.case = T), "pre", 
                    ifelse(grepl("wash", GROUP_NAME, ignore.case = T), "pre", "post")) %>% factor(levels = c("pre", "post")),
        TREATMENT = factor(TREATMENT),
        TRT = sub(".+ ", "", GROUP_NAME)[seq(2, 80, 2)] %>% rep(each = 2) %>% factor()
    ) %>%
    column_to_rownames("PARENT_SAMPLE_NAME")
fdata <- read_xlsx(file, sheet = 2) %>%
    as.data.frame()
rownames(fdata) <- fdata$CHEM_ID
colnames(pdata) <- tolower(colnames(pdata))
colnames(fdata) <- tolower(colnames(fdata))

metabolites <- HTSet(edata, fdata, pdata)
featureNames(metabolites) <- fdata$chemical_name
saveRDS(metabolites, 'metabolites.rds')
