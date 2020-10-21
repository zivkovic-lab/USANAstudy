setwd(dirname(parent.frame(2)$ofile))
library(HTSet)
library(readxl)
library(tidyverse)

edata <- read.csv("../data-raw/taxatable-filtered-absolute.tsv", 
                    header = T, sep = "\t", comment.char = "") 
taxa <- edata$X.tax
fdata <- data.frame(
            otuID = paste0("otu", (1:length(taxa))),
            taxa = taxa
) %>% 
            separate(
                        taxa, 
                        into = c("kingdom", "phylum", "class", "order", "family","genus", "species", "strain"),
                        sep = ";" ,
                        convert = T
            ) %>%
            mutate_at(
                        c("kingdom", "phylum", "class", "order", "family","genus", "species", "strain"),
                        ~sub(".__", "", .)
            ) %>%
            column_to_rownames("otuID")
edata <- select(edata, !X.tax) 
rownames(edata) <- rownames(fdata)
colnames(edata) <- sub("X", "", colnames(edata))

metadata = read_xlsx(
            "../data-raw/USANA Study Clinical Meta Data Sheet.xlsx", 
            range = "A1:K81",
            sheet = 1
)
bmi = read_xlsx(
    "../data-raw/USANA Study Clinical Meta Data Sheet.xlsx", 
    range = "A1:I101",
    sheet = 2
) %>%
    dplyr::rename('subject' = '...1') %>%
    mutate(subject = rep(subject[seq(1,100, 5)], each = 5)) %>%
    pivot_longer(cols = 3:9, names_to = "timepoint", values_to = "value") %>%
    pivot_wider(id_cols = c('subject', 'timepoint'), names_from = "Timepoint", values_from = "value") 

pdata <- data.frame(
    id = colnames(edata),
    subject = sub("[A-Z]", "", colnames(edata)),
    timepoint_code = rep(LETTERS[1:7], 20),
    first_trt = rep(metadata$treatment[seq(1, 80, 4)], each = 7),
    second_trt = rep(metadata$treatment[seq(3, 80, 4)], each = 7),
    gender = rep(metadata$gender[seq(1, 80, 4)], each = 7),
    age = rep(metadata$age[seq(1, 80, 4)], each = 7)
) %>%
    cbind(bmi[,3:7]) %>%
    mutate(
        annotation = apply(., 1, function(x){
            switch(
                x[3],
                "A" = "baseline",
                "B" = paste(x[4], "_2_wks"),
                "C" = paste(x[4], "_4_wks"),
                "D" = "washout_2_wks",
                "E" = "washout_4_wks",
                "F" = paste(x[5], "_2_wks"),
                "G" = paste(x[5], "_4_wks")
            )
        })
    ) %>%
    mutate(
        timepoint = rep(c("pre", "mid", "post", "mid", "pre", "mid", "post"), 20),
        trt = ifelse(.$timepoint_code %in% c('A', 'B', 'C'), first_trt, 
                     ifelse(.$timepoint_code == 'D', 'washout', second_trt))
    ) %>%
    mutate(
        timepoint = factor(timepoint, levels = c("pre", "mid", "post")),
        trt = factor(trt)
    ) %>%
    select(!(first_trt:second_trt)) %>%
    column_to_rownames("id")
edata <- as.matrix(edata)

taxa_count = HTSet(edata = edata, fdata = fdata, pdata = pdata)
saveRDS(taxa_count, file = "taxa.rds")