setwd(dirname(parent.frame(2)$ofile))

pkgs=c('tidyverse', 'HTSet', 'readxl')
for(pkg in pkgs){
            suppressPackageStartupMessages(library(pkg, character.only=TRUE))
}

pdata = read_xlsx(
            "../data-raw/USANA Study Clinical Meta Data Sheet.xlsx", 
            range = "A1:H81",
            sheet = 1
) %>%
            mutate_at(vars(subject), function(x)as.character(x)) %>%
            mutate_at(vars(timepoint), function(x)factor(x, levels = c("pre", "post"))) %>%
            mutate(rowname = paste0(subject, `timepoint coded`)) %>%
            column_to_rownames("rowname") 
glucose = read_xlsx(
            "../data-raw/USANA Study Clinical Meta Data Sheet.xlsx", 
            range = "A1:G21",
            sheet = 3
) %>%
            select("Glucose", "A", "B", "C", "D") %>%
            rename(subject = Glucose) %>%
            pivot_longer(
                        cols = c("A", "B", "C", "D"), 
                        names_to = "timepoint", 
                        values_to = "glucose"
            ) %>%
            transmute(
                        rowname = paste0(subject, timepoint),
                        glucose = glucose
            ) %>%
            column_to_rownames("rowname") %>%
            t()

insulin = read_xlsx(
            "../data-raw/USANA Study Clinical Meta Data Sheet.xlsx", 
            range = "A1:E21",
            sheet = 4
) %>%
            rename(subject = Insulin) %>%
            pivot_longer(
                        cols = c("A", "B", "C", "D"), 
                        names_to = "timepoint", 
                        values_to = "insulin"
            ) %>%
            transmute(
                        rowname = paste0(subject, timepoint),
                        insulin = insulin
            ) %>%
            column_to_rownames("rowname") %>%
            t()

lipidA = read_xlsx(
            "../data-raw/USANA Study Clinical Meta Data Sheet.xlsx", 
            range = "A1:G21",
            sheet = 11
) %>%
            rename(subject = ...1) %>%
            mutate(timepoint = rep("A", nrow(.)))
lipidB = read_xlsx(
            "../data-raw/USANA Study Clinical Meta Data Sheet.xlsx", 
            range = "A1:G21",
            sheet = 12
) %>%
            rename(subject = ...1) %>%
            mutate(timepoint = rep("B", nrow(.)))
lipidC = read_xlsx(
            "../data-raw/USANA Study Clinical Meta Data Sheet.xlsx", 
            range = "A1:G21",
            sheet = 13
) %>%
            rename(subject = ...1) %>%
            mutate(timepoint = rep("C", nrow(.)))
lipidD = read_xlsx(
            "../data-raw/USANA Study Clinical Meta Data Sheet.xlsx", 
            range = "A1:G21",
            sheet = 14
) %>%
            rename(subject = ...1) %>%
            mutate(timepoint = rep("D", nrow(.)))
lipid = bind_rows(lipidA, lipidB, lipidC, lipidD) %>%
            mutate(rowname = paste0(subject, timepoint)) %>%
            select(-c("subject", "timepoint")) %>%
            column_to_rownames("rowname") %>%
            t()
lipid = lipid[, rownames(pdata)]

# identical(colnames(insulin), colnames(lipid))
# identical(colnames(glucose), colnames(lipid))

edata = rbind(glucose, insulin) %>%
            rbind(lipid) 

fdata = data.frame(
            rowname = rownames(edata),
            unit = c("mg/dL", "uU/mL", "mg/dL", "mg/dL", "mg/dL", "", "mg/dL", "mg/dL")
) %>%
            column_to_rownames("rowname") 

lipids = HTSet(edata, fdata, pdata)

saveRDS(lipids, "lipids.rds")