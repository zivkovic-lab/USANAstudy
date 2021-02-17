setwd(dirname(parent.frame(2)$ofile))
pkgs <- c('tidyverse', 'data.table')
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}


## prepare the files
path <- "../raw-data/Food Pro Data Export to Mac"
files <- list.files(path, pattern = "Visit", full.names = TRUE)


## read 140 files as table and store the tables in a list
nutrients <- vector("list", length = length(files))
ids <- vector("character", length = length(files))
for(i in seq_along(files)){
    file <- files[i]
    id <- basename(files[i]) %>% 
        sub("\\.txt", "", .) %>% 
        str_split(" ") %>%
        unlist() %>%
        str_c(., collapse = "")
    raw_data <- fread(file, nrows=254, skip = "Multi-Column", header = T)
    data <- rbind(raw_data[,1:2], raw_data[,4:5])[-1,]
    ids[i] <- id
    nutrients[[i]] <- data
}

## assign the names of list
names(nutrients) <- ids


saveRDS(nutrients, "diet.rds")
data <- readRDS("diet.rds")