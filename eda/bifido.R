setwd(dirname(parent.frame(2)$ofile))
paths <- list.files("../data-raw", pattern="annotations", full.names = TRUE)
# sample.names <- sapply(strsplit(basename(paths), "_"), `[`, 2)
# genes <- vector("list", length(sample.names))
# names(genes) <- sample.names
# for (i in 1:length(paths)){
#     genes[[i]] <- read.csv(paths[i], header = T, sep = "\t", skip = 3)
# }

colname <- read.csv(paths[2], header = T, sep = "\t", skip = 3) %>%
    colnames()
path <- list.files("../data-raw", pattern="201A", full.names = TRUE)
gene1 <- read.csv(path, header = F, sep = "\t", comment.char = "#")
colnames(gene1) <- colname
koTable <- gene1 %>%
    group_by(KEGG_ko) %>%
    summarise(count = n())