setwd(dirname(parent.frame(2)$ofile))

# packages ----------------------------------------------------------------

pkgs <- c("tidyverse", "HTSet", "DESeq2", "ggsci", "edgeR", "clusterProfiler", "enrichplot")
for(pkg in pkgs){
    suppressPackageStartupMessages(library(pkg, character.only = T))
}
functions <- readRDS('../data/functions.rds')

diffExp <- function(dataset, engine){
    dataset = subset_samples(dataset, dataset$pdata$timepoint %in% c("pre", "post"))
    dataset = subset_samples(dataset, dataset$pdata$trt %in% c('A', 'B'))
    dataset$pdata$timepoint = droplevels(dataset$pdata$timepoint)
    dataset$pdata$trt = factor(dataset$pdata$trt)
    dataset = dataset[rowSums(dataset$edata) > 200,]
    design = model.matrix(~ trt * timepoint + subject, data = dataset$pdata)
    model_fit(dataset, design, coef = "trtB:timepointpost", engine = engine)
}
#---
# use DESeq2 for ko data, pvalues skew to 1
# use edgeR for ko data, pvalues skew to 1 but better than the results from DESeq2
# rowSums > 200:
# use edgeR for ko data, pvalues are uniform distribution
# use deseq2 for ko data, the distribution of pvalues is better
# for enzyme data, the distribution of pvalues skew to 0 with edgeR, 
# skew to 1 with deseq2
#---
# Pre-post
koRes <- diffExp(functions$ko, 'DESeq2')
enzymeRes <- diffExp(functions$enzyme, 'edgeR')
# Pre-mid
diffExp2 <- function(dataset, engine){
    dataset = subset_samples(dataset, dataset$pdata$timepoint %in% c("pre", "mid"))
    dataset = subset_samples(dataset, dataset$pdata$trt %in% c('A', 'B'))
    dataset$pdata$timepoint = droplevels(dataset$pdata$timepoint)
    dataset$pdata$trt = factor(dataset$pdata$trt)
    dataset = dataset[rowSums(dataset$edata) > 200,]
    design = model.matrix(~ trt * timepoint + subject, data = dataset$pdata)
    model_fit(dataset, design, coef = "trtB:timepointmid", engine = engine)
}
tmp1 <- diffExp2(functions$ko, 'DESeq2')
# tmp2 <- diffExp2(functions$ko, 'edgeR')


genelist <- koRes$results$logFC
names(genelist) <- rownames(koRes$results)
genelist <- base::sort(genelist, decreasing = T)
gseM <- gseMKEGG(
    genelist,
    organism = 'ko',
    pvalueCutoff = 1
)
gseP <- gseKEGG(
    genelist,
    organism = 'ko',
    pvalueCutoff = 1
)


boxplot <- function(dataset, feature, tp){
    dataset <- dataset %>%
        HTSet::subset_samples((dataset$pdata$trt %in% c('A', 'B'))&
                                  (dataset$pdata$timepoint %in% tp)) 
    selected <- dataset[feature,] 
    HTSet::plot_boxplot(selected, 
                        feature = feature, 
                        x = 'timepoint', cols = "trt", 
                        line_by = "subject", color_by = "subject") +
        scale_y_log10()
}

violinplot <- function(dataset, feature, timepoint = c("pre", "post")){
    selected <- dataset[feature,]
    df <- data.frame(
        counts = selected$edata[1,],
        trt = selected$pdata$trt,
        tp = selected$pdata$timepoint,
        suj = selected$pdata$subject
    ) %>%
        filter((tp %in% timepoint)&(trt %in% c('A', 'B'))) %>%
        mutate(
            tr = droplevels(tp),
            trt = droplevels(trt)
        )
    ggplot(df, aes(tp, counts)) +
        geom_violin() +
        geom_point(aes(color = suj)) +
        geom_line(aes(group = suj, color = suj)) +
        scale_y_log10() +
        facet_wrap(~ trt) + 
        theme_bw()
}

save(enzymeRes, koRes, gseM, gseP, file = '../data/functions.rda')
# functions$ko$fdata['K01226',]
# boxplot(functions$enzyme, '4.2.1.3 ACO, acnA; aconitate hydratase', c('pre', 'post'))
# gseaplot2(gseP, geneSetID = 4, title = gseP$Description[4])

# library("pathview")
# pathview(
#     gene.data  = toGeneList(koRes),
#     species = "ko",
#     pathway.id = "ko00620",
# )

abundance <- apply(functions$ko$edata, 2, function(x){x/sum(x)*100})
df <- abundance[c('K13788', 'K00925', 'K17735', 'K00929'),] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("id") %>%
    cbind(functions$ko$pdata[, c('trt', 'timepoint', 'subject')]) %>%
    pivot_longer(cols = 2:5, names_to = 'ko', values_to = "abundance")
ggplot(df, aes(x = ko, y = abundance, color = trt, shape = timepoint)) +
    geom_jitter(width = 0.1)+
    scale_y_log10() +
    theme_bw()