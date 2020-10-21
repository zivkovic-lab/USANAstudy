setwd(dirname(parent.frame(2)$ofile))
pkgs <- c('tidyverse', 'ggplot2', 'HTSet', "limma")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

scfa <- readRDS("../data/scfa.rds")

fitScfa <- function(contrast = c("pre", "post")) {
    scfa = subset_samples(scfa, scfa$pdata$trt %in% c('A', 'B'))
    scfa$pdata$trt = droplevels(scfa$pdata$trt)
    scfa = subset_samples(scfa, scfa$pdata$timepoint %in% contrast)
    scfa$pdata$timepoint = droplevels(scfa$pdata$timepoint)
    design = model.matrix(~ trt * timepoint + subject, data = scfa$pdata)
    fit = lmFit(log(scfa$edata+0.1), design = design) %>%
        eBayes() %>%
        topTable(coef = colnames(design)[ncol(design)], number = Inf, sort.by = 'none')
    return(fit)
}

sboxplot <- function(feature, timepoint = c("pre", "post")){
    # feature is the rowname of fit table
    dataset = subset_samples(scfa, scfa$pdata$trt %in% c('A', 'B'))
    dataset$pdata$trt = droplevels(dataset$pdata$trt)
    dataset <- subset_samples(dataset, dataset$pdata$timepoint %in% timepoint)
    dataset$pdata$timepoint <- droplevels(dataset$pdata$timepoint)
    df = dataset$edata %>%
        t() %>%
        cbind(dataset$pdata[c("trt", "timepoint", "subject")])
    ggplot(df, aes_string("timepoint", feature))+
        geom_boxplot()+
        geom_point(aes(color = subject))+
        geom_line(aes(group = subject, color = subject))+
        facet_wrap(~trt)+
        theme_bw()+
        labs(title = feature)
}

sboxplot("Acetic.Acid")
