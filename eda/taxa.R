setwd(dirname(parent.frame(2)$ofile))

# packages ----------------------------------------------------------------

pkgs <- c("tidyverse", "HTSet", "DESeq2", "ggsci", "pheatmap", "edgeR", "zinbwave")
for(pkg in pkgs){
            suppressPackageStartupMessages(library(pkg, character.only = T))
}

taxa <- readRDS("../data/taxa.rds")

# edgeR -------------------------------------------------------------------

fitTaxa <- function(taxaLevel, engine) {
    # taxaLevel can be "kingdom", "phylum", "class", "order", 
    # "family","genus", "species", "strain"
    genus = summarize_feature(taxa, taxaLevel)
    genus = subset_samples(genus, genus$pdata$timepoint %in% c("pre", "post"))
    genus$pdata$timepoint = droplevels(genus$pdata$timepoint)
    genus$pdata$trt = factor(genus$pdata$trt)
    genus = genus[rowSums(genus$edata) > 40,]
    # genus$pdata$trt = factor(genus$pdata$trt, levels = c('B', 'A'))
    design = model.matrix(~ trt * timepoint + subject, data = genus$pdata)
    fit = model_fit(genus, design = design, engine = engine, coef = "trtB:timepointpost")
    return(fit)
}

# zinbwave ----------------------------------------------------------------

# library(DESeq2)
# genus = summarize_feature(taxa, 'genus')
# genus = subset_samples(genus, genus$pdata$timepoint %in% c("pre", "post"))
# genus$pdata$timepoint = droplevels(genus$pdata$timepoint)
# genus$pdata$trt = factor(genus$pdata$trt)
# 
# design = model.matrix(~ trt * timepoint + subject, data = genus$pdata)
# dds <- DESeq2::DESeqDataSetFromMatrix(genus$edata, genus$pdata, design)
# dds <- DESeq(dds, sfType="poscounts", useT=TRUE, minmu=1e-6)
# res <- lfcShrink(dds, coef = 'trtB.timepointpost',
#                  type = "normal")
# head(res)
# or NB regression
# res2 = model_fit(genus, design = design, engine = "edgeR", coef = "trtB:timepointpost")
# bad model

# zinb <- data.frame(
#     count = genus$edata['Klebsiella',],
#     trt = genus$pdata$trt,
#     tp = genus$pdata$timepoint,
#     subj = genus$pdata$subject
# )
# ggplot(zinb, aes(count, fill = trt)) +
#     geom_histogram() +
#     scale_x_log10() +
#     facet_grid(trt ~ ., margins=TRUE, scales="free_y")
# m1 <- zeroinfl(count ~ trt*tp+subj | subj,
#                data = zinb, dist = "negbin")
# summary(m1)
# m2 <- glm.nb(count ~ trt*tp+subj, data = zinb)
# vuong(m1, m2)
# nb fits better!!!

# Visulization ------------------------------------------------------------

boxplot <- function(feature, taxaLevel, timepoint = c("pre", "post")){
    # FIX ME! input of feature and timepoint
    dataset = summarize_feature(taxa, taxaLevel)
    dataset <- subset_samples(dataset, dataset$pdata$timepoint %in% timepoint)
    dataset$pdata$timepoint <- droplevels(dataset$pdata$timepoint)
    HTSet::plot_boxplot(dataset, 
                        feature = feature, 
                        x = 'timepoint', cols = "trt", 
                        line_by = "subject", color_by = "subject") +
        scale_y_log10()
}

violinplot <- function(feature, taxaLevel, timepoint = c("pre", "post")) {
    dataset = summarize_feature(taxa, taxaLevel)
    dataset <- subset_samples(dataset, dataset$pdata$timepoint %in% timepoint)
    dataset$pdata$timepoint <- droplevels(dataset$pdata$timepoint)
    df <- data.frame(
        counts = dataset$edata[feature,],
        trt = dataset$pdata$trt,
        tp = dataset$pdata$timepoint,
        suj = dataset$pdata$subject
    ) 
    ggplot(df, aes(tp, counts)) +
        geom_violin() +
        geom_point(aes(color = suj)) +
        geom_line(aes(group = suj, color = suj)) +
        scale_y_log10() +
        facet_wrap(~ trt)
}

genus = summarize_feature(taxa, 'genus')
genus = subset_samples(genus, genus$pdata$timepoint %in% c("pre", "post"))
genus$pdata$timepoint = droplevels(genus$pdata$timepoint)
genus$pdata$trt = factor(genus$pdata$trt)
abundance <- apply(genus$edata, 2, function(x){x/sum(x)*100})
df <- abundance[c('Bifidobacterium', 'Lactobacillus', 'Akkermansia'),] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("id") %>%
    cbind(genus$pdata[, c('trt', 'timepoint', 'subject')]) %>%
    pivot_longer(cols = 2:4, names_to = 'genus', values_to = "abundance")
ggplot(df, aes(x = genus, y = abundance, color = trt, shape = timepoint)) +
    geom_jitter(width = 0.1)+
    scale_y_log10() +
    theme_bw()
