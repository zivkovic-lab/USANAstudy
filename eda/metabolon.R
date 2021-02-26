setwd(dirname(parent.frame(2)$ofile))

pkgs <- c("tidyverse", "HTSet", "ggsci", "limma", "factoextra", "pheatmap")
for(pkg in pkgs){
    suppressPackageStartupMessages(library(pkg, character.only = T))
}

metabolites <- readRDS("../data/metabolites.rds")

# The normality of the data
plotHist <- function(transform = "log"){
    if (transform == "log") {
        abund <- apply(metabolites$edata, 1, function(x) median(scale(log(x+1))))
    } else {
        abund <- apply(metabolites$edata, 1, function(x) median(scale(x)))
    }
    data.frame(
        chem = featureNames(metabolites),
        abund = abund
    ) %>% 
        ggplot(aes(x = abund, y = ..density..)) +
        geom_histogram(color = "white", bins = 25, fill = "steelblue") +
        geom_density()
}
plotQQ <- function(transform = "log"){
    if (transform == "log") {
        abund <- apply(metabolites$edata, 1, function(x) median(scale(log(x+1))))
    } else {
        abund <- apply(metabolites$edata, 1, function(x) median(scale(x)))
    }
    qqnorm(abund, pch = 1, frame = FALSE)
    qqline(abund, col = "steelblue", lwd = 2)
}

# -------------------------------------------------------------------------
# Linear Model
# -------------------------------------------------------------------------
# remove the chem with variance = 0:
metabolites <- metabolites[apply(metabolites$edata, 1, var) != 0,]
calDE <- function(transform = "log"){
    if (transform == "log") {
        edata <- log(metabolites$edata)
    } else {
        edata <- metabolites$edata
    }
    model <- model.matrix(~tp*trt+subject_id, data = metabolites$pdata)
    print(colnames(model))
    fit <- lmFit(edata, model) %>% eBayes()
    topTable(fit, coef = "tppost:trtB", number = Inf, sort.by = "none")
}
plotVolc <- function(topTable){
    ggplot(topTable, aes(logFC, -log(P.Value))) +
        geom_point()+
        geom_hline(yintercept = -log(0.05)) +
        theme_bw()
}
plotHistP <- function(topTable){
    ggplot(topTable, aes(P.Value)) +
        geom_histogram(fill = "steelblue", color = "black", binwidth = 0.05)+
        geom_vline(xintercept = 0.05, color = "red", linetype = 2) +
        theme_bw()
}
plotBox <- function(topTable, rows_selected = 1){
    if (identical(rownames(topTable), featureNames(metabolites))) {
        # print(rownames(metabolites$edata)[rows_selected])
        dt <- data.frame(
            Response = metabolites$edata[rows_selected,],
            Timepoint = metabolites$pdata$tp,
            Treatment = metabolites$pdata$trt
        ) 
        print(dt)
            ggplot(dt, aes(x = Timepoint, y = Response))+
            geom_boxplot(aes(fill = Treatment))+
            geom_point(aes(group = Treatment)) +
            facet_wrap(~Treatment)+
            theme_bw()
    } else {
        print("The rownames of DE table and original edata table are not identical.")
    }
}
# -------------------------------------------------------------------------
# PCA and Heatmap
# -------------------------------------------------------------------------
DEgenes <- calDE() %>%
    dplyr::filter(P.Value<0.05) %>%
    rownames()
dat4PCA <- metabolites[DEgenes,]
df <- dat4PCA$edata %>%
    t()
res.pca <- prcomp(df, center = T, scale. = T)
fviz_screeplot(res.pca, addlabels = TRUE)
# Extract the results for variables
var <- get_pca_var(res.pca)
# Graph of variables: default plot
fviz_pca_var(res.pca, col.var = "black")
plotPCA <- function(color_by = "group_name"){
    fviz_pca_ind(res.pca,
                 label = "none", # hide individual labels
                 habillage = dat4PCA$pdata[[color_by]], # color by groups
                 addEllipses = TRUE # Concentration ellipses
    )
}
plotHeatmap <- function(ann_col = "group_name", ann_row = "super_pathway"){
    mat <- scale(t(dat4PCA$edata)) %>% t()
    mat[mat > 2] <- 2
    mat[mat < -2] <- -2
    # Generate annotations for rows and columns
    annotation_col <- data.frame(
        group = dat4PCA$pdata[,ann_col]
    )
    rownames(annotation_col) <- sampleNames(dat4PCA)
    annotation_row <- data.frame(
        pathway = dat4PCA$fdata[,ann_row]
    )
    rownames(annotation_row) <- featureNames(dat4PCA)
    # Specify colors
    ann_colors = list(
        group = pal_igv()(length(unique(annotation_col$group))),
        pathway = pal_ucscgb()(length(unique(annotation_row$pathway)))
    )
    names(ann_colors$group) <- unique(annotation_col$group)
    names(ann_colors$pathway) <- unique(annotation_row$pathway)
    pheatmap(mat, show_colnames = FALSE, 
             annotation_col = annotation_col, 
             annotation_row = annotation_row, 
             annotation_colors = ann_colors)
}
# -------------------------------------------------------------------------
# Pathway Enrichment
# -------------------------------------------------------------------------
# TBD