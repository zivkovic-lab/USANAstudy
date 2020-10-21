setwd(dirname(parent.frame(2)$ofile))
pkgs <- c('tidyverse', 'ggplot2', 'HTSet', 'phyloseq', "ggsci")
for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = T))
}
taxa <- readRDS("../data/taxa.rds")

toPhyloseq <- function(taxaLevel){
    dataset <- summarize_feature(taxa, taxaLevel)
    otu <- otu_table(dataset$edata, taxa_are_rows = TRUE)
    sam <- dataset$pdata %>%
        rownames_to_column() %>%
        mutate(
            bmi = ifelse(`BMI (kg/m^2)` < 18.5, 'underweight', 
                         ifelse(`BMI (kg/m^2)` < 24.9, 'normal',
                                ifelse(`BMI (kg/m^2)` < 29.9, 'overweight', 'obese'))),
            age2 = ifelse(age>=30, '>=30', '<30')
        ) %>%
        mutate(
            bmi = factor(bmi, levels = c('underweight', 'normal', 'overweight', 'obese')),
            age2 = factor(age2),
            annotation = factor(annotation, levels = c('baseline', 'A _2_wks', 'A _4_wks', 
                                                       'washout_2_wks', 'washout_4_wks',
                                                       'B _2_wks', 'B _4_wks'))
        ) %>%
        column_to_rownames() %>%
        sample_data()
    tax <- rownames(dataset$edata)
    names(tax) <- rownames(dataset$edata)
    tax <- as.matrix(tax) %>% tax_table()
    phyloseq(otu, tax, sam)
}
taxa <- toPhyloseq('species')

# alpha diversity ---------------------------------------------------------

# Should I normalize my data before alpha-diversity analysis
# No. Generally speaking, the answer is no. Most alpha diversity methods will be 
#most effective when provided with the originally-observed count values.
richness <- estimate_richness(taxa, measures = c("Shannon")) 

timecourse <- function(adjust = T){
    df <- cbind(richness, sample_data(taxa))
    if (adjust) {
        df$annotation <- factor(
            df$annotation, 
            levels = c('baseline', 'A _2_wks', 'A _4_wks', 'washout_2_wks', 
                       'washout_4_wks', 'B _2_wks', 'B _4_wks')
        )
        x = "annotation"
        title = "Alpha Diversity Over Time (adjusted)"
    } else {
        x = "timepoint_code"
        title = "Alpha Diversity Over Time"
    }
    ggplot(df, aes_string(x = x, y = 'Shannon')) +
        stat_boxplot(geom = "errorbar", width = 0.5) +
        geom_boxplot(size = 0.5) +
        geom_point(aes(color = subject)) +
        geom_line(aes(color = subject, group = subject)) +
        labs(title = title) +
        theme_bw()
}
anovaTable <- function(contrast = c("pre", "post")){
    df <- cbind(richness, sample_data(taxa)) %>%
        filter(trt %in% c('A', 'B')) %>%
        filter(timepoint %in% contrast) %>%
        mutate(
            trt = droplevels(trt),
            timepoint = droplevels(timepoint)
        )
    model <- lm(Shannon~trt+timepoint+trt:timepoint+subject, data = df)
    anova(model)
}
boxplot <- function(tp = c("pre", "post")){
    df <- cbind(richness, sample_data(taxa)) %>%
        filter(trt %in% c('A', 'B')) %>%
        filter(timepoint %in% tp) %>%
        mutate(
            trt = droplevels(trt),
            timepoint = droplevels(timepoint)
        )
    ggplot(df, aes(timepoint, Shannon)) +
        geom_boxplot() +
        geom_point(aes(color = subject)) +
        geom_line(aes(color = subject, group = subject)) +
        theme_bw() +
        facet_wrap(~trt)
}

# beta diversity ----------------------------------------------------------

betaDiversity <- function(dist_methods = "bray", ordinate_method = "PCoA"){
    taxa_relative <- transform_sample_counts(taxa, function(x) x / sum(x) )
    taxa_relative <- filter_taxa(taxa_relative, function(x) mean(x) > 1e-5, TRUE)
    # batch <- taxa@sam_data$subject %>% factor()
    # ex_b_edata <- removeBatchEffect(log(taxa_relative@otu_table@.Data+1e-10), batch) %>%
    #     exp()
    # taxa_relative@otu_table@.Data <- ex_b_edata
    # Calculate distance matrix
    iDist <- distance(taxa_relative, method=dist_methods)
    # Calculate ordination
    ord  <- ordinate(taxa_relative, ordinate_method, distance=iDist)
}
# plot_scree(betaDiversity(), "Scree plot for species, Bray-Curtis/PCoA")

plotBetaDiversity <- function(ord, color_by = "annotation", tp = LETTERS[1:7]){
    taxa_relative <- transform_sample_counts(taxa, function(x) x / sum(x) )
    taxa_relative <- filter_taxa(taxa_relative, function(x) mean(x) > 1e-5, TRUE)
    taxa_relative <- subset_samples(taxa_relative, timepoint_code %in% tp)
    plot_ordination(taxa_relative, ord, "samples", color=color_by) +
        geom_point(size = 2) + 
        geom_line(aes(group = subject)) +
        # stat_ellipse() + 
        theme_bw()
    # ggplot(df, aes_string(colnames(df)[1], colnames(df)[2], color = color_by))+
    #     geom_point()+
    #     theme_bw()
}    
