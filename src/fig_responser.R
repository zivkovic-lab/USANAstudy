setwd(dirname(parent.frame(2)$ofile))

# sourcing -----------------------------------------------------------------

source("global.R")

# packages ----------------------------------------------------------------

pkgs <- c("dplyr", "tidyverse", "Metabase", "limma", "ggpubr", "ggsci")
for(pkg in pkgs){
            suppressPackageStartupMessages(library(pkg, character.only = T))
}

# data --------------------------------------------------------------------

lipids <- readRDS("../data/lipids.rds")

trtA <- subset_samples(lipids, lipids$sample_table$treatment == 'A')
trtB <- subset_samples(lipids, lipids$sample_table$treatment == 'B')

preA <- subset_samples(trtA, trtA$sample_table$timepoint == 'pre')
postA <- subset_samples(trtA, trtA$sample_table$timepoint == 'post')
preB <- subset_samples(trtB, trtB$sample_table$timepoint == 'pre')
postB <- subset_samples(trtB, trtB$sample_table$timepoint == 'post')

deltaA <- preA$conc_table - postA$conc_table
deltaB <- preB$conc_table - postB$conc_table

# lm model ----------------------------------------------------------------

# design <- model.matrix(
#             ~ treatment * timepoint + subject, 
#             data = as(lipids$sample_table, "data.frame")
# )
# fit <- lmFit(lipids$conc_table, design) %>% 
#             eBayes() %>%
#             topTable(coef = ncol(design), number = Inf)

features <- rownames(deltaA)
responseA <- matrix(NA, nrow = 8, ncol = 20)
rownames(responseA) <- features
colnames(responseA) <- colnames(deltaA)
for(i in features) {
            if(rowMeans(deltaA)[i] <= 0){responseA[i, ] = deltaA[i,] < 0}
            else {responseA[i, ] = deltaA[i,] > 0 }
}

responseB <- matrix(NA, nrow = 8, ncol = 20)
rownames(responseB) <- features
colnames(responseB) <- colnames(deltaB)
for(i in features) {
            if(rowMeans(deltaB)[i] <= 0){responseB[i, ] = deltaB[i,] < 0}
            else {responseB[i, ] = deltaB[i,] > 0 }
}

# plot with responders ----------------------------------------------------

responsesOnly <- function(x){
            dfA = data.frame(
                        value = trtA$conc_table[x,],
                        subject = trtA$sample_table[, 'subject'],
                        timepoint = trtA$sample_table[, 'timepoint'],
                        response = factor(
                                    rep(responseA[x,], each = 2), 
                                    labels = c("non-responder", "responder")
                        ),
                        treatment = 'A'
            ) %>%
            filter(response == 'responder') 
            
            dfB <- data.frame(
                        value = trtB$conc_table[x,],
                        subject = trtB$sample_table[, 'subject'],
                        timepoint = trtB$sample_table[, 'timepoint'],
                        response = factor(
                                    rep(responseB[x,], each = 2), 
                                    labels = c("non-responder", "responder")
                        ),
                        treatment = 'B'
            ) %>%
            filter(response == 'responder') 
            
            df <- bind_rows(dfA, dfB)
            ## which model should we use?
            model = lm(value ~ timepoint * treatment, data = df)
            pvalue = summary(model)
            
            p <- ggplot(df, aes(x = timepoint, y = value)) +
                        stat_boxplot(geom = "errorbar", width = 0.5)+
                        geom_boxplot(size = 0.5) +
                        geom_point(aes(color = subject)) +
                        geom_line(aes(group = subject, color = subject), size = 0.5) +
                        facet_wrap(~treatment) +
                        labs(
                                    title = paste0(x, '_responders', ' P=', round(pvalue$coefficients[4, 4], digits = 2)),
                                    y = lipids$feature_data[x,1]
                        ) +
                        theme_boxplot() + 
                        scale_color_manual(values = colorRampPalette(pal_lancet()(9))(20))+
                        theme(
                                    legend.text = element_text(size = 7),
                                    legend.key.height = unit(0.5, "cm"),
                                    legend.key.width = unit(0.5, "cm")
                        )
            
            return(p)
}

# plot with non-responders ------------------------------------------------

nonresponsesOnly <- function(x){
            dfA = data.frame(
                        value = trtA$conc_table[x,],
                        subject = trtA$sample_table[, 'subject'],
                        timepoint = trtA$sample_table[, 'timepoint'],
                        response = factor(
                                    rep(responseA[x,], each = 2), 
                                    labels = c("non-responder", "responder")
                        ),
                        treatment = 'A'
            ) %>%
                        filter(response == 'non-responder') 
            
            dfB <- data.frame(
                        value = trtB$conc_table[x,],
                        subject = trtB$sample_table[, 'subject'],
                        timepoint = trtB$sample_table[, 'timepoint'],
                        response = factor(
                                    rep(responseB[x,], each = 2), 
                                    labels = c("non-responder", "responder")
                        ),
                        treatment = 'B'
            ) %>%
                        filter(response == 'non-responder') 
            
            df <- bind_rows(dfA, dfB)
            ## which model should we use?
            model = lm(value ~ timepoint * treatment, data = df)
            pvalue = summary(model)
            
            p <- ggplot(df, aes(x = timepoint, y = value)) +
                        stat_boxplot(geom = "errorbar", width = 0.5)+
                        geom_boxplot(size = 0.5) +
                        geom_point(aes(color = subject)) +
                        geom_line(aes(group = subject, color = subject), size = 0.5) +
                        facet_wrap(~treatment) +
                        labs(
                                    title = paste0(x, '_nonresponders', ' P=', round(pvalue$coefficients[4, 4], digits = 2)),
                                    y = lipids$feature_data[x,1]
                        ) +
                        theme_boxplot() +
                        scale_color_manual(values = colorRampPalette(pal_lancet()(9))(20))+
                        theme(
                                    legend.text = element_text(size = 7),
                                    legend.key.height = unit(0.5, "cm"),
                                    legend.key.width = unit(0.5, "cm")
                        )
            return(p)
}

# save plots --------------------------------------------------------------

features <- rownames(lipids$feature_data)
for (feature in features) {
            ggsave(paste0(feature, "_responders", ".png"), responsesOnly(feature), device = 'png',
                   path = '../img/',
                   scale = 1.5, width = 3, height = 3, units = "in")
            ggsave(paste0(feature, "_nonresponders", ".png"), nonresponsesOnly(feature), device = 'png',
                   path = '../img/',
                   scale = 1.5, width = 3, height = 3, units = "in")
}

ggsave(paste0("chol_HDL_responders", ".png"), responsesOnly(features[6]), device = 'png',
       path = '../img/',
       scale = 1.5, width = 3, height = 3, units = "in")
ggsave(paste0("chol_HDL_nonresponders", ".png"), nonresponsesOnly(features[6]), device = 'png',
       path = '../img/',
       scale = 1.5, width = 3, height = 3, units = "in")

plot_responseA <- function(x){
            # x is one of the features from glucose/insulin/lipid panels
            df = data.frame(
                        value = trtA$conc_table[x,],
                        subject = trtA$sample_table[, 'subject'],
                        timepoint = trtA$sample_table[, 'timepoint'],
                        response = factor(
                                    rep(responseA[x,], each = 2),
                                    labels = c("non-responder", "responder")
                        )
            )
            unit = lipids$feature_data[x ,1]
            p <- ggplot(df, aes(x = timepoint, y = value)) +
                        stat_boxplot(geom = "errorbar", width = 0.5)+
                        geom_boxplot(size = 0.5) +
                        geom_point(aes(color = subject)) +
                        geom_line(aes(group = subject, color = subject), size = 0.5) +
                        facet_wrap(~ response) +
                        labs(
                                    title = paste("Treatment A", x, sep = " "),
                                    y = unit
                        ) +
                        theme_boxplot() +
                        scale_color_manual(values = colorRampPalette(pal_lancet()(9))(20))
                        theme(
                                    legend.text = element_text(size = 7),
                                    legend.key.height = unit(0.5, "cm"),
                                    legend.key.width = unit(0.5, "cm")
                        )
            return(p)
}

features <- rownames(trtA$feature_data)
for (feature in features) {
            ggsave(paste0(feature, "_A", ".png"), plot_responseA(feature), device = 'png',
                   path = '../img/',
                   scale = 1.5, width = 3, height = 3, units = "in")
}

ggsave(paste0("chol_HDL_A", ".png"), plot_responseA(features[6]), device = 'png',
       path = '../img/',
       scale = 1.5, width = 3, height = 3, units = "in")

plot_responseB <- function(x){
            # x is one of the features from glucose/insulin/lipid panels
            df = data.frame(
                        value = trtB$conc_table[x,],
                        subject = trtB$sample_table[, 'subject'],
                        timepoint = trtB$sample_table[, 'timepoint'],
                        response = factor(
                                    rep(responseB[x,], each = 2),
                                    labels = c("non-responder", "responder")
                        )
            )
            unit = lipids$feature_data[x ,1]
            p <- ggplot(df, aes(x = timepoint, y = value)) +
                        stat_boxplot(geom = "errorbar", width = 0.5)+
                        geom_boxplot(size = 0.5) +
                        geom_point(aes(color = subject)) +
                        geom_line(aes(group = subject, color = subject), size = 0.5) +
                        facet_wrap(~ response) +
                        labs(
                                    title = paste("Treatment B", x, sep = " "),
                                    y = unit
                        ) +
                        theme_boxplot() +
                        scale_color_manual(values = colorRampPalette(pal_lancet()(9))(20))
            theme(
                        legend.text = element_text(size = 7),
                        legend.key.height = unit(0.5, "cm"),
                        legend.key.width = unit(0.5, "cm")
            )
            return(p)
}

features <- rownames(trtB$feature_data)
for (feature in features) {
            ggsave(paste0(feature, "_B", ".png"), plot_responseB(feature), device = 'png',
                   path = '../img/',
                   scale = 1.5, width = 3, height = 3, units = "in")
}

ggsave(paste0("chol_HDL_B", ".png"), plot_responseB(features[6]), device = 'png',
       path = '../img/',
       scale = 1.5, width = 3, height = 3, units = "in")

write.csv(responseA, "../img/responseA.csv")
write.csv(responseB, "../img/responseB.csv")
