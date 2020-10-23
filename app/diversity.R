taxa <- readRDS("data/taxa.rds")

toPhyloseq <- function(taxa, taxaLevel){
    dataset <- HTSet::summarize_feature(taxa, taxaLevel)
    otu <- phyloseq::otu_table(dataset$edata, taxa_are_rows = TRUE)
    sam <- dataset$pdata %>%
        tibble::rownames_to_column() %>%
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
        tibble::column_to_rownames() %>%
        phyloseq::sample_data()
    tax <- rownames(dataset$edata)
    names(tax) <- rownames(dataset$edata)
    tax <- as.matrix(tax) %>% phyloseq::tax_table()
    phyloseq::phyloseq(otu, tax, sam)
}
 
anovaTable <- function(richness, taxa2, contrast = c("pre", "post")){
    df <- cbind(richness, phyloseq::sample_data(taxa2)) %>%
        filter(trt %in% c('A', 'B')) %>%
        filter(timepoint %in% contrast) %>%
        mutate(
            trt = droplevels(trt),
            timepoint = droplevels(timepoint)
        )
    model <- lm(Shannon~trt+timepoint+trt:timepoint+subject, data = df)
    anova(model)
}

timecourse <- function(richness, taxa2, adjust = T){
    df <- cbind(richness, phyloseq::sample_data(taxa2))
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
dboxplot <- function(richness, taxa2, tp = c("pre", "post")){
    df <- cbind(richness, phyloseq::sample_data(taxa2)) %>%
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

betaDiversity <- function(taxa2, dist_methods, ordinate_method){
    taxa_relative <- phyloseq::transform_sample_counts(taxa2, function(x) x / sum(x) )
    taxa_relative <- phyloseq::filter_taxa(taxa_relative, function(x) mean(x) > 1e-5, TRUE)
    # Calculate distance matrix
    iDist <- phyloseq::distance(taxa_relative, method=dist_methods)
    # Calculate ordination
    ord  <- phyloseq::ordinate(taxa_relative, ordinate_method, distance=iDist)
}
# plot_scree(betaDiversity(), "Scree plot for species, Bray-Curtis/PCoA")

plotBetaDiversity <- function(taxa2, ord, color_by, timecode){
    taxa_relative <- phyloseq::transform_sample_counts(taxa2, function(x) x / sum(x) )
    taxa_relative <- phyloseq::filter_taxa(taxa_relative, function(x) mean(x) > 1e-5, TRUE)
    logic <- phyloseq::sample_data(taxa_relative)$timepoint_code %in% timecode
    taxa_relative2 <- phyloseq::prune_samples(logic, taxa_relative)
    phyloseq::plot_ordination(taxa_relative2, ord, "samples", color=color_by) +
        geom_point(size = 2) + 
        # geom_line(aes(group = subject)) +
        # stat_ellipse() + 
        theme_bw()
    # ggplot(df, aes_string(colnames(df)[1], colnames(df)[2], color = color_by))+
    #     geom_point()+
    #     theme_bw()
} 

diversityUI <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            box(
                title = "ANOVA Test for Alpha Diversity",
                width = 6,
                solidHeader = T,
                status = "primary",
                selectInput(
                    ns("taxaLevel"), "Taxonomic Level", 
                    choices = c("phylum", "class", "order", "family","genus", "species", "strain"),
                    selected = "species"
                ),
                radioButtons(
                    ns('contrast'),
                    'Select the timepoint to contrast',
                    choices = c('post-pre', 'mid-pre'),
                    selected = 'post-pre',
                    inline = T
                ),
                dataTableOutput(ns("anova"))
            ),
            tabBox(
                title = "Alpha Diversity Visulization",
                width = 6,
                tabPanel(
                    "Time Series", 
                    plotlyOutput(ns("timeplot"))
                ),
                tabPanel(
                    "Boxplot", 
                    inputPanel(
                        radioButtons(
                            ns('tp1'), 'Timepoint',
                            choices = c('pre-post', 'pre-mid-post'),
                            selected = 'pre-post',
                            inline = T
                        )
                    ),
                    plotlyOutput(ns("boxplot"))
                )
            )
        ),
        fluidRow(
            box(
                title = "Control Panel for Beta Diversity Visulization",
                width = 4,
                solidHeader = T,
                status = "primary",
                inputPanel(
                    selectInput(
                        ns("distance"), "Distance Method", 
                        choices = phyloseq::distanceMethodList$vegdist,
                        selected = "bray"
                    ),
                    selectInput(
                        ns("ordinate"), "Ordinate Method", 
                        choices = c("PCoA", "MDS"),
                        selected = "PCoA"
                    ),
                    selectInput(
                        ns("legend"), "The color of legends represents:", 
                        choices = c("annotation" = "annotation", 
                                    "timepoint" = "timepoint", 
                                    "trt" = "trt", 
                                    "subject" = "subject",
                                    "gender" = "gender", 
                                    "age" = "age2", 
                                    "BMI" = "bmi"),
                        selected = "annotation"
                    ),
                    checkboxGroupInput(
                        ns("tp2"), "Select the timepoints presented in the distance plot:", 
                        choices = LETTERS[1:7],
                        selected = LETTERS[1:7],
                        inline = T
                    ),
                    actionButton(
                        ns('ok'),
                        'OK'
                    )
                )
            ),
            box(
                title = "Beta Diversity Visulization",
                width = 8,
                status = "primary",
                h3('Click "OK" to generate a new plot'),
                plotlyOutput(ns('beta'))
            )
        )
    )
}
diversityServer <- function(id) {
    moduleServer(
        id,
        function(input, output, session) {
            taxa2 <- reactive({toPhyloseq(taxa, input$taxaLevel)})
            richness <- reactive({
                phyloseq::estimate_richness(taxa2(), measures = c("Shannon"))
            })
            output$anova <- renderDataTable({
                contrast <- str_split(input$contrast, "-")[[1]]
                anovaTable(richness(), taxa2(), contrast) %>%
                    datatable(
                        rownames = TRUE,
                        selection = list(mode = 'single', selected = 1)
                    ) %>%
                    formatSignif(1:5, digits = 3)
            })
            output$timeplot <- renderPlotly({
                timecourse(richness(), taxa2())
            })
            output$boxplot <- renderPlotly({
                tp1 <- str_split(input$tp1, "-")[[1]]
                dboxplot(richness(), taxa2(), tp1)
            })
            beta <- eventReactive(input$ok, {
                ord <- betaDiversity(taxa2(), input$distance, input$ordinate)
                plotBetaDiversity(taxa2(), ord, input$legend, input$tp2)
            })
            output$beta <- renderPlotly({
                beta()
            })
        }
    )
}

# server = function(input, output, session) {
#     diversityServer("a")
# }
# 
# shinyApp(ui = diversityUI("a"), server = server)