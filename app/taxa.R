taxa <- readRDS("data/taxa.rds")

fitTaxa <- function(taxaLevel) {
    # taxaLevel can be "kingdom", "phylum", "class", "order", 
    # "family","genus", "species", "strain"
    genus = summarize_feature(taxa, taxaLevel)
    genus = subset_samples(genus, genus$pdata$timepoint %in% c("pre", "post"))
    genus$pdata$timepoint = droplevels(genus$pdata$timepoint)
    genus$pdata$trt = factor(genus$pdata$trt)
    genus = genus[rowSums(genus$edata)>40,]
    genus = genus[apply(genus$edata, 1, function(x){sum(x != 0)})>4,]
    # genus$pdata$trt = factor(genus$pdata$trt, levels = c('B', 'A'))
    design = model.matrix(~ trt * timepoint + subject, data = genus$pdata)
    fit = model_fit(genus, design = design, engine = "DESeq2", coef = "trtB:timepointpost")
    return(fit)
}


tboxplot <- function(feature, taxaLevel, timepoint = c("pre", "post")){
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

tviolinplot <- function(feature, taxaLevel, timepoint = c("pre", "post")) {
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

deUI <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            6,
            box(
                title = "Differential Expression",
                width = NULL,
                solidHeader = T,
                status = "primary",
                selectInput(
                    ns("taxaLevel"), "Taxonomic Level", 
                    choices = c("phylum", "class", "order", "family","genus", "species", "strain"),
                    selected = "genus"
                ),
                dataTableOutput(ns("de"))
            )
        ),
        column(
            6,
            box(
                title = "Boxplot",
                width = NULL,
                solidHeader = T,
                status = "primary",
                collapsible = T,
                inputPanel(
                    radioButtons(
                        ns('tp'),
                        'Timepoint',
                        choices = c("pre-post", "pre-mid-post"),
                        selected = 'pre-post',
                        inline = T
                    ),
                    radioButtons(
                        ns('type'),
                        'Plot Type',
                        choices = c("Boxplot", "Violin plot"),
                        selected = 'Boxplot',
                        inline = T
                    )
                ),
                plotlyOutput(ns("plot")),
            )
        )
    )
}

deServer <- function(id) {
    moduleServer(
        id,
        function(input, output, session) {
            res <- reactive(fitTaxa(input$taxaLevel))
            output$de <- renderDataTable({
                fitRes <- res()$results
                datatable(fitRes, 
                          rownames = TRUE,
                          selection = list(mode = 'single', selected = 1)
                ) %>%
                    formatSignif(columns = 1:5, digits = 3)
            })
            output$plot <- renderPlotly({
                fitRes <- res()$results
                # print(input$de_rows_selected)
                feature <- rownames(fitRes)[input$de_rows_selected]
                tp <- str_split(input$tp, "-")[[1]]
                if (input$type == 'Boxplot'){
                    tboxplot(feature, input$taxaLevel, tp)
                } else {
                    tviolinplot(feature, input$taxaLevel, tp)
                }
            })
        }
    )
}

# server = function(input, output, session) {
#     deServer("a")
# }
# 
# shinyApp(ui = deUI("a"), server = server)