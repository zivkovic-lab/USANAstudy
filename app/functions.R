load("data/functions.rda")
functions <- readRDS('data/functions.rds')

fboxplot <- function(dataset, feature, tp){
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

fviolinplot <- function(dataset, feature, timepoint = c("pre", "post")){
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
funcdeUI <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            width = 6,
            box(
                title = "Differential Expression",
                width = NULL,
                solidHeader = T,
                status = "primary",
                radioButtons(
                    ns('dataset'),
                    "Dataset",
                    c('KEGG Orthology (includes both enzyme and non-enzyme genes' = 'ko',
                      'Enzyme' = 'enzyme'),
                    'ko',
                    inline = T
                ),
                dataTableOutput(ns("de"))
            )
        ),
        column(
            width = 6,
            box(
                title = "Plots",
                width = NULL,
                solidHeader = T,
                status = "primary",
                inputPanel(
                    radioButtons(
                        ns('tp'),
                        'Timepoint',
                        c('pre-post', 'pre-mid-post'),
                        'pre-post',
                        inline = T
                    ),
                    radioButtons(
                        ns('type'),
                        'Plot Type',
                        c('Boxplot', 'Violin Plot'),
                        'Boxplot',
                        inline = T
                    )
                ),
                plotlyOutput(ns("plot"))
            )
        )
    )
}
funcdeServer <- function(id){
    moduleServer(
        id,
        function(input, output, sesion){
            output$de <- renderDataTable({
                if (input$dataset == 'ko'){
                    annotation <- functions$ko$fdata[rownames(koRes$results),]
                    koRes$results %>%
                        cbind(annotation = annotation) %>%
                        datatable(
                            selection = list(mode = 'single', selected = 1)
                        ) %>%
                        formatSignif(columns = 1:5)
                } else {
                    datatable(
                        enzymeRes$results, 
                        selection = list(mode = 'single', selected = 1)
                    ) %>%
                    formatSignif(columns = 1:5)
                }
                
            })
            output$plot <- renderPlotly({
                tp <- str_split(input$tp, "-")[[1]]
                if (input$dataset == 'ko'){
                    dataset = functions$ko
                    feature <- rownames(koRes$results)[input$de_rows_selected]
                } else {
                    dataset = functions$enzyme
                    feature <- rownames(enzymeRes$results)[input$de_rows_selected]
                }
                if (input$type == 'Boxplot') {
                    fboxplot(dataset, feature, tp)
                } else {
                    fviolinplot(dataset, feature, tp)
                }
            })
        }
    )
}

# server = function(input, output, session) {
#     funcdeServer("a")
# }
# 
# shinyApp(ui = funcdeUI("a"), server = server)