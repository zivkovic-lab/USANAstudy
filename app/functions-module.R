load("data/functions.rda")

funcmoduleUI <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            box(
                title = "KEGG Module Enrichment Analysis (gsea)",
                width = 6,
                solidHeader = T,
                status = "primary",
                dataTableOutput(ns("gse"))
            ),
            box(
                title = "Dot Plot for GSEA",
                width = 6,
                solidHeader = T,
                status = "primary",
                sliderInput(
                    ns("categories"), 
                    "The number of categories", 
                    5, 
                    nrow(gseM@result),
                    10, step = 1
                ),
                plotOutput(ns("dot_gse"))
            )
        ),
        fluidRow(
            box(
                title = "Visulization of GSEA Results",
                width = 12,
                solidHeader = T,
                status = "primary",
                h3(textOutput(ns("title")), align='center'),
                plotOutput(ns("gseaPlot"))
            )
        )
    )
}

funcmoduleServer <- function(id) {
    moduleServer(
        id,
        function(input, output, session) {
            output$gse <- renderDataTable({
                datatable(gseM@result, 
                          rownames = FALSE,
                          selection = list(mode = 'single', selected = 1)
                ) %>%
                    formatSignif(columns = 4:8, digits = 3)
            })
            output$dot_gse <- renderPlot({
                dotplot(gseM, showCategory = input$categories)
            })
            output$gseaPlot <- renderPlot({
                gseaplot2(gseM, geneSetID = input$gse_rows_selected, title = gseM$Description[input$gse_rows_selected])
            })
        }
    )
}
