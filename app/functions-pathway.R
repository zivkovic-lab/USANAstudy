load("data/functions.rda")

pathwayUI <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            box(
                title = "KEGG Pathway Enrichment Analysis (gsea)",
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
                    nrow(gseP@result),
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

pathwayServer <- function(id) {
    moduleServer(
        id,
        function(input, output, session) {
            output$gse <- renderDataTable({
                datatable(gseP@result, 
                          rownames = FALSE,
                          selection = list(mode = 'single', selected = 1)
                ) %>%
                    formatSignif(columns = 4:8, digits = 3)
            })
            output$dot_gse <- renderPlot({
                dotplot(gseP, showCategory = input$categories)
            })
            output$gseaPlot <- renderPlot({
                gseaplot2(gseP, geneSetID = input$gse_rows_selected, title = gseP$Description[input$gse_rows_selected])
            })
        }
    )
}
