metabolites <- readRDS("data/metabolites.rds")
metabolites <- metabolites[apply(metabolites$edata, 1, var) != 0,]

calDE <- function(transform = "log"){
    if (transform == "log") {
        edata <- log(metabolites$edata)
    } else {
        edata <- metabolites$edata
    }
    model <- model.matrix(~tp*trt+subject_id, data = metabolites$pdata)
    fit <- lmFit(edata, model) %>% eBayes()
    topTable(fit, coef = "tppost:trtB", number = Inf, sort.by = "none")
}
plotVolc <- function(topTable){
    topTable %>%
        rownames_to_column("Name") %>%
    ggplot(aes(logFC, -log(P.Value))) +
        geom_point(aes(Name = Name), shape = 21)+
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
    if (identical(rownames(topTable), HTSet::featureNames(metabolites))) {
        # print(rownames(metabolites$edata)[rows_selected])
        dt <- data.frame(
            Response = metabolites$edata[rows_selected,],
            Timepoint = metabolites$pdata$tp,
            Treatment = metabolites$pdata$trt,
            Subject_id = metabolites$pdata$subject_id
        ) 
        # print(dt)
        ggplot(dt, aes(x = Timepoint, y = Response))+
            geom_boxplot()+
            geom_point(aes(group = Treatment, color = Subject_id), shape = 21) +
            geom_line(aes(group = Subject_id, color = Subject_id)) +
            facet_wrap(~Treatment)+
            labs(title = rownames(topTable)[rows_selected]) +
            theme_bw()
    } else {
        print("The rownames of DE table and original edata table are not identical.")
    }
}

metabDEUI <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            6,
            box(
                title = "Differential Expression",
                width = NULL,
                solidHeader = T,
                status = "primary",
                radioButtons(
                    ns("transform"), "Select data transformation type", 
                    choices = c("log", 
                                "none"),
                    selected = "log",
                    inline = T
                ),
                dataTableOutput(ns("de"))
            ),
            box(
                title = "P-Value Distribution",
                width = NULL,
                solidHeader = T,
                status = "primary",
                plotlyOutput(ns("hist"))
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
                plotlyOutput(ns("box"))
            ),
            box(
                title = "Volcano Plot",
                width = NULL,
                solidHeader = T,
                status = "primary",
                plotlyOutput(ns("volc"))
            )
        )
    )
}
metabDEServer <- function(id) {
    moduleServer(
        id,
        function(input, output, session) {
            DEtbl <- reactive({
                calDE(input$transform)
            })
            output$de <- renderDataTable({
                datatable(DEtbl(), 
                          rownames = TRUE,
                          selection = list(mode = 'single', selected = 1)
                ) %>%
                    formatSignif(columns = 1:6, digits = 3)
            })
            output$box <- renderPlotly({
                plotBox(DEtbl(), input$de_rows_selected)
            })
            output$hist <- renderPlotly({
                plotHistP(DEtbl())
            })
            output$volc <- renderPlotly({
                plotVolc(DEtbl())
            })
        }
    )
}
# server = function(input, output, session) {
#     metabDEServer("a")
# }
# 
# shinyApp(ui = metabDEUI("a"), server = server)