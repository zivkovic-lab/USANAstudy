scfa <- readRDS("data/scfa.rds")

fitScfa <- function(contrast = c("pre", "post")) {
    scfa = subset_samples(scfa, scfa$pdata$trt %in% c('A', 'B'))
    scfa$pdata$trt = droplevels(scfa$pdata$trt)
    scfa = subset_samples(scfa, scfa$pdata$timepoint %in% contrast)
    scfa$pdata$timepoint = droplevels(scfa$pdata$timepoint)
    design = model.matrix(~ trt * timepoint + subject, data = scfa$pdata)
    fit = lmFit(log(scfa$edata+0.1), design = design) %>%
        eBayes() %>%
        topTable(coef = colnames(design)[ncol(design)], number = Inf, sort.by = 'none')
    return(fit)
}

sboxplot <- function(feature, timepoint = c("pre", "post")){
    # FIX ME! input of feature and timepoint
    dataset = subset_samples(scfa, scfa$pdata$trt %in% c('A', 'B'))
    dataset$pdata$trt = droplevels(dataset$pdata$trt)
    dataset <- subset_samples(dataset, dataset$pdata$timepoint %in% timepoint)
    dataset$pdata$timepoint <- droplevels(dataset$pdata$timepoint)
    df = dataset$edata %>%
        t() %>%
        cbind(dataset$pdata[c("trt", "timepoint", "subject")])
    ggplot(df, aes_string("timepoint", feature))+
        geom_boxplot()+
        geom_point(aes(color = subject))+
        geom_line(aes(group = subject, color = subject))+
        facet_wrap(~trt)+
        theme_bw()+
        labs(title = feature)
}

scfadeUI <- function(id) {
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
                    ns("contrast"), "Select the two timepoints to compare", 
                    choices = c("post-pre", 
                                "mid-pre"),
                    selected = "post-pre",
                    inline = T
                ),
                dataTableOutput(ns("de"))
            ),
            box(
                title = "Histogram",
                width = NULL,
                solidHeader = T,
                status = "primary",
                radioButtons(
                    ns("transform"), "Data Transform", 
                    choices = c("Original", "log(x+0.1)"),
                    selected = "Original",
                    inline = T
                ),
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
                radioButtons(
                    ns("tp"), "Timepoint", 
                    choices = c("pre-post", "pre-mid-post"),
                    selected = "pre-post",
                    inline = T
                ),
                plotlyOutput(ns("plot"))
            )
        )
    )
}

scfadeServer <- function(id) {
    moduleServer(
        id,
        function(input, output, session) {
            res <- reactive({
                contrast <- str_split(input$contrast, "-")[[1]]
                fitScfa(contrast)
            })
            output$de <- renderDataTable({
                datatable(res(), 
                          rownames = TRUE,
                          selection = list(mode = 'single', selected = 1)
                ) %>%
                    formatSignif(columns = 1:6, digits = 3)
            })
            output$plot <- renderPlotly({
                feature <- rownames(res())[input$de_rows_selected]
                tp <- str_split(input$tp, "-")[[1]]
                sboxplot(feature, tp)
            })
            output$hist <- renderPlotly({
                if (input$transform == "Original") {
                    df <- t(scfa$edata)
                } else {
                    df <- log(scfa$edata + 1) %>%
                        t()
                }
                feature <- rownames(res())[input$de_rows_selected]
                ggplot(df) +
                    geom_histogram(aes_string(feature), bins = 40) +
                    labs(title = feature) +
                    theme_bw()
            })
        }
    )
}
