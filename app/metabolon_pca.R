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

plotPCA <- function(HTSetObj = dat4PCA, color_by = "group_name"){
    df <- base::t(HTSetObj$edata)
    res.pca <- prcomp(df, center = T, scale. = T)
    fviz_pca_ind(res.pca,
                 label = "none", # hide individual labels
                 habillage = HTSetObj$pdata[[color_by]], # color by groups
                 addEllipses = TRUE # Concentration ellipses
    )
}
plotHeatmap <- function(HTSetObj = dat4PCA, ann_col = "group_name", ann_row = "super_pathway"){
    mat <- scale(t(HTSetObj$edata)) %>% t()
    mat[mat > 2] <- 2
    mat[mat < -2] <- -2
    # Generate annotations for rows and columns
    annotation_col <- data.frame(
        group = HTSetObj$pdata[,ann_col]
    )
    rownames(annotation_col) <- HTSet::sampleNames(HTSetObj)
    annotation_row <- data.frame(
        pathway = HTSetObj$fdata[,ann_row]
    )
    rownames(annotation_row) <- HTSet::featureNames(HTSetObj)
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
metabPCAUI <- function(id) {
    ns <- NS(id)
    tagList(
        h1("Note: Chemical compounds with P-Value < 0.05 are used for PCA and heatmap."),
        column(
            4,
            box(
                title = "PCA",
                width = NULL,
                solidHeader = T,
                status = "primary",
                selectInput(
                    ns("color"), "Color by", 
                    choices = colnames(metabolites$pdata)[c(10,11,19,24,25)],
                    selected = "group_name",
                ),
                plotlyOutput(ns("pca"))
            )
        ),
        column(
            8,
            box(
                title = "Heatmap",
                width = NULL,
                solidHeader = T,
                status = "primary",
                selectInput(
                    ns("ann_col"), "Column Annotation", 
                    choices = colnames(metabolites$pdata)[c(10,11,19,24,25)],
                    selected = "group_name",
                ),
                selectInput(
                    ns("ann_row"), "Row Annotation", 
                    choices = c("super_pathway", "sub_pathway"),
                    selected = "super_pathway",
                ),
                plotOutput(ns("heat"))
            )
        )
    )
}
metabPCAServer <- function(id) {
    moduleServer(
        id,
        function(input, output, session) {
            DEgenes <- calDE("log") %>%
                dplyr::filter(P.Value<0.05) %>%
                rownames()
            dat4PCA <- metabolites[DEgenes,]
            output$pca <- renderPlotly({
                plotPCA(HTSetObj = dat4PCA, input$color)
            })
            output$heat <- renderPlot({
                print(input$ann_col)
                plotHeatmap(HTSetObj = dat4PCA, ann_col = input$ann_col, ann_row = input$ann_row)
            })
        }
    )
}
# server = function(input, output, session) {
#     metabPCAServer("a")
# }
# 
# shinyApp(ui = metabPCAUI("a"), server = server)