library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(shinydashboard)
library(dplyr)
library(tibble)
library(ggplot2)
library(plotly)
library(ggsci)
library(factoextra)
library(pheatmap)
library(HTSet)
library(DT)
library(ggfortify)
# library(shinyauthr)
# library(shinyjs)
# library(shiny.router)
library(enrichplot)
library(stringr)
library(edgeR)
library(DESeq2)

dashboardHeader = eval(parse(text = deparse(shinydashboard::dashboardHeader)[-5]))
