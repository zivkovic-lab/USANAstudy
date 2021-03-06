# source("global.R")
# source("taxa.R")
# source("functions.R")
# source("functions-pathway.R")
# source("functions-module.R")
# source("scfa.R")
# source('diversity.R')

dashboardUI <- function(id) {
        ns <- NS(id)
        header <- dashboardHeader(
                title = "USANA Fiber Study",
                tags$li(
                        class = "nav-item",
                        tags$a(
                                class = "btn btn-danger action-button",
                                id = "logout-button",
                                type = "button",
                                icon("sign-out-alt"), "Log out"
                        )
                )
        )
        
        sidebar <- dashboardSidebar(
                sidebarMenu(
                        # first menu
                        menuItem("Taxonomic Analysis", tabName = "taxa",
                                 menuSubItem("Alpha/beta Diversity", tabName = 'taxa-diversity'),
                                 menuSubItem("Differential Expression", tabName = "taxa-de")
                                 
                        ),
                        # second menu
                        menuItem("Functional Analysis", tabName = "func",
                                 menuSubItem("Differential Expression", tabName = "func-de"),
                                 menuSubItem("KEGG Pathway", tabName = "func-pathway"),
                                 menuSubItem("KEGG Module", tabName = "func-module")
                        ),
                        menuItem("Metabolites", tabName = "metabolites",
                                  menuSubItem("Short Chain Fatty Acids", tabName = "scfa")
                        )
                )
        )
        
        body <- dashboardBody(
                shinyjs::useShinyjs(),
                tags$link(href="styles.css", rel="stylesheet"),
                tabItems(
                        # the first subtab content
                        tabItem(
                                tabName = "taxa-diversity",
                                diversityUI(ns("taxa-diversity"))
                        ),
                        tabItem(
                                tabName = "taxa-de",
                                deUI(ns("taxa-de"))
                        ),
                        # the second subtab content
                        tabItem(
                                tabName = "func-de",
                                funcdeUI(ns("func-de"))
                        ),
                        # the third subtab content
                        tabItem(
                                tabName = "func-pathway",
                                pathwayUI(ns("func-pathway"))
                        ),
                        tabItem(
                                tabName = "func-module",
                                funcmoduleUI(ns("func-module"))
                        ),
                        tabItem(
                                tabName = "scfa",
                                scfadeUI(ns("scfa"))
                        )
                )
        )
        dashboardPage(header, sidebar, body)
}

dashboardServer <- function(id) {
        moduleServer(
                id,
                function(input, output, session) {
                        deServer("taxa-de")
                        diversityServer("taxa-diversity")
                        funcdeServer("func-de")
                        pathwayServer("func-pathway")
                        funcmoduleServer("func-module")
                        scfadeServer("scfa")
                }
        )
}
# server = function(input, output, session) {
#         dashboardServer("a")
# }
# 
# shinyApp(ui = dashboardUI("a"), server = server)