source("global.R")
source("taxa.R")
source("functions.R")
source("functions-pathway.R")
source("functions-module.R")
source("scfa.R")
source('diversity.R')
source("dashboard.R")
source("authr.R")
source("user_base.R")

router <- make_router(
    default = route(
        "/login",
        ui = authrUI("NULL"),
        server = function(input, output, session) {
            callModule(authrServer, "NULL")
        }
    ),
    route(
        "dashboard",
        ui = dashboardUI("dashboard"),
        server = function(input, output, session) {
            dashboardServer("dashboard")
        }
    )
)

ui <- fluidPage(
    title = "USANA",
    router$ui
)

# router <- make_router(
#     route(
#         "/",
#         ui = authrUI("auth"),
#         server = function(input, output, session) {
#             callModule(authrServer, "auth")
#         }
#     ),
#     route(
#         "dashboard",
#         ui = dashboardUI("dashboard"),
#         server = function(input, output, session) {
#             dashboardServer("dashboard")
#         }
#     )
# )
# 
# ui <- fluidPage(
#     title = "Router demo",
#     router$ui
# )

# server <- function(input, output, session) {
#     router$server(input, output, session)
# }

server <- shinyServer(function(input, output, session) {
    credentials = callModule(authrServer, "NULL")
    router$server(input, output, session)
    observe({
        if(is.null(credentials())){
            change_page("/login", mode = "push")
        }else if (is.null(credentials()$user_auth)) {
            change_page("/login", mode = "push")
        }else if(credentials()$user_auth){
            change_page("/dashboard", mode = "push")
        } else {
            print(credentials()$user_auth)
            change_page("/login", mode = "push")
        }
    })
    observe({
        if(is_page("dashboard")) {
            session$onFlushed(function(){
                shinyjs::addClass(selector = "html body", class = "skin-blue")
            }, once = FALSE)
        }
    })
})

shinyApp(ui, server)