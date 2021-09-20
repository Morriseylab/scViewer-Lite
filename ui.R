library(shiny)
library(shinyBS)
library(plotly)
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)
library(shinythemes)

options(shiny.sanitize.errors = FALSE)
options(shiny.maxRequestSize=600*1024^2) 



ui <- fluidPage(theme = shinytheme("sandstone"),
                tags$style(HTML("
                        #primariy2 .box.box-solid.box-primary>.box-header {
                        color:#fff;
                        background:#666666
                        }

                        .box.box-solid.box-primary {
                        border-bottom-color:#666666;
                        border-left-color:#666666;
                        border-right-color:#666666;
                        border-top-color:#666666;
                        }
                        
                        .shiny-notification {
                          height: 100px;
                          width: 800px;
                          position:fixed;
                          top: calc(50% - 50px);;
                          left: calc(50% - 400px);;
                        }
                        
                        ")),
                navbarPage("Morrisey Lab Lung Niche Exploration",
                   #shinythemes::themeSelector(), 
                   ###### Here : insert shinydashboard dependencies ######
                   header = tagList(
                     useShinydashboard()
                   ),
######################################################################################################################################
                   
                   tabPanel("About",
                            fluidRow(column(4,h1()),
                            column(6,box(solidHeader = FALSE,width=9,align="left",
                           includeHTML("include.html"))),
                           column(3,h1()))
                           
                   ),
######################################################################################################################################
tabPanel("Explore Data",
         sidebarLayout(
           sidebarPanel(h3("Controls",align="center"),br(),width = 4,
                        uiOutput("projectlist"),       
                            fluidRow(
                              column(6,selectInput("categorya2", "Select plot A variable",c('Celltype' = "clust", 'Gene Expression' = "geneexp"),selected = "clust")),
                              column(6,selectInput("categoryb2", "Select plot B variable",c('Celltype' = "clust", 'Gene Expression' = "geneexp"),selected = "clust"))
                            ),
                            fluidRow(
                              column(6,conditionalPanel(
                                condition = "input.categorya2 == 'geneexp'",uiOutput("gene1aui")
                                   )),
                              column(6,conditionalPanel(
                                condition = "input.categoryb2 == 'geneexp'",uiOutput("gene2aui")
                              ))),
                        fluidRow(
                         column(6,conditionalPanel(
                            condition = "input.categorya2 == 'clust' || input.categorya2 == 'var'",
                            checkboxInput("checklabel1", label = "Check for cell group labelling (A)", value = TRUE)
                            )),
                          column(6,conditionalPanel(
                            condition = "input.categoryb2 == 'clust' || input.categoryb2 == 'var'",
                            checkboxInput("checklabel2", label = "Check for cell group labelling (B)", value = TRUE)
                          ))),
                        sliderInput("pointa2", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                        downloadButton('downloadtsneplot', 'Download Plot')
              ),
           mainPanel(
             tags$div(class = "another-box", id = "primariy2",
             box(title = textOutput('title'),solidHeader = TRUE,width=12,status='primary',
                 plotOutput("comptsne2", height = 1000)),
              ))
         )
)

)
)