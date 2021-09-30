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
                navbarPage("LungMap Reference Data",
                   #shinythemes::themeSelector(), 
                   ###### Here : insert shinydashboard dependencies ######
                   header = tagList(
                     useShinydashboard()
                   ),

######################################################################################################################################
tabPanel("Explore Data",
         sidebarLayout(
           sidebarPanel(h3("Controls",align="center"),br(),width = 2,
                        box(title="Plot A",width=12,
                        uiOutput("umapa"), #Dimensionality reduction method of left plot
                        selectInput("categorya2", "Select plot A variable",c('Celltype' = "clust", 'Gene Expression' = "geneexp"),selected = "clust"),
                        conditionalPanel(condition = "input.categorya2 == 'geneexp'",uiOutput("gene1aui")),
                        conditionalPanel(condition = "input.categorya2 == 'clust'",
                         checkboxInput("checklabel1", label = "Check for cell group labelling (A)", value = FALSE)),
                        checkboxInput("subsa", label = "Check to highlight cells (A)", value = FALSE),
                        conditionalPanel(condition = "input.subsa ==true",uiOutput("subsaui")) #generate ident list for left plot
                        ),br(),
                        box(title="Plot B",width=12,
                        uiOutput("umapb"), #Dimensionality reduction method of left plot
                        selectInput("categoryb2", "Select plot B variable",c('Celltype' = "clust", 'Gene Expression' = "geneexp"),selected = "clust"),
                        conditionalPanel(condition = "input.categoryb2 == 'geneexp'",uiOutput("gene2aui")),
                        conditionalPanel(condition = "input.categoryb2 == 'clust'",
                                         checkboxInput("checklabel2", label = "Check for cell group labelling (B)", value = FALSE)),
                        checkboxInput("subsb", label = "Check to highlight cells (B)", value = FALSE),
                        conditionalPanel(condition = "input.subsb ==true",uiOutput("subsbui")) #generate ident list for right plot
                        ),br(),
                        sliderInput("pointa2", "Point Size:",min = 0, max = 5, value = 1,step=.25),
                        downloadButton('downloadtsneplot', 'Download Plot')
              ),
           mainPanel(width=10,
             tags$div(class = "another-box", id = "primariy2",
             box(title = textOutput('title'),solidHeader = TRUE,width=12,status='primary',
                 plotOutput("comptsne2", height = 1000)),
              ))
         )
)

)
)