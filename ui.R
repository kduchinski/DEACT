library(shiny)
library(mlbench)
library(plotly)
library(dplyr)


# ui.R definition
ui <- fluidPage(
  navbarPage("Dataset Exploration Tool",
             tabPanel("Genes Affected in Both Conditions", {
               fluidRow(
                 column(3, 
                   fileInput("file", "Upload data"),
                   p("Use the lasso or box select tools to choose points of interest."),
                   p("The selected data will be previewed below. Press the button to download."),
                   downloadButton('downloadData', 'Download selected data'),
                   dataTableOutput("StatsTable")
                   ),
                 column(1,
                   selectInput("ChangeChoice", NULL, list("FPKM Difference", "log2 Fold Change"), selected = "fpkm diff", multiple = FALSE),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   actionButton("quad1Button", "Add to dataset"),
                   actionButton("quad3Button", "Add to dataset"),
                   actionButton("quad4Button", "Add to dataset"),
                   actionButton("quad2Button", "Add to dataset")),
                 column(7, (plotlyOutput("Plot", height = "600px")), offset = 1)
                 ) 
             },
             dataTableOutput("DataTable"),
             checkboxGroupInput("infChoice", "Include the following infinite fold change values:",
                                c("+Inf for Loss" = 'posInfLoss',
                                  "-Inf for Loss" = 'negInfLoss',
                                  "+Inf for Gain" = 'posInfGain',
                                  "-Inf for Gain" = 'negInfGain'),
                                inline = T),
             radioButtons("showInf", NULL, c("Show all infinite points", "Show selected infinite points"), inline = T),
             dataTableOutput("InfTableBoth")
             ), tabPanel("Full Dataset", {
               sidebarLayout(
                 sidebarPanel(
                   p("Use the lasso or box select tools to choose points of interest."),
                   p("The selected data will be previewed below."),
                   p("Press the button below to add these points to the full selection on the other page."),
                   actionButton('addButton', "Add selection to dataset")
                 ),
                 mainPanel(
                   plotlyOutput("Plot1", height = "600px")
                 ))
             },
             dataTableOutput("extraPoints"),
             checkboxInput("posInfLoss", label = "+Inf for Loss", value = F),
             checkboxInput("negInfLoss", label = "-Inf for Loss", value = F),
             checkboxInput("posInfGain", label = "+Inf for Gain", value = F),
             checkboxInput("negInfGain", label = "-Inf for Gain", value = F),
             radioButtons("showInf", NULL, c("Show all infinite points" = "all", "Show selected infinite points" = "some"), selected = "all", inline = T),
             dataTableOutput("InfTable")
            ), windowTitle = "DET"
  )
)