library(shiny)
library(mlbench)
library(plotly)
library(dplyr)


# ui.R definition
ui <- fluidPage(
  navbarPage("Dataset Exploration And Curation Tool",
             tabPanel("Genes Affected in Both Conditions", {
               fluidRow(
                 column(3, 
                        fileInput("file", "Upload data"),
                        p("Use the lasso or box select tools to choose points of interest."),
                        p("The selected data will be previewed below. Press the button to download."),
                        downloadButton('downloadData', 'Download selected data'),
                        dataTableOutput("StatsTable")
                 ),
                 column(2,
                        selectInput("ChangeChoice", NULL, list("FPKM Difference", "log2 Fold Change"), selected = "FPKM Difference", multiple = FALSE, width = "150px"),
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
                        p(""),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        actionButton("quad1Button", "Add to dataset", width = "150px"),
                        actionButton("quad3Button", "Add to dataset", width = "150px"),
                        actionButton("quad4Button", "Add to dataset", width = "150px"),
                        actionButton("quad2Button", "Add to dataset", width = "150px")),
                 column(6, (plotlyOutput("Plot", height = "600px")), offset = 0.5)
               ) 
             },
             conditionalPanel(
               condition = "input.ChangeChoice == 'log2 Fold Change'",
               p("Genes with infinite fold changes will not be graphed. Include the following infinite fold change values:"),
               radioButtons("showInf", NULL, c("All", "Both Only"), selected = "Both Only", inline = T),
               checkboxInput("posInfLoss", label = "+Inf for Loss", value = F),
               checkboxInput("negInfLoss", label = "-Inf for Loss", value = F),
               checkboxInput("posInfGain", label = "+Inf for Gain", value = F),
               checkboxInput("negInfGain", label = "-Inf for Gain", value = F)
             ),
             dataTableOutput("DataTable")
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
             dataTableOutput("extraPoints")
             ), windowTitle = "DEACT"
  )
)