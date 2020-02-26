#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Percentile analysis"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose an input file",
                multiple = F),
      selectInput("Gene_regions", label = "Gene regions:",
                  choices = c("extended","normal")),
      uiOutput("Significance_threshold"
      ),

      uiOutput("geneset"),
      sliderInput('plotHeight', 'Bar which does nothing, use if bored', 
                  min = 100, max = 2000, value = 1000),
      radioButtons("filetype", "File type:",
                   choices = c("csv", "tsv")),
      downloadButton('downloadData', 'Download')
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
       plotOutput("Centile_plot",click = "plot1_click", brush = brushOpts(
         id = "plot1_brush"
       )),
    
    fluidRow(
      column(width = 6,
             h4("Points near click"),
             verbatimTextOutput("click_info")
      ),
      column(width = 6,
             h4("Brushed points"),
             dataTableOutput("brush_info")
      )
    )
  )
))
)
