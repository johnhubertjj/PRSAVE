#ModulesUI

innerModUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(fluidRow(
    uiOutput(ns("inner_slider")),
    plotOutput(ns("inner_plot"))
  ))
}

page_one_UI <- function(id) {
  
  ns <- NS(id)

  tabPanel("Gene-set Analysis Viewer",
           # Sidebar layout with input and output definitions ----
           
           sidebarLayout(
             
             ## Write shiny UI across 4 parameters in the data table
             sidebarPanel(
               fileInput("file1", "Choose an input file",
                         multiple = F),
               radioButtons("step_regress", label = "Add stepwise regression for DSM grouping?",
                            c("Yes", "No"),
                            selected = "No"),
               uiOutput(ns("Significance_threshold")
               ),
               checkboxGroupInput("Gene_regions", label = "Length of Gene regions:",
                                  choices = c("Genome-wide","Gene-set")),
               uiOutput(ns("DSM")),
               uiOutput(ns("geneset")),
               sliderInput('plotHeight', 'Bar which does nothing, use if bored', 
                           min = 100, max = 2000, value = 1000)
             ),
             
             mainPanel(
               tabsetPanel(id = "tabs",
                           
                           tabPanel("Plots",
                                    conditionalPanel(
                                      condition = "input.tabs == 'Plots' & input.step_regress == 'No'",
                                      plotOutput(ns('PvalPlot')),
                                      plotOutput(ns('Beta_plot')),
                                      plotOutput(ns('R2_plot'))),
                                    conditionalPanel(condition = "input.tabs == 'Plots' & input.step_regress == 'Yes'",
                                                     plotOutput(ns('PvalPlot_step')),
                                                     plotOutput(ns('Beta_plot_step')),
                                                     plotOutput(ns('R2_plot_step')))),
                           tabPanel("Table", dataTableOutput(ns("summary_table"))),
                           tabPanel("Input variables",
                                    rclipboardSetup(),
                                    
                                    
                                    textAreaInput(ns("text_2"),label = "Full message will appear here:",width = "500px",height = "100px", resize = "both",
                                                  placeholder = "Twitter handles will appear here at the end of your message depending on the options selected to the left (eg: Pint of Science is Great! @virustinkerer)")
                                    ,
                                    # UI ouputs for the copy-to-clipboard buttons
                                    uiOutput(ns("clip")))
               )
             )
           )
  )
}