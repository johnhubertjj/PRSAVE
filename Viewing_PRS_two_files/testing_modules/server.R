#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  #source("reactivity_scripts.R", local = TRUE)
  source("PRSicely_reactive_input_from_PRSet.R", local = TRUE)
  
  output$Significance_threshold <- renderUI({ 
    significance_threshold.input <- as.numeric(My_data()$significance_threshold.input)
    checkboxGroupInput(session$ns("Significance_threshold"), label = "PRS P Value Threshold:",
                       choices = significance_threshold.input[order(significance_threshold.input)], selected = significance_threshold.input)
  })
  
  output$DSM <- renderUI({ 
    DSM.input <- My_data()$DSM.input
    selectInput(session$ns("DSM"), label = "DSM type:",
                choices = DSM.input)
  })
  
  output$geneset <- renderUI({ 
    Gene.sets.input <- My_data()$Gene.sets.input
    checkboxGroupInput(session$ns("geneset"), label = "Geneset PRS to include:",
                       choices = Gene.sets.input, selected = Gene.sets.input)
  })
  
  # debouncing algorithm -> MUST GO HERE because of the rendering UI scripts above and reactivity scripts have not been processed yet
  sigthreshold_debounce <- reactive({ input$Significance_threshold }) %>% debounce(1000)
  gene_set_debounce <- reactive({ input$geneset }) %>% debounce(1000)
  Gene_region_debounce <- reactive({ input$Gene_regions }) %>% debounce(1000)
  
  pvalplot_modules <- callModule(results_mod, "data_file")
  
  output$PvalPlot <- renderPlot({
    pvalplot_modules()
  
  })
  
  })
