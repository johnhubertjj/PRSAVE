#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize=1000*1024^2) 
  source("reactivity_scripts.R", local = TRUE)
  
  
  output$Significance_threshold <- renderUI({ 
    significance_threshold.input <- as.numeric(My_data()$significance_threshold.input)
    selectInput("Significance_threshold", label = "PRS P Value Threshold:",
                       choices = significance_threshold.input[order(significance_threshold.input)], selected = significance_threshold.input)
  })
  
  output$geneset <- renderUI({ 
    Gene.sets.input <- My_data()$Gene.sets.input
    selectInput("geneset", label = "Geneset PRS to include:",
                       choices = Gene.sets.input, selected = Gene.sets.input)
  })


output$Centile_plot <- renderPlot({
  
  My_data()
  
  if (is.null(input$Significance_threshold)) {
    return(NULL)
  }    
  if (is.null(input$geneset)) {
    return(NULL)
  }    
  if (is.null(input$Gene_regions)) {
    return(NULL)
  }    
  
  p <- ggplot(Part_2()$Result_table, aes(x= Percentile, y = Odds_Ratio, ymin = Minimum, ymax = Maximum)) 
  #p <- ggplot(data = Result_table, aes(x= Decile, y = Odds_Ratio)) 
  p <- p + geom_point()
  p <- p + geom_pointrange(aes(colour = cut(Odds_Ratio, c(-Inf, 2, 5, 10, Inf))))
  p <- p + scale_x_discrete(breaks = seq(0,100, by = 5))
  p <- p + ylab("Odds Ratio")
  
  Title_name <- gsub("_", " ", input$geneset)
  p <- p + ggtitle(paste0(Title_name))
  
  p <- p + scale_color_manual(name = "Odds Ratio", 
                              values = c("(-Inf,2]" = "black",
                                         "(2,5]" = "#e98b56",
                                         "(5,10]" = "#e9d556", 
                                         "(10, Inf]" = "#D72340"),
                              labels = c("<= 2", "2 < Odds Ratio <= 5","5 < Odds Ratio <= 10", "> 10"))
  p
})


#observe({
  # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
  # were a base graphics plot, we'd need those.
#  DF_near_points <- nearPoints(Part_2()$Result_table, input$plot1_click, addDist = TRUE, xvar = "Percentile", yvar = "Odds_Ratio")
#  DF_near_points_PRS <- Part_2()$PRS_table
#  setnames(DF_near_points, old = "Percentile" , new = paste0(input$Gene_regions,"_geneset_SCORE_",input$geneset,"_",input$Significance_threshold,"_decile"))
#  new_DF <- merge(DF_near_points,DF_near_points_PRS, by = paste0(input$Gene_regions,"_geneset_SCORE_",input$geneset,"_",input$Significance_threshold,"_decile"))
  
#  output$click_info <- renderDT({
#    new_DF
#  })
  
#})
observe({
output$brush_info <- renderDataTable({
  df <- brushedPoints(Part_2()$Result_table, input$plot1_brush, xvar = "Percentile", yvar = "Odds_Ratio")
  df$Percentile <- as.numeric(as.character(df$Percentile))
  df <- as.data.table(df)
  df_PRS <- Part_2()$PRS_table
  #colnames(df_PRS)
  index_PRS <- grep("decile", x = names(df_PRS))
  df_PRS2 <- df_PRS[,index_PRS, with = F]
  setnames(df_PRS2, "Percentile")
  test <- match(df_PRS2$Percentile,df$Percentile )
  test2 <- !is.na(test)
  df_PRS[test2]
  
  #setnames(df_PRS, old = paste0(input$Gene_regions,"_geneset_SCORE_",input$geneset,"_",input$Significance_threshold,"_decile"), new = "Odds_ratio")
  
  
  #df_PRS[ df_PRS[,index_PRS, with = F] %in% df2$Percentile]
})

output$downloadData <- downloadHandler(
  
  # This function returns a string which tells the client
  # browser what name to use when saving the file.
  filename = function() {
    filename_input <- paste0(input$Gene_regions,"_geneset_SCORE_",input$geneset,"_",input$Significance_threshold,"_decile")
    paste(filename_input, input$filetype, sep = ".")
  },
  
  # This function should write data to a file given to it by
  # the argument 'file'.
  content = function(file) {
    sep <- switch(input$filetype, "csv" = ",", "tsv" = "\t")
    
    # Write to a file specified by the 'file' argument
    write.table(df_PRS[test2], file, sep = sep,
                row.names = FALSE)
  }
)
})



})
