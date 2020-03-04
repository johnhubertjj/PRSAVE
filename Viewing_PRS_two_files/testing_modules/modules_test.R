# test_script_modules

page_one_UI <- function(id, Label = "page_one_UI") {
  
  ns <- NS(id)
  

               tagList(
               fileInput("file1", "Choose an input file",
                         multiple = F),
               uiOutput(ns("Significance_threshold")
               ),
               checkboxGroupInput("Gene_regions", label = "Length of Gene regions:",
                                  choices = c("Genome-wide","Gene-set")),
               uiOutput(ns("DSM")),
               uiOutput(ns("geneset"))
               )
               
}             
            

                                    
                              
                                    
                                    
                               
                                    
                      

results_mod <- function(input,output, session) {
  


##### TAB1  ####  
output$PvalPlot <- renderPlot({
  
  My_data()
  
  if (is.null(sigthreshold_debounce())) {
    return(NULL)
  }    
  if (is.null(gene_set_debounce())) {
    return(NULL)
  }    
  if (is.null(Gene_region_debounce())) {
    return(NULL)
  }    
  
  
  
  # Plot the resulting table for comparisons
  p <- ggplot(part_2(), aes(x=score, y=logp, fill = Type, group=Significance_thresholds))
  
  p <- p +
    geom_point(aes(colour = Type))
  
  
  p <- p + scale_x_discrete(labels=levels(part_2()$alterations))
  
  p <- p + facet_grid(. ~ as.double(Significance_thresholds),scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(size = 10))
  
  p <- p + geom_hline(aes(yintercept=1.30103), colour = "red", linetype= "solid", alpha = 0.25)
  p <- p + scale_fill_brewer(palette = "Paired")
  p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
  p <- p + ggtitle(part_2()$.id[1])
  #p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5,vjust = 0.5,angle = 90))
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  p <- p + ylab(label = expression(-'log'[10]*'(p)'))
  p <- p + xlab(label = "Polygenic risk score")
  p
  #ggplotly(p) %>% 
  #  layout(height = input$plotHeight, autosize=TRUE)
  
  # Possible improvements:
  # Implement in switch from whole genome to gene-sets
  # Implement data-table of the raw results
  # Implement output file of the plots
  # Colour rows for significant values
  # Incorporate into its own app
  # 
  
})}