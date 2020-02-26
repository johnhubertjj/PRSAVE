# Define UI

Gene_set_UI <- function(id) {
  ns <- NS(id)
  
tabPanel("Gene-set Analysis Viewer",
         # Sidebar layout with input and output definitions ----
         
         sidebarLayout(
           
           ## Write shiny UI across 4 parameters in the data table
           
           sidebarPanel(
             fileInput(ns("file1"), "Choose an input file",
                       multiple = F),
             radioButtons("step_regress", label = "Add stepwise regression for DSM grouping?",
                          c("Yes", "No"),
                          selected = "No"),
             uiOutput(ns("Significance_threshold")
             ),
             checkboxGroupInput("Gene_regions", label = "Length of Gene regions:",
                                choices = c("extended","normal","Full")),
             uiOutput(ns("DSM")),
             uiOutput(ns("geneset")),
             sliderInput('plotHeight', 'Bar which does nothing, use if bored', 
                         min = 100, max = 2000, value = 1000),
             actionButton("add", "Add 'Dynamic' tab"),
             actionButton("remove", "Remove 'Foo' tab")
           ),
           
           mainPanel(
             tabsetPanel(id = "tabs",
                         
                         tabPanel("Plots",
                                  conditionalPanel(
                                    condition = "input.tabs == 'Plots' & input.step_regress == 'No'",
                                    plotOutput(ns('PvalPlot')),
                                    plotOutput(ns('Beta_plot')),
                                    plotOutput(ns('R2_plot')),
                                  conditionalPanel(condition = "input.tabs == 'Plots' & input.step_regress == 'Yes'",
                                                   plotOutput(ns('PvalPlot_step')),
                                                   plotOutput(ns('Beta_plot_step')),
                                                   plotOutput(ns('R2_plot_step')))),
                         tabPanel("Table", dataTableOutput(ns("summary_table"))),
                         tabPanel("Input variables",
                                  rclipboardSetup(),
                                  
                                  
                                  textAreaInput("text_2",label = "Full message will appear here:",width = "500px",height = "100px", resize = "both",
                                                placeholder = "Twitter handles will appear here at the end of your message depending on the options selected to the left (eg: Pint of Science is Great! @virustinkerer)")
                                  ,
                                  # UI ouputs for the copy-to-clipboard buttons
                                  uiOutput(ns("clip")))
             )
           )
         )
)
)
}


# Define function to create UI interfaces




# Define function to create the tab
Gene_set_PRS <- function(input, output, session, data, part_2){
  
  output$Significance_threshold <- renderUI({ 
    significance_threshold.input <- data()$significance_threshold.input
    checkboxGroupInput("Significance_threshold", label = "PRS P Value Threshold:",
                       choices = significance_threshold.input, selected = significance_threshold.input)
  })
  
  output$DSM <- renderUI({ 
    DSM.input <- data()$DSM.input
    selectInput("DSM", label = "DSM type:",
                choices = DSM.input)
  })
  
  output$geneset <- renderUI({ 
    Gene.sets.input <- data()$Gene.sets.input
    checkboxGroupInput("geneset", label = "Geneset PRS to include:",
                       choices = Gene.sets.input, selected = Gene.sets.input)
  })
  
  observeEvent(input$add, {
    insertTab(inputId = "tabs",
              tabPanel("Dynamic", "This a dynamically-added tab"),
              target = "Table"
    )
  })
  
  observeEvent(input$remove, {
    removeTab(inputId = "tabs", target = "bar")
  })
  
  ##### TAB1  ####  
  output$PvalPlot <- renderPlot({
    
    data()
    
    if (is.null(input$Significance_threshold)) {
      return(NULL)
    }    
    if (is.null(input$geneset)) {
      return(NULL)
    }    
    if (is.null(input$Gene_regions)) {
      return(NULL)
    }    
    
    Sample_analysis_2 <- as.data.table(part_2())
    Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
    Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
    Sample_analysis_2
    
    
    # Plot the resulting table for comparisons
    p <- ggplot(Sample_analysis_2, aes(x=score, y=logp, fill = Type, group=Significance_thresholds))
    
    p <- p +
      geom_point(aes(colour = Type))
    
    
    p <- p + scale_x_discrete(labels=levels(Sample_analysis_2$alterations))
    
    p <- p + facet_grid(. ~ as.double(Significance_thresholds),scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    
    p <- p + geom_hline(aes(yintercept=1.30103), colour = "red", linetype= "solid", alpha = 0.25)
    p <- p + scale_fill_brewer(palette = "Paired")
    p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
    p <- p + ggtitle(Sample_analysis_2$.id[1])
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
    
  })
  output$Beta_plot <- renderPlot({
    
    data()
    if (is.null(input$Significance_threshold)) {
      return(NULL)
    }    
    if (is.null(input$geneset)) {
      return(NULL)
    }    
    if (is.null(input$Gene_regions)) {
      return(NULL)
    }    
    
    Sample_analysis_2 <- as.data.table(part_2())
    Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
    Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
    Sample_analysis_2
    
    # Put in the code below above, removing all of the excess alterations work to create the pdf plots...
    
    p <- ggplot(Sample_analysis_2, aes(x=score, y=estimate, fill = Type, group=Significance_thresholds))
    
    p <- p +
      geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25) +
      geom_point(aes(colour = Type))
    
    p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
    p <- p + facet_grid(. ~ as.double(Significance_thresholds), scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
    p <- p + scale_fill_brewer(palette = "Paired")
    p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
    p <- p + ggtitle(Sample_analysis_2$.id[1])
    p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
    p <- p + ylab(label = "BETA")
    p <- p + xlab(label = "Polygenic risk score")
    p
    
    # ggplotly(p) %>% 
    # layout(height = input$plotHeight, autosize=TRUE)
    
    # Possible improvements:
    # Implement in switch from whole genome to gene-sets
    # Implement data-table of the raw results
    # Implement output file of the plots
    # Colour rows for significant values
    # Incorporate into its own app
    # 
    
  })
  output$R2_plot <- renderPlot({
    
    data()
    
    if (is.null(input$Significance_threshold)) {
      return(NULL)
    }    
    if (is.null(input$geneset)) {
      return(NULL)
    }    
    if (is.null(input$Gene_regions)) {
      return(NULL)
    }    
    
    
    Sample_analysis_2 <- as.data.table(part_2())
    Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
    Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
    Sample_analysis_2
    
    p <- ggplot(Sample_analysis_2, aes(x=score, y=r2_dir, fill = Type, group=Significance_thresholds))
    p <- p +
      geom_bar(stat = "identity", aes(colour = Type), position = "dodge") +
      geom_text(data=subset(Sample_analysis_2, p < 0.05),
                aes(x=score,y=r2_dir,label=p_value_text, hjust=ifelse(sign(r2_dir)>0, 0, 0)), angle = 90, position = position_dodge(width = 1), size = 2.9)
    
    #Problem with labels with a workaround
    # I use the score column in the format of factors and reference each relevant dataset for ggplot.
    # However this relies on having 0.05 and 0.5 in the value name.
    # scale_x_discrete accepts functions, but I also need to convert SCORE_0.05 and Score_0.5 into a "Whole_genome_PRS" which is almost impossible to write"
    # However as the labels function accepts key:value pairs, I wrote a vector in R that maps the original names of the pathways to "human readable" format using names function in R
    # This should work for most instances
    
    p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
    p <- p + scale_y_continuous(expand = expand_scale(mult = c(0.2,.6)))
    p <- p + facet_grid(. ~ as.double(Significance_thresholds), scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    p <- p + theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1,vjust = 0.5))
    p <- p + ggtitle(Sample_analysis_2$.id[1])
    p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
    p <- p + ylab(label = "R2_dir (%)")
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
    
  })
  



output$PvalPlot_step <- renderPlot({
  
  data()
  
  if (is.null(input$Significance_threshold)) {
    return(NULL)
  }    
  if (is.null(input$geneset)) {
    return(NULL)
  }    
  if (is.null(input$Gene_regions)) {
    return(NULL)
  }    
  
  
  Sample_analysis_2 <- as.data.table(part_2())
  Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
  Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
  Sample_analysis_2
  
  # Plot the resulting table for comparisons
  p <- ggplot(Sample_analysis_2, aes(x=score, y=logp, fill = Type, group=Significance_thresholds))
  
  p <- p +
    geom_point(aes(colour = Type))
  
  p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
  
  p <- p + facet_grid(. ~ as.double(Significance_thresholds),scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(size = 10))
  
  p <- p + geom_hline(aes(yintercept=1.30103), colour = "red", linetype= "solid", alpha = 0.25)
  p <- p + scale_fill_brewer(palette = "Paired")
  p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
  p <- p + ggtitle(Sample_analysis_2$.id[1])
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
  
})
output$Beta_plot_step <- renderPlot({
  
  data()
  
  if (is.null(input$Significance_threshold)) {
    return(NULL)
  }    
  if (is.null(input$geneset)) {
    return(NULL)
  }    
  if (is.null(input$Gene_regions)) {
    return(NULL)
  }    
  
  
  Sample_analysis_2 <- as.data.table(part_2())
  Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
  Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
  Sample_analysis_2
  
  
  if(any(grep(input$DSM, Sample_analysis_2$predictors_retained))){
    rows_to_be_highlighted <- grep(input$DSM, Sample_analysis_2$predictors_retained)
    Sample_analysis_2[, colour_change := "black"]
    Sample_analysis_2[rows_to_be_highlighted,colour_change := "green"]
    testing_sample_analysis <- Sample_analysis_2[order(Sample_analysis_2$Significance_thresholds,Sample_analysis_2$alterations)]
    d1 <- split(testing_sample_analysis$colour_change, ceiling(seq_along(testing_sample_analysis$colour_change)/length(levels(testing_sample_analysis$alterations))))
    
    # Put in the code below above, removing all of the excess alterations work to create the pdf plots...
    
    p <- ggplot(Sample_analysis_2, aes(x=score, y=estimate, fill = Type, group=Significance_thresholds))
    
    p <- p +
      geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25) +
      geom_point(aes(colour = Type))
    
    p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
    p <- p + facet_grid(. ~ as.double(Significance_thresholds), scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
    p <- p + scale_fill_brewer(palette = "Paired")
    p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
    p <- p + ggtitle(Sample_analysis_2$.id[1])
    p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
    p <- p + ylab(label = "BETA")
    p <- p + xlab(label = "Polygenic risk score")
    pt <- ggplotGrob(p)
    
    integer <- length(input$Significance_threshold)*2 +1
    for (i in 1:length(input$Significance_threshold)){
      pt$grobs[[integer+i]]$children[[2]]$grobs[[2]]$children[[1]]$gp$col <- d1[[i]]
    }
    
    grid::grid.draw(pt)
  }else{
    p <- ggplot(Sample_analysis_2, aes(x=score, y=estimate, fill = Type, group=Significance_thresholds))
    
    p <- p +
      geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25) +
      geom_point(aes(colour = Type))
    
    p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
    p <- p + facet_grid(. ~ as.double(Significance_thresholds), scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
    p <- p + scale_fill_brewer(palette = "Paired")
    p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
    p <- p + ggtitle(Sample_analysis_2$.id[1])
    p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
    p <- p + ylab(label = "BETA")
    p <- p + xlab(label = "Polygenic risk score")
    p
  }
  # ggplotly(p) %>% 
  # layout(height = input$plotHeight, autosize=TRUE)
  
  # Possible improvements:
  # Implement in switch from whole genome to gene-sets
  # Implement data-table of the raw results
  # Implement output file of the plots
  # Colour rows for significant values
  # Incorporate into its own app
  # 
  
})
output$R2_plot_step <- renderPlot({
  
  data()
  
  if (is.null(input$Significance_threshold)) {
    return(NULL)
  }    
  if (is.null(input$geneset)) {
    return(NULL)
  }    
  if (is.null(input$Gene_regions)) {
    return(NULL)
  }    
  
  
  Sample_analysis_2 <- as.data.table(part_2())
  Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
  Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
  Sample_analysis_2
  
  p <- ggplot(Sample_analysis_2, aes(x=score, y=r2_dir, fill = Type, group=Significance_thresholds))
  p <- p +
    geom_bar(stat = "identity", aes(colour = Type), position = "dodge") +
    geom_text(data=subset(Sample_analysis_2, p < 0.05),
              aes(x=score,y=r2_dir,label=p_value_text, hjust=ifelse(sign(r2_dir)>0, 0, 0)), angle = 90, position = position_dodge(width = 1), size = 2.9)
  
  #Problem with labels with a workaround
  # I use the score column in the format of factors and reference each relevant dataset for ggplot.
  # However this relies on having 0.05 and 0.5 in the value name.
  # scale_x_discrete accepts functions, but I also need to convert SCORE_0.05 and Score_0.5 into a "Whole_genome_PRS" which is almost impossible to write"
  # However as the labels function accepts key:value pairs, I wrote a vector in R that maps the original names of the pathways to "human readable" format using names function in R
  # This should work for most instances
  
  p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
  p <- p + scale_y_continuous(expand = expand_scale(mult = c(0.2,.6)))
  p <- p + facet_grid(. ~ as.double(Significance_thresholds), scales = "free_x", space = "free_x") +
    theme(strip.text.x = element_text(size = 10))
  p <- p + theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1,vjust = 0.5))
  p <- p + ggtitle(Sample_analysis_2$.id[1])
  p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
  p <- p + ylab(label = "R2_dir (%)")
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
  
}) 

# TAB2
output$summary_table <- renderDataTable({
  
  My_data()
  
  # These are required in case no tick boxes are selected
  if (is.null(input$Significance_threshold)) {
    return(NULL)
  }    
  if (is.null(input$geneset)) {
    return(NULL)
  }    
  if (is.null(input$Gene_regions)) {
    return(NULL)
  }    
  
  # Select columns you wish to output
  cols <- c("estimate", "SE","r.squared","p")
  
  ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
  sample_analysis <- My_data()$Full_data %>%
    filter(samples.i. == input$DSM,
           Gene_regions %in% input$Gene_regions,
           Significance_thresholds %in% input$Significance_threshold,
           Genesets %in% input$geneset
    )  %>%
    select(c(score:SE,p,r.squared)) %>%
    arrange(p)
  
  ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
  Sample_analysis_2 <- as.data.table(sample_analysis)
  Sample_analysis_2[, (cols) := lapply(.SD, formatC, digits = 3, format = "g"), .SDcols = cols]
  
  ## leave datatable function in for "prettyfying" the final result    
  datatable(data = Sample_analysis_2,
            options = list(pageLength = 10),
            rownames = F)
  
  # Possible improvements:
  # colour rows for significant values
  # incorporate into its own app
  # 
})

# TAB 3
observe({
  
  
  if (!is.null(input$Significance_threshold)) {
    Sig_print <- paste(input$Significance_threshold, collapse = ",")
    
  } else{
    Sig_print <- "placeholder"
  }  
  
  if (!is.null(input$geneset)) {
    geneset_print <- paste("\"",input$geneset,"\"", collapse = ",", sep = "")
  }else{
    geneset_print <- "placeholder"
  }    
  
  if (!is.null(input$Gene_regions)) {
    Generegion_print <- paste("\"",input$Gene_regions,"\"", collapse = ",", sep="")
  } else{
    Generegion_print <- "placeholder"
  }   
  
  if (!is.null(input$DSM)) {
    DSM_print <- paste0("\"",input$DSM,"\"")
  }else{
    DSM_print <- "placeholder"
  }
  
  if (!is.null(input$file1$datapath)){
    data_path_print <- paste(input$file1$datapath)
  }else{
    data_path_print <- "path_not_found"
  }
  
  output_text <- paste0( "input <- list() \n",
                         "\n input$Significance_threshold <- c(", Sig_print,")",
                         "\n input$geneset <- c(", geneset_print,")", 
                         "\n input$Gene_regions <- c(", Generegion_print,")", 
                         "\n input$DSM <- ", DSM_print, 
                         "\n data_print_path <- ", data_path_print)
  
  
  
  updateTextInput(session,"text_2",value = output_text)
  
  #number_of_characters <- paste(text_output_speaker_2, collapse = " ")
  #number_of_characters <- nchar(input$text_2)
  #output$length_text_left <- renderText(280 - number_of_characters)
})

output$clip <- renderUI({
  rclipButton("clipbtn", "Copy to Clipboard", input$text_2, icon("clipboard"))
})
}