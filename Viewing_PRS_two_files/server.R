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

options(shiny.maxRequestSize=30*1024^2)
#source("reactivity_scripts.R", local = TRUE)
source("PRSicely_reactive_input_from_PRSet.R", local = TRUE)
  
    output$Significance_threshold <- renderUI({ 
      significance_threshold.input <- as.numeric(My_data()$significance_threshold.input)
      checkboxGroupInput("Significance_threshold", label = "PRS P Value Threshold:",
                         choices = significance_threshold.input[order(significance_threshold.input)], selected = significance_threshold.input)
    })
    
    output$DSM <- renderUI({ 
      DSM.input <- My_data()$DSM.input
      selectInput("DSM", label = "DSM type:",
                  choices = DSM.input)
    })
    
    output$geneset <- renderUI({ 
      Gene.sets.input <- My_data()$Gene.sets.input
      checkboxGroupInput("geneset", label = "Geneset PRS to include:",
                         choices = Gene.sets.input, selected = Gene.sets.input)
    })
    
  # debouncing algorithm -> MUST GO HERE because of the rendering UI scripts above and reactivity scripts have not been processed yet
    sigthreshold_debounce <- reactive({ input$Significance_threshold }) %>% debounce(1000)
    gene_set_debounce <- reactive({ input$geneset }) %>% debounce(1000)
    Gene_region_debounce <- reactive({ input$Gene_regions }) %>% debounce(1000)
    
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
      
    })
  output$Beta_plot <- renderPlot({
    
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
    
      
      # Put in the code below above, removing all of the excess alterations work to create the pdf plots...
      
      p <- ggplot(part_2(), aes(x=score, y=estimate, fill = Type, group=Significance_thresholds))
      
      p <- p +
        geom_errorbar(aes(ymin = SE_higher, ymax = SE_lower), position = "dodge", width = 0.25) +
        geom_point(aes(colour = Type))
      
      p <- p + scale_x_discrete(labels= levels(part_2()$alterations))
      p <- p + facet_grid(. ~ as.double(Significance_thresholds), scales = "free_x", space = "free_x") +
        theme(strip.text.x = element_text(size = 10))
      p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
      p <- p + scale_fill_brewer(palette = "Paired")
      p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
      p <- p + ggtitle(part_2()$.id[1])
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
    
      p <- ggplot(part_2(), aes(x=score, y=r2_dir, fill = Type, group=Significance_thresholds))
      p <- p +
        geom_bar(stat = "identity", aes(colour = Type), position = "dodge") +
        geom_text(data=subset(part_2(), P < 0.05),
                  aes(x=score,y=r2_dir,label=p_value_text, hjust=ifelse(sign(r2_dir)>0, 0, 0)), angle = 90, position = position_dodge(width = 1), size = 2.9)
      
      #Problem with labels with a workaround
      # I use the score column in the format of factors and reference each relevant dataset for ggplot.
      # However this relies on having 0.05 and 0.5 in the value name.
      # scale_x_discrete accepts functions, but I also need to convert SCORE_0.05 and Score_0.5 into a "Whole_genome_PRS" which is almost impossible to write"
      # However as the labels function accepts key:value pairs, I wrote a vector in R that maps the original names of the pathways to "human readable" format using names function in R
      # This should work for most instances
      
      p <- p + scale_x_discrete(labels= levels(part_2()$alterations))
      p <- p + scale_y_continuous(expand = expand_scale(mult = c(0.2,.6)))
      p <- p + facet_grid(. ~ as.double(Significance_thresholds), scales = "free_x", space = "free_x") +
        theme(strip.text.x = element_text(size = 10))
      p <- p + theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1,vjust = 0.5))
      p <- p + ggtitle(part_2()$.id[1])
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
    
    ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
    Full_data <- My_data()$Full_data %>%
      filter(samples.i. == input$DSM,
             Gene_regions %in% input$Gene_regions,
             Significance_thresholds %in% input$Significance_threshold,
             Genesets %in% input$geneset)
    
    
    
    ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
    
    ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
    ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
    ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.
    
    Sample_analysis_2 <- as.data.table(Full_data)
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
    
    ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
    Full_data <- My_data()$Full_data %>%
      filter(samples.i. == input$DSM,
             Gene_regions %in% input$Gene_regions,
             Significance_thresholds %in% input$Significance_threshold,
             Genesets %in% input$geneset)
    
    
    
    ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
    
    ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
    ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
    ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.
    
    Sample_analysis_2 <- as.data.table(Full_data)
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
    
    ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
    Full_data <- My_data()$Full_data %>%
      filter(samples.i. == input$DSM,
             Gene_regions %in% input$Gene_regions,
             Significance_thresholds %in% input$Significance_threshold,
             Genesets %in% input$geneset)
    
    
    
    ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
    
    ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
    ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
    ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.
    
    Sample_analysis_2 <- as.data.table(Full_data)
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
      cols <- c("estimate", "SE","R2","P")
      
      ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
      sample_analysis <- My_data()$Full_data %>%
        filter(samples.i. == input$DSM,
               Gene_regions %in% input$Gene_regions,
               Significance_thresholds %in% input$Significance_threshold,
               Genesets %in% input$geneset
        )  %>%
        select(c(Genesets,Gene_regions,Significance_thresholds,estimate,SE,P,R2,Num_SNP)) %>%
        arrange(P)
      
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
  
##### TAB_2 ####

  output$Significance_threshold_2 <- renderUI({
    significance_threshold_2_input <- as.numeric(My_data_2()$significance_threshold_2_input)
    checkboxGroupInput("Significance_threshold_2", label = "PRS P Value Threshold:",
                       choices = significance_threshold_2_input, selected = "0.05")
  })
  
  output$DSM_2 <- renderUI({
    DSM_2_input <- My_data_2()$DSM_2_input
    selectInput("DSM_2", label = "DSM type:",
                choices = DSM_2_input)
  })
  
  output$GWAS_to_include <- renderUI({
    GWAS_to_choose_from <- My_data_2()$GWAS_to_choose_from
    checkboxGroupInput("GWAS_to_include", label = "Choose GWAS for input into PCA",
                       choices = GWAS_to_choose_from)
  })
  
  
output$PCA_plot  <- renderPlot({
My_data_2() 
  
  if (is.null(input$Significance_threshold_2)) {
    return(NULL)
  }    
  if (is.null(input$GWAS_to_include)) {
    return(NULL)
  }    
  

g <- ggbiplot(Combined_PRS_part_2()$PCA_analysis, ellipse = F, choices = c(1,2), varname.size = input$varname.size, varname.adjust = input$varname.adjust,
              circle = F, alpha = input$alpha)

#g <- g + scale_color_manual(name="CLOZUK Phenotype",label = c("Controls","Cases"), values = c("#E87D72","#4EA8EC")) 
##4EA8EC
#g <- g + scale_colour_brewer(palette = "Set3")
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top',
               legend.box.background = element_blank(),
               legend.background = element_blank(),
               legend.key = element_blank(),
               axis.line = element_line(colour = "black"),
               panel.background = element_blank(),
               legend.title = element_text(face = "bold"))
g
})

output$R2_plot_2 <- renderPlot({
  
  My_data_2() 
  
  if (is.null(input$Significance_threshold_2)) {
    return(NULL)
  }    
  if (is.null(input$GWAS_to_include)) {
    return(NULL)
  }    

p <- ggplot(Combined_PRS_part_2()$rows_of_polygenic_risk_scores_to_use, aes(x=score, y=r2_dir, fill = Number_of_GWAS))
p <- p +
  geom_bar(stat = "identity", aes(colour = Number_of_GWAS), position = "dodge") +
  geom_text(data=subset(Combined_PRS_part_2()$rows_of_polygenic_risk_scores_to_use, p < 0.05),
            aes(x=score,y=r2_dir,label=p_value_text, hjust=ifelse(sign(r2_dir)>0, 0, 0)), angle = 90, position = position_dodge(width = 1), size = 2.9)

#Problem with labels with a workaround
# I use the score column in the format of factors and reference each relevant dataset for ggplot.
# However this relies on having 0.05 and 0.5 in the value name.
# scale_x_discrete accepts functions, but I also need to convert SCORE_0.05 and Score_0.5 into a "Whole_genome_PRS" which is almost impossible to write"
# However as the labels function accepts key:value pairs, I wrote a vector in R that maps the original names of the pathways to "human readable" format using names function in R
# This should work for most instances

p <- p + scale_x_discrete(labels= Combined_PRS_part_2()$rows_of_polygenic_risk_scores_to_use$alterations)
p <- p + scale_y_continuous(expand = expand_scale(mult = c(0.2,.6)))
p <- p + theme(axis.text.x = element_text(angle = 60, size = 10, hjust = 1,vjust = 1))
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


  
  
output$Significance_threshold_used <- renderText({

  My_data_2() 
  
  if (is.null(input$Significance_threshold_2)) {
    return(NULL)
  }    
  if (is.null(input$GWAS_to_include)) {
    return(NULL)
  }    
  
   test <- paste("Significance threshold for ", Combined_PRS_part_2()$rows_of_polygenic_risk_scores_to_use$Number_of_GWAS, " is ", Combined_PRS_part_2()$rows_of_polygenic_risk_scores_to_use$Significance_thresholds, collapse = "\n")
})

output$Number_of_individuals_used <- renderText({
  My_data_2() 
  
  if (is.null(input$Significance_threshold_2)) {
    return(NULL)
  }    
  if (is.null(input$GWAS_to_include)) {
    return(NULL)
  }    
  
    paste0("Number of People considered in analysis = ", nrow(Combined_PRS_part_2()$Polygenic_risk_scores_for_analysis))
})
# output$P <- renderPlot()
# output$Beta_plot_2 <- renderPlot()
# output$R2_plot_2
  ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
  
  ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
  ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
  ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.

  
# Apply groupings to pick out the best threshold from the p-values mentioned and use those profiles in the PCA of the polygenic risk scores
  

####  

  
  
  })
