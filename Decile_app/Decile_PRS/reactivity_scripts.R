
# Actually, it only requires one drop-down box here...

My_data <- reactive({ req(input$file1)
  
  ## Read in data
  Full_data <- fread(input$file1$datapath)
  #Full_data <- fread("~/Documents/shiny_PRS/Decile_app/Decile_PRS/centile_table_example.txt")
  
  ## Create new columns parsing the identifiers in the Full_data score column and input to the shiny app
  colnames_PRS <- colnames(Full_data)
  genesets_index <- grep(pattern = ".*geneset_SCORE_(.*)_.*", x = colnames_PRS)
  whole_genome_genic_positions <- grep(x = Full_data$Gene_regions, pattern = "genic.genome_SCORE_whole_genome",perl = T)
  
  genesets_all <- colnames_PRS[genesets_index]
  
  ## Create arguments to shiny app
  Gene.sets.input <- unique(gsub(pattern = ".*_SCORE_(.*)_.*", replacement = "\\1", x = genesets_all, perl = T))
  Gene.sets.input <- c(Gene.sets.input,"genic.genome","All.genome")
  Gene_regions <- unique(gsub(pattern = "^(.*)_geneset_SCORE_.*", replacement = "\\1", x = genesets_all, perl = T))
  significance_threshold.input <- unique(gsub(pattern = ".*_(.*$)", replacement = "\\1", x = genesets_all, perl = T))
  
  
  Full_data <- list(Full_data = Full_data ,Gene.sets.input = Gene.sets.input, significance_threshold.input = significance_threshold.input, Gene_regions.input = Gene_regions)
  Full_data# Potential solution to multiple
})

Part_2 <- reactive({ 
  
 Current_table <- My_data()$Full_data
#Current_table <- Full_data$Full_data
  
  index1 <- grep(input$geneset, x = names(Current_table))
  pheno <- grep("PHENO|FID|IID", x = names(Current_table)) 
  index_pheno <- c(index1, pheno)
  Current_table <- Current_table[,index_pheno, with = F]

  index2 <- grep(input$Significance_threshold, x = names(Current_table))
  pheno <- grep("PHENO|FID|IID", x = names(Current_table)) 
  index_pheno <- c(index2, pheno)
  Current_table <- Current_table[,index_pheno, with = F]
  
  index3 <- grep(input$Gene_regions, x = names(Current_table))
  pheno <- grep("PHENO|FID|IID", x = names(Current_table)) 
  index_pheno <- c(index3, pheno)
  Current_table <- Current_table[,index_pheno, with = F]
  
  Columns_use_shiny_names <- paste0(input$Gene_regions,"_geneset_SCORE_",input$geneset,"_",input$Significance_threshold)
  Columns_use_shiny_names_decile <- paste0(Columns_use_shiny_names, "_decile")
  
  Current_table[, c(Columns_use_shiny_names_decile) := lapply(.SD, function(x){ntile(x, 100)}), .SDcols = Columns_use_shiny_names]
  if (all(Current_table$PHENO == 1 | Current_table$PHENO == 0) == FALSE){
    Current_table[,PHENO := PHENO - 1]
  }
  
  Result_table <- matrix(data = rep(0,30), nrow = 100, ncol = 3)
  Result_table[1,] <- c(1,0,0)
  
  
  for (i in 2:100){
      All_deciles_controls <- Current_table[PHENO == 0, .(eval(parse( text = paste0(input$Gene_regions,"_geneset_SCORE_",input$geneset,"_",input$Significance_threshold,"_decile"))))]
      All_deciles_cases <- Current_table[PHENO == 1, .(eval(parse(text = paste0(input$Gene_regions,"_geneset_SCORE_",input$geneset,"_",input$Significance_threshold,"_decile"))))]

    
    All_deciles_controls_equals_one <- length(which(All_deciles_controls == 1))
    All_deciles_controls_equals_two <- length(which(All_deciles_controls == i))
    All_deciles_cases_equals_one <- length(which(All_deciles_cases == 1))
    All_deciles_cases_equals_two <- length(which(All_deciles_cases == i))
    
    test <- matrix(data = c(All_deciles_cases_equals_two, All_deciles_controls_equals_two, All_deciles_cases_equals_one, All_deciles_controls_equals_one), byrow = T, ncol = 2, nrow = 2)
    test_result <- oddsratio(test)
    Result_table[i,] <- test_result$measure[2,]
  }
  
  Result_table <- as.data.frame(Result_table)
  Result_table <- cbind(c(1:100), Result_table)
  colnames(Result_table) <- c("Percentile","Odds_Ratio", "Minimum", "Maximum")
  Result_table$Percentile <- as.factor(Result_table$Percentile)
  
  Result_table <- list(Result_table = Result_table, PRS_table = Current_table)
})


