# Prsice conversion

library(tidyverse)
# Code that will convert the prsice input to a standardised output
# Input lists

Full_data <- fread("~/PRSice_shiny_app/test_run_thesis_phenotype.prsice")
Prsicely_input <- fread("~/Documents/Cognition_paper_automated/Data/Paper_results/Collated_PRS_analysis_forIQ3_COGS.txtCOGS_SCZ_pathway_supersets_analysis_PRS_analysis_run_2_update_both_IQ_and_SCZ_conditional_26-10-2018_at_14.02.txt")

  input <- list() 
  
  input$Significance_threshold <- c(0.2,0.5,1)
  input$geneset <- c("placeholder")
  input$Gene_regions <- c("placeholder")
  input$DSM <- "placeholder"
  input$file1$datapath <- Prsicely_input 


Prsice_conversion <- function(Prcise_output){
  # Reactive script
  
  My_data <- reactive({ req(input$file1)
    
    ## Read in data
    Full_data <- fread(input$file1$datapath)
    
    setnames(Full_data, old = c("Set","Threshold"), new = c("Genesets", "Significance_thresholds"))
    Full_data[, Gene_regions := Genesets]
    
    Genome_wide_PRS <- which(Full_data$Gene_regions == "Base")
    Gene_set_PRS <- which(Full_data$Gene_regions != "Base")
    
    Full_data[,Gene_regions := "NA"]
    Full_data[Genome_wide_PRS, Gene_regions := "Genome-wide"]
    Full_data[Gene_set_PRS, Gene_regions := "Gene-set"]
 
    
    ## Create arguments to shiny app
    Gene.sets.input <- unique(Full_data$Genesets)
    significance_threshold.input <- unique(Full_data$Significance_thresholds)
    DSM.input <- "Everything"
    
    ## Identify which rows in the data table contain whole genome information
    whole_genome_genic_positions_Full_data <- NULL
    whole_genome_all_genome_positions_Full_data <- Genome_wide_PRS
    
    whole_genome_plot_all_positions <- Genome_wide_PRS
    
    
    Full_data[whole_genome_plot_all_positions, Type := "Whole_genome"]
    Full_data[!whole_genome_plot_all_positions, Type:= "Pathway"]
    
# Got here, need to figure out how to manipulate the P-value columns

    if(any(Full_data$p == 0) == T){
      Full_data[,P_altered := p]
      Full_data[P_altered == 0, P_altered := 1e-300]
      Full_data$logp <- -log10(Full_data$P_altered)
    }else{
      Full_data$logp <- -log10(Full_data$p)
    }
    
    Full_data$SE_higher <- Full_data$estimate + Full_data$SE
    Full_data$SE_lower <- Full_data$estimate - Full_data$SE
    Full_data$r2_dir <- 100 * (as.numeric(Full_data$r.squared) *
                                 (sign(as.numeric(Full_data$estimate))))
    Full_data$p_value_text <- paste("p =", scientific(Full_data$p, digits = 2), sep = " ")
    
    if(any(Full_data$p_value_text == "p = 0.0e+00") == TRUE){
      
      Full_data[p_value_text == "p = 0.0e+00", p_value_text := "p < 1e-300"]
    }
    
    ## Add alterations column to create "human readable" formats of the data
    alterations <- Full_data$score
    alterations[-whole_genome_plot_all_positions] <- str_replace(string = alterations[-whole_genome_plot_all_positions], pattern = ".*SCORE_(.*)_.*",replacement = "\\1")
    alterations[-whole_genome_plot_all_positions] <- str_replace(string = alterations[-whole_genome_plot_all_positions], pattern = "long_term_potentiation",replacement = "LTP")
    alterations[-whole_genome_plot_all_positions] <- str_replace(string = alterations[-whole_genome_plot_all_positions], pattern = "action_potential",replacement = "AP")
    alterations[-whole_genome_plot_all_positions] <- str_replace(string = alterations[-whole_genome_plot_all_positions], pattern = "depolarization",replacement = "DP")
    
    # IQ alterations
    alterations[-whole_genome_plot_all_positions] <- tolower(alterations[-whole_genome_plot_all_positions])
    alterations[-whole_genome_plot_all_positions] <- str_replace(string = alterations[-whole_genome_plot_all_positions], pattern = "nervous_system",replacement = "NS")
    alterations[-whole_genome_plot_all_positions] <- str_replace(string = alterations[-whole_genome_plot_all_positions], pattern = "regulation",replacement = "reg")
    
    alterations[whole_genome_genic_positions_Full_data] <- "Whole Genome PRS GENIC"
    alterations[whole_genome_all_genome_positions_Full_data] <- "Whole Genome PRS ALL"
    
    Full_data[, alterations := alterations]
    
    
    Full_data <- list(Full_data = Full_data ,Gene.sets.input = Gene.sets.input, significance_threshold.input = significance_threshold.input, DSM.input = DSM.input)
    Full_data# Potential solution to multiple
  })
  
  part_2 <- reactive({ 
    
    Current_table <- My_data()$Full_data %>%
      
      filter(samples.i. == input$DSM,
             Gene_regions %in% Gene_region_debounce(),
             Significance_thresholds %in% sigthreshold_debounce(),
             Genesets %in% gene_set_debounce())
    
    Sample_analysis_2 <- as.data.table(Current_table)
    Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
    Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
    Sample_analysis_2
    
    ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
    
    ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
    ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
    ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.
    
    
  })
  
  My_data_2 <- reactive({
    
    ## Require input from files
    req(input$file2)
    req(input$file3)
    
    #file <- rep("a",2)
    #file[1] <- paste0(Location_of_analysis,Name_of_analysis,"_",time_of_analysis,".txt")
    #file[2] <- paste0(Location_of_analysis,Name_of_analysis,"_PRS_Profiles_for_shiny",time_of_analysis,".txt")
    #File_two_data <- lapply(file, fread, header=TRUE) 
    
    ## Discern which data table has linear regression information and which has the polygenic risk score profiles
    #File_two_data <- lapply(input$file2$datapath, fread, header=TRUE) 
    #colnames_first <- colnames(File_two_data[[1]])
    
    #if (colnames_first[1] == "FID") { 
    #  Full_data_2 <- File_two_data[[2]]
    #  Polygenic_risk_scores <- File_two_data[[1]]
    #}else{
    #  Full_data_2 <- File_two_data[[1]]
    #  Polygenic_risk_scores <- File_two_data[[2]]
    #}
    
    Full_data_2 <- fread(input$file2$datapath)
    Polygenic_risk_scores <- fread(input$file3$datapath)
    
    Full_data_2[, Significance_thresholds := gsub(pattern = ".*_(.*$)", replacement = "\\1", x = Full_data_2$score, perl = T)]
    
    changeCols <- c("estimate", "SE", "zvalue", "p", "r.squared.Nagel", "lower", "upper","Significance_thresholds")
    Full_data_2[,(changeCols) := lapply(.SD, as.numeric), .SDcols = changeCols]
    
    # Read in significance level input
    significance_threshold_2_input <- unique(Full_data_2$Significance_thresholds)
    DSM_2_input <- unique(Full_data_2$samples.i.)
    
    ## Identify which rows in the data table contain whole genome information
    Full_data_2[, Number_of_GWAS := gsub(pattern = "(^.*)_WHOLE.*", replacement = "\\1", x = Full_data_2$score, perl = T)]
    Number_of_GWAS_vector <- length(unique(Full_data_2$Number_of_GWAS))
    GWAS_to_choose_from <- unique(Full_data_2$Number_of_GWAS)
    
    Full_data_2$logp <- -log10(Full_data_2$p)
    Full_data_2$SE_higher <- Full_data_2$estimate + Full_data_2$SE
    Full_data_2$SE_lower <- Full_data_2$estimate - Full_data_2$SE
    Full_data_2$r2_dir <- 100 * (as.numeric(Full_data_2$r.squared) *
                                   (sign(as.numeric(Full_data_2$estimate))))
    Full_data_2$p_value_text <- paste("p =", scientific(Full_data_2$p, digits = 2), sep = " ")
    
    
    File_two_data <- list (Full_data_2 = Full_data_2,
                           Polygenic_risk_scores = Polygenic_risk_scores,
                           significance_threshold_2_input = significance_threshold_2_input,
                           DSM_2_input = DSM_2_input,
                           GWAS_to_choose_from = GWAS_to_choose_from)
    File_two_data
  })
  
  Combined_PRS_part_2 <- reactive({
    
    Full_data_2 <- My_data_2()$Full_data_2
    Polygenic_risk_scores <- My_data_2()$Polygenic_risk_scores
    
    
    PCA_data <- Full_data_2  %>%
      filter(samples.i. == input$DSM_2,
             Significance_thresholds %in% input$Significance_threshold_2,
             Number_of_GWAS %in% input$GWAS_to_include
      ) %>%
      arrange(Number_of_GWAS,p)
    
    rows_of_polygenic_risk_scores_to_use <- PCA_data [match(unique(PCA_data$Number_of_GWAS), PCA_data$Number_of_GWAS),]
    
    Names_for_PCA_analysis <- c("FID","IID",rows_of_polygenic_risk_scores_to_use$score)
    
    # Here is where the changes are.
    # I need to read in a Polygenic risk score profile with:
    # - PRS for all traits (pre-standardised for each other with population covariates)
    # - Phenotype
    # - covariates for regression (age + sex)
    
    # Here I remove missing values related to the FINAL_TRS phenotype, might need some extra work for other instances
    Polygenic_risk_scores_for_analysis_phenotype_changes <- Polygenic_risk_scores[FINAL_TRS != 9]
    
    # Here I specify what DSM is under each trait, so I remove individuals depending on what the input is selected
    # Not sure if the `if` statement will work here, but I should check
    
    if (input$DSM_2 == "broad"){
      rows_to_keep <- Polygenic_risk_scores_for_analysis_phenotype_changes[DSM_diagnosis_new == 1 | DSM_diagnosis_new == 2 | DSM_diagnosis_new == 3 | DSM_diagnosis_new == 4]
    }
    
    # select out the test for the regression
    rows_to_keep_TRS <- rows_to_keep[,test := FINAL_TRS]
    
    # Segment the table to only columns that you want to keep # Works as standardised code as before
    # Need to feed in the row cut-offs for both DSM and missing values to the initial regression results ASAP, currently does not cut down data properly within the app
    # 
    Polygenic_risk_scores_for_analysis <- rows_to_keep_TRS[,Names_for_PCA_analysis, with = F]
    test1 <- Names_for_PCA_analysis
    test1 <- gsub(pattern = "(^.*)_WHOLE.*", replacement = "\\1",x = test1)
    setnames(Polygenic_risk_scores_for_analysis, old = Names_for_PCA_analysis, new = test1)
    Columns_to_include <- test1[3:length(test1)]
    
    alterations <- rows_of_polygenic_risk_scores_to_use$score
    alterations_2 <- str_replace(string = alterations, pattern = "^(.*?)_.*SCORE_(.*_.*)_(.*)$",replacement = "\\1_\\2_\\3")
    rows_of_polygenic_risk_scores_to_use <- as.data.table(rows_of_polygenic_risk_scores_to_use)
    rows_of_polygenic_risk_scores_to_use[, alterations := alterations_2]
    
    PCA_analysis <- prcomp(Polygenic_risk_scores_for_analysis [,c(Columns_to_include), with = F],center = T)
    
    # Assuming that the PCA rows match the identifiers in the input (PCA matrix has no identifier column)
    regression <- cbind(rows_to_keep_TRS,PCA_analysis$x)
    
    fit <- lrm(test ~ PC1 + male_sex + Age_at_Interview, data = regression ,se.fit = T)
    fit2 <- glm(test ~ PC1 + male_sex + Age_at_Interview, data = regression , family = "binomial")
    
    tmp <- coef(summary(fit2))
    
    SE <- sqrt(diag(vcov(fit)))
    SE <- SE[2]
    Coef <- fit$coefficients[2]
    
    CI <- confint(fit2, level=0.95)
    CI <- CI[2,]
    # upper CI for beta
    CI_upper <- exp((Coef)+(1.96*SE))
    # lower CI for beta
    CI_lower <- exp((Coef)-(1.96*SE))
    tmp2 <- fit$stats
    Nagel_r2 <- tmp2[c("R2")]
    
    tmp3 <- c(input$DSM_2,input$Significance_threshold_2,tmp[2,], Nagel_r2, CI)
    
    
    
    regression_result <- "Still_need_to_complete"
    
    # Adjusted the Polygenic risk scores for analysis above, not sure the shiny app will still work
    Next_part_data <- list( Polygenic_risk_scores_for_analysis = Polygenic_risk_scores_for_analysis,
                            Columns_to_include = Columns_to_include,
                            rows_of_polygenic_risk_scores_to_use = rows_of_polygenic_risk_scores_to_use,
                            PCA_analysis = PCA_analysis
    )
  })
}