#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


# Define UI for application that draws a histogram
shinyUI( navbarPage("Polygenic Risk Score Analysis Viewer Environment",
            fluid = T,
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
              uiOutput("Significance_threshold"
                       ),
              checkboxGroupInput("Gene_regions", label = "Length of Gene regions:",
                                 choices = c("extended","normal","Full")),
              uiOutput("DSM"),
              uiOutput("geneset"),
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
                          plotOutput('PvalPlot'),
                         plotOutput('Beta_plot'),
                         plotOutput('R2_plot')),
                        conditionalPanel(condition = "input.tabs == 'Plots' & input.step_regress == 'Yes'",
                                         plotOutput('PvalPlot_step'),
                                         plotOutput('Beta_plot_step'),
                                         plotOutput('R2_plot_step'))),
                tabPanel("Table", dataTableOutput("summary_table")),
                tabPanel("Input variables",
                         rclipboardSetup(),
                         

                         textAreaInput("text_2",label = "Full message will appear here:",width = "500px",height = "100px", resize = "both",
                                       placeholder = "Twitter handles will appear here at the end of your message depending on the options selected to the left (eg: Pint of Science is Great! @virustinkerer)")
                         ,
                         # UI ouputs for the copy-to-clipboard buttons
                         uiOutput("clip"))
              )
            )
          )
          ),
          tabPanel("Combined_PRS",
          
            sidebarLayout(         
              ## Write shiny UI across 4 parameters in the data table
              sidebarPanel(
                fileInput("file2", "Choose regression input file",
                      multiple = F),
                fileInput("file3", "Choose Polygenic risk score file",
                          multiple = F),
                uiOutput("GWAS_to_include"),
                uiOutput("Significance_threshold_2"
                ),
                uiOutput("DSM_2"),
              
              sliderInput(inputId = "alpha", 
                          label = "alpha of PCA plot", 
                          min = 0, 
                          max = 1,
                          value = 0.2),
              sliderInput(inputId = "varname.size", 
                          label = "Size of variable name", 
                          min = 0, 
                          max = 10,
                          value = 3),
              sliderInput(inputId = "varname.adjust", 
                          label = "distance of variable name from factor loading", 
                          min = 0, 
                          max = 10,
                          value = 2)
              ),
            mainPanel(
              tabsetPanel(id = "tabs_2",
                          
                          tabPanel("Plots_2",
                                     plotOutput('PCA_plot'),
                                     plotOutput('R2_plot_2'),
                                   verbatimTextOutput('Significance_threshold_used'),
                                     textOutput('Number_of_individuals_used')),
                                    
  
                          tabPanel("Table_2", dataTableOutput("summary_table_2"))
                          
            )
          )
          )
        )
) 
)
