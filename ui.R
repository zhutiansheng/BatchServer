sidebar <- dashboardSidebar(
  sidebarMenu(id="mysidebar",
    menuItem("Home", tabName = "home", icon = icon("home")),
    menuItem("Data Input", icon = icon("th"), tabName = "dataInput",
             badgeLabel = "matrix", badgeColor = "green"),
    menuItem("Evaluatation", icon = icon("angle-right"), tabName = "evaluatation",
             badgeColor = "yellow",
             menuSubItem("PVCA", tabName = "pvca", icon = icon("angle-left")),
             menuSubItem("UMAP", tabName = "umap", icon = icon("angle-left"))
             ),
    menuItem("Correction", icon = icon("line-chart"), tabName = "elimination",
             badgeColor = "blue",
             menuSubItem("ComBat", tabName = "combat", icon = icon("angle-left"))
             #,menuSubItem("RandomForest", tabName = "rf", icon = icon("angle-left"))
             ),
    menuItem("Readme", tabName = "readme",badgeLabel = "help", badgeColor = "red", icon=icon("mortar-board")),
    menuItem("About", tabName = "about", icon = icon("question"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "home",
            #h2("Welcome to Batch Server home"),
            includeMarkdown("help/main.Rmd")
            #HTML('<p class="MsoNormal">Batch effects are unwanted sources of variation irrelevant to biological variation inevitably introduced to the samples during experimental handling which would obscure the biological signal. Batch effects are one of the biggest challenges faced by high throughput omics science, especially in the context of large cohort of thousands of samples. Existing batch effect-correcting tools focus mainly on the development of methods that are not convenience of use, usually requiring extensive coding experiences, sometimes even need to know the prior distribution of the data. Moreover, few tools offer both evaluation and correction of batch effects. We developed an open-source web server-based batch effect correction tool, namely BatchServer, which enables users to interactively evaluate and correct batch effects of a variety of omics data.')
    ),
    
    tabItem(tabName = "dataInput",
            # Input: Select a file ----
            fileInput("myd", "Please choose your data file",
                      multiple = FALSE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv",
          "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
          ".xlsx")),
            # Input: Checkbox if file has header ----
            #checkboxInput("header", "Header", TRUE),
            
            # Input: Select separator ----
            radioButtons("sep", "Separator",
                         choices = c(Comma = ",",
                                     Semicolon = ";",
                                     Tab = "\t",
                                     'xls/xlsx' = "xlsx"),
                         selected = ",",inline = T),
          radioButtons("missing_replace_input",
                       "Missing value replacement",
                       choices = c(
                         "None" = "none",
                         "1" = '1',
                         '0' = "0",
                         "10% of minimum" = '0.1',
                         "minimum" = "minimum"
                       ),inline = TRUE,selected = "0"),  
          checkboxInput("qn", "Quantile normalization", TRUE),
          checkboxInput("log2", "Log2 transform", TRUE), 
            # Horizontal line ----
            tags$hr(),
            # Input: Select a file ----
            fileInput("sample_info", "Please choose your sample information file",
                      multiple = FALSE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv",
          "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
          ".xlsx")),
            # Input: Checkbox if file has header ----
            #checkboxInput("sample_header", "Header", TRUE),
            
            # Input: Select separator ----
            radioButtons("sample_sep", "Separator",
                         choices = c(Comma = ",",
                                     Semicolon = ";",
                                     Tab = "\t",
                                     'xls/xlsx' = "xlsx"),
                         selected = ",",inline = T),
              
            # Horizontal line ----
            tags$hr(),
            actionButton("input_submit", "Submit", class = "btn-primary"),
            verbatimTextOutput("upload_note")
    ),
    
    tabItem(tabName = "pvca",
            h3("PVCA"),
            h4('PVCA assess the batch sources by fitting all "sources" as random effects including two-way interaction terms in the Mixed Model(depends on lme4 package) to selected principal components, which were obtained from the original data correlation matrix. Pierre Bushel (2019). pvca: Principal Variance Component Analysis (PVCA). R package version 1.24.0.'),
            selectInput("pvca_effect_name","Select contributing effect column name(s)",
                        choices = effect_name,multiple = T),
            sliderInput("pvca_threshold", "Set the percentile value of the minimum amount of the variabilities that the selected principal components need to explain",
                        min = 0, max = 1,
                        value = 0.7, step = 0.1),
            actionButton("pvca_submit", "Submit", class = "btn-primary"),
            tags$hr(),
            tabsetPanel(
            tabPanel(
              "Pieplot",
              plotOutput("draw_pie"),
              uiOutput("pvca_pie_ui")
              ),
            tabPanel(
              "Barplot",
              column(12,plotOutput("draw_pvca"),
                     uiOutput("pvca_ui"))
            )
            )
    ),
    tabItem(tabName = "umap",            
            h3("UMAP"),
            h4("Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique that can be used for visualisation similarly to t-SNE, but also for general non-linear dimension reduction. Please note that missing value is not allowed in data matrix in umap."),      
            hr(),
            radioButtons("missing_replace",
              "Missing value replacement",
              choices = c(
                "None" = "none",
                "1" = '1',
                '0' = "0",
                "10% of minimum" = '0.1',
                "minimum" = "minimum"
              ),inline = TRUE,selected = "0"),
	    sliderInput("n_neighbors", "number of nearest neighbors:",
                        min = 1, max = 100,
                        value = 15),
            radioButtons("metric", "distances method",
                         choices = c(euclidean = "euclidean",
                                     manhattan = "manhattan",
                                     cosine = "cosine",
                                     pearson = "pearson"),
                         selected = "euclidean",inline = T),
            sliderInput("n_epochs", "number of  iterations:",
                        min = 1, max = 1000,
                        value = 200),
            radioButtons("init", "initial coordinates",
                         choices = c(spectral = "spectral",
                                     random = "random"),
                         selected = "spectral",inline = T),
            sliderInput("min_dist", "minimumn dist in the final layout",
                        min = 0.01, max = 1,
                        value = 0.1, step = 0.01),
            numericInput("alpha", "initial value of 'learning rate'", 1,
                         0.1, 10, 0.1),
            numericInput("gamma", "learning rate", 1,
                         0.1, 10, 0.1),
            numericInput("negative_sample_rate", "non-neighbor points are used per point and
                         per iteration during layout optimization", 5,
                         1, 100, 1),
            
            actionButton("umap_submit", "Calculate", class = "btn-primary"),
            verbatimTextOutput("umap_note"),
            tags$hr(),
            selectInput("umap_effect_name","Select contributing effect column name(s)",
                        choices = effect_name,multiple = F),
            plotlyOutput("draw_umap"),
            uiOutput("umap_ui")
    ),
    tabItem(tabName = "combat",
            h3("ComBat"),
            h4("The ComBat function adjusts for known batches using an empirical Bayesian
framework. So known batch variable is required in your dataset. Here you should pay attention to the [parametric estimate method] choice, which was improved compare to the original ComBat method. The option [automatic] will automatically decide to set parametric estimate method to parametric or nonparametric according to the data distribution."),
            hr(),
            selectInput("batch_effect_name","Select known batch effect column name",
                        choices = NULL,multiple = F),
            
            tags$div(title="Surrogate variables are covariates constructed directly from high-dimensional data (like gene expression/RNA sequencing/methylation/brain imaging data) that can be used in subsequent analyses to adjust for unknown, unmodeled, or latent sources of noise.",
            selectInput("adjust_variables","Select surrogate variable(s)",
                        choices = NULL,multiple = T)
            ),
            radioButtons("par.prior", "Parametric estimate method",
                         choices = c(automatic= "auto",
                                     parameter = "parameter",
                                     noparameter = "noparameter"),
                         selected = "auto",inline = T),
            # radioButtons("fit.method", "Fitness method",
            #              choices = c("maximum likelihood" = "mle",
            #                          "moment matching" =  "mme",
            #                          "quantile matching" = "qme",
            #                          "maximizing goodness-of-fit estimation" = "mge"),
            #              selected = "mle",inline = T),
            radioButtons("mean.only", "Only adjusts the
mean of the batch effects across batches (default adjusts the mean and variance)",
                         choices = c(No = FALSE,
                                     Yes = TRUE),
                         selected = FALSE,inline = T),
            actionButton("elimination_submit", "Submit", class = "btn-primary"),
            verbatimTextOutput("combat_log"),
            uiOutput("combat_ui")
    ),
    # tabItem(tabName = "rf",
    #         h3("Description:"),
    #         h4("Remove most importances batch related variables using Random Forest "),
    #         hr(),
    #         selectInput("batch_effect_name_rf","Select Known Batch Effect Column Name",
    #                     choices = NULL,multiple = F),
    #         numericInput("ntree", "Number of trees to grow", 500,
    #                      1, 5000, 1),
    #         numericInput("nodesize", "Minimum size of terminal nodes", 5,
    #                      1, 500, 1),
    #         numericInput("topN", "Number of top effect batch related variables to delete", 5,
    #                      1, 500, 1),
    #         actionButton("rf_submit", "Submit", class = "btn-primary"),
    #         verbatimTextOutput("rf_log"),
    #         uiOutput("rf_ui")
    # ),
    tabItem(tabName = "readme",
            h1("Readme"),
            h4("Batch effects are unwanted data variations that may obscure biological signals, leading to bias or errors in subsequent data analyses. Effective evaluation and elimination of batch effects are necessary for omics data analysis. In order to facilitate the evaluation and correction of batch effects, here we present BatchSever, an open-source R/Shiny based user-friendly interactive graphical web platform for batch effects analysis. As the autohrs of original ComBat have extensively investigated, if the experiment has not been properly designed, or if the batch design information is missing, no effective batch correction could be performed. Unbalanced batch-group design and inappropriate missing value imputation will pose challenges to effective batch effect correction."),
            h2("Test data download"),
            downloadButton("testData_download", "dataMatrix", class = "btn-primary"),
            downloadButton("sampleData_download", "sampleInfo", class = "btn-primary"),

            includeMarkdown("help/readme.Rmd")

            ),
    tabItem(tabName = "about",
            h3("Software author:"),
            HTML("Tiansheng Zhu; tszhu @ fudan.edu.cn"),
            h3("License:"),            
            HTML("BatchServer is an open-source software implemented in pure R language and the source code is freely available at https://github.com/zhutiansheng/BatchServer. 
Now Batch Server is supported by both school of computer science of Fudan University (zhou's lab: admis.fudan.edu.cn) and school of life sciences of Westlake University (guo's lab: www.guomics.com). The software is published by ''")         
    )
  )
)
dashboardPage(skin = "purple",
                    dashboardHeader(title = "BatchServer"),
                    sidebar,
                    body
)
