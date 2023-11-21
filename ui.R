#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with hccTAAb here: http://nscc.zzu.edu.cn/hccTAAb/

header <- dashboardHeader(titleWidth = 250,
                          ### changing logo
                          #tags$img(src = "ltd.png", height = "30px"),
                          title = shinyDashboardLogo(
                            theme = "blue_gradient",
                            boldText = tags$span(style = "font-family: Arial; font-size: 22px;", "hccTAAb"),
                            mainText = " Atlas",
                            badgeText = "v1.2"),
                          userOutput("user")
)

sidebar <- dashboardSidebar(
  width = 250,
  sidebarMenu(
    menuItem("Home", tabName = "home", icon = icon("house")),
    menuItem("TAAb Information", tabName = "infor", icon = icon("folder-open")),
    menuItem("Analysis",icon = icon("chart-column"),startExpanded = TRUE, 
             menuSubItem("Expression Pattern",tabName = "expression"),
             menuSubItem("Survival Analysis",tabName = "survival"),
             menuSubItem("Immune Infiltration",tabName = "immune"),
             menuSubItem("Similarity Analysis",tabName = "similarity"),
             menuSubItem("DNA Methylation",tabName = "methylation"),
             menuSubItem("DNA Mutation",tabName = "mutation")),
    menuItem("Datasets", tabName = "datasets", icon = icon("server"),startExpanded = FALSE,
             menuSubItem("Public Datasets",tabName = "public_datasets"),
             menuSubItem("Submit your Data",tabName = "submit")),
    menuItem("Help",icon = icon("hire-a-helper"),startExpanded = FALSE, 
             menuSubItem("Help",tabName = "help"),
             menuSubItem("About",tabName = "about"),
             menuSubItem("Update Logs",tabName = "update"))
  )
)

body <- dashboardBody(
  tabItems(
    ## home
    tabItem("home", 
            hr(),
            fluidPage(
              tags$div(style="margin: 0 auto; text-align: center; font-family: Kano; font-weight: bold; font-size: 55px; padding-bottom: 10px;",  
                       p(
                         strong("Welcome to hccTAAb Atlas", style="color:#3C8DBC;")
                       )),
              tags$div(style="margin: 0 auto; text-align: center; font-family: Times New Roman; font-size: 25px; padding-bottom: 15px;",
                       "hccTAAb Atlas is a comprehensive and accessible web-based resource that offers systematic knowledge on tumor-associated autoantibodies (TAAbs) 
                       in hepatocellular carcinoma (HCC). The hccTAAb Atlas provides detailed information about TAAbs and offers a range of functions, including expression 
                       analysis, survival analysis, similarity analysis, immune infiltration analysis, DNA methylation and mutation analysis. 
                       This tool serves as a valuable resource that addresses current gaps in our understanding of tumor-associated autoantibody functionality in HCC."), 
        
              ## carousel
              tags$div(#style = "margin-top: -20px;", # 调整上方外边距
                       carousel(id = "mycarousel",width = 12,
                                carouselItem(tags$img(src = "0_flow diagram.png",
                                                      style = "display: block; margin: 0 auto; width: 55%; height: 55%;")),
                                carouselItem(tags$img(src = "0_map.png",
                                                      style = "display: block; margin: 0 auto; width: 43%; height: 43%;")),
                                carouselItem(tags$img(src = "0_TAAb list.png",
                                                      style = "display: block; margin: 0 auto; width: 60%; height: 65%;"))
                                #carouselItem(tags$img(src = "3_immune.png",
                                                      #style = "display: block; margin: 0 auto; width: 52%; height: 52%;")),
                                #carouselItem(tags$img(src = "4_Survival_Correlation.png",
                                                      #style = "display: block; margin: 0 auto; width: 52%; height: 52%;")),
                                #carouselItem(tags$img(src = "5_DNA mutation_methylation.png",
                                                      #style = "display: block; margin: 0 auto; width: 55%; height: 55%;"))
                       ))
              
            )),
    
    ## TAAb information
    tabItem("infor",
            ## plotly figure
            box(title = tags$b("Compendium of Identified TAAbs for Detection of Hepatocellular Carcinoma (HCC)"), 
                width = 12, collapsible = T, solidHeader = T,
                status = "primary",
                plotlyOutput("plotly_info", height = "700px")),
            ## table
            fluidPage(width = 12, dataTableOutput("table_taabinfo"))
            #box(title = tags$b("Detailed Information Related to Each Tumor-Associated Autoantibody (TAAb)"), 
                #width = 12, collapsible = T, solidHeader = T,
                #status = "primary",
                #dataTableOutput("table_taabinfo"))
    ),
    
    ## expression
    tabItem("expression", 
            fluidPage(selectInput(inputId = "taab_name1", label = h4(strong("Select Interested Indicator")), 
                                  choices = idtaa, selected = "ATIC"),
                      fluidRow(column(2, actionButton("action1", label = "Calculator", icon = icon("calculator"))))
            ),
            br(),
            box(title = tags$b("Gene Expression Levels for the Selected TAAb"), 
                width = 12, collapsible = T, solidHeader = T, status = "primary",
                plotOutput("expr_gene", height = "36vh"),downloadButton("download_plot_expr_gene", 
                                                                    label = "Download", class = "text-primary")),
            #fluidPage(width = 12, plotOutput("expr_gene",height = 430)),
            br(),
            #fluidPage(width = 12, uiOutput("expr_protein")), ## UALCAN 
            box(title = tags$b("Protein Expression Levels for the Selected TAAb"), 
                width = 5, collapsible = T, solidHeader = T, status = "primary",
                plotOutput("expr_protein",height = "34vh"),downloadButton("download_plot_expr_protein", 
                                                                       label = "Download", class = "text-primary")),
            br(),
            box(title = tags$b("Table of Protein Expression for the Selected TAAb using HPA Pathology Data"), 
                width = 7, collapsible = T, solidHeader = T, status = "primary",
                dataTableOutput("expr_protein_hpa_pathology")),
            box(title = tags$b("Table of Protein Expression for the Selected TAAb using HPA Normal Data"), 
                width = 7, collapsible = T, solidHeader = T, status = "primary",
                dataTableOutput("expr_protein_hpa_normal"))
            #box(width = 12, plotOutput("expr_plot", height = 430, width = 1500))
            #box(title = "Exp_Plot",status = "success",solidHeader = TRUE, width = 12,plotOutput("expr_plot", height = 430, width = 1500))
    ),
    
    ## survival
    tabItem("survival", 
            fluidPage(
              selectInput(inputId = "taab_name2", label = h4(strong("Select Interested Indicator")), choices = idtaa, selected = "ATIC"),
              fluidRow(column(width = 2, actionButton("action2", label = "Calculator", icon = icon("calculator"))))),
            br(),
            box(title = tags$b("Survival Analysis (TCGA-LIHC)"), 
                width = 6, collapsible = T, solidHeader = T, status = "primary",
                fluidRow(
                  column(width = 6, plotOutput("survival_tcga",height = "460px"),
                         downloadButton("download_plot_sur_tcga", 
                                        label = "Download", class = "text-primary")),
                  column(width = 6, verbatimTextOutput("cox_summary1"))
                )
            ),
            
            box(title = tags$b("Survival Analysis (ICGC-LIRI)"), 
                width = 6, collapsible = T, solidHeader = T, status = "primary",
                fluidRow(
                  column(width = 6, plotOutput("survival_icgc",height = "460px"),
                         downloadButton("download_plot_sur_icgc", 
                                        label = "Download", class = "text-primary")),
                  column(width = 6, verbatimTextOutput("cox_summary2"))
                )
            )
    ),
    
    ## immune
    tabItem("immune", 
            fluidPage(
              fluidRow(column(width = 3, selectInput("taab_name3", label = h4(strong("Select Interested Indicator")), 
                                                     choices = idtaa, selected = "ATIC"))),
                       #column(width = 3, selectInput("immune_cells", label = h4(strong("Select Immune cells")), selected = "B_cells_naive",
                                                     #choices = unique(expr_Immune$immune_cell), multiple = TRUE))),
              fluidRow(column(width = 2, actionButton("action3", label = "Calculator", icon = icon("calculator"))))),
            br(),
            #fluidRow(column(width = 12, 
            #       box(title = tags$b("Barplot of ssGSEA Analysis for 29 Immunological Features accroding to the Selected TAAb Expression"), 
            #          width = 12, 
            #          collapsible = TRUE, solidHeader = TRUE, status = "primary",
            #         plotOutput("gsva_plot",height = 580)))),
            
            box(title = tags$b("Barplot of ssGSEA Analysis for 29 Immunological Features Accroding to the Selected TAAb Expression"), 
                width = 12, 
                collapsible = TRUE, solidHeader = TRUE, status = "primary",
                plotOutput("gsva_plot",height = 578),
                downloadButton("download_plot_immune_ssGSEA", label = "Download", class = "text-primary")),
            br(),
            #fluidRow(column(width = 12, plotOutput("immunecells_plot")))
            #fluidPage(plotOutput("immunecells_plot"))
            box(title = tags$b("Correlation Plot of 22 Immune Cell Scores Quantified by CIBERSORT with Selected TAAb Expression"), 
                width = 12, 
                collapsible = TRUE, solidHeader = TRUE, status = "primary",
                plotOutput("immunecells_plot",height = 1115),
                downloadButton("download_plot_immune_CIBERSORT", label = "Download", class = "text-primary"))
    ),
    
    ## Similarity
    tabItem("similarity", 
            fluidPage(
              fluidRow(column(width = 3, selectInput(inputId = "taab_name4", label = h4(strong("Select Interested Indicator")), choices = idtaa, selected = "ATIC")),
                       column(width = 3, sliderInput(inputId = "mumber_top",label = h4(strong("Select the Number to Display")), min = 1, max = 600, value = 200, step = 5))
              ),
              fluidRow(column(width = 2, actionButton("action4", label = "Calculator", icon = icon("calculator"))))),
            br(),
            box(title = tags$b("List of Positively Correlated Genes"), 
                width = 6, collapsible = T, solidHeader = T,
                status = "primary",
                dataTableOutput("positive_table")),
            box(title = tags$b("List of Negatively Correlated Genes"), 
                width = 6, collapsible = T, solidHeader = T,
                status = "primary",
                dataTableOutput("negative_table"))
    ),
    
    ## Methylation
    tabItem("methylation",
            fluidPage(width = 12, 
              box(title = tags$b("Table of DNA Methylation Differential Analysis Results"), 
                  status = "primary",width = 12,
                  dataTableOutput("methy_table"))
            )
    ),
    
    ## Mutation
    tabItem("mutation", 
            box(title = tags$b("Mutation Summary"), 
                width = 12, height = 8, collapsible = T, solidHeader = T,
                status = "primary",
                fluidPage(fluidRow(column(width = 4, plotOutput("mutation_plot1",height = 500)),
                                   column(width = 5, plotOutput("mutation_plot2",height = 500)),
                                   column(width = 3, plotOutput("mutation_plot3",height = 500))))),
            hr(),
            fluidPage(selectInput("taab_name5", label = h4(strong("Select Interested Indicator")), choices = idtaa, selected = "ATIC"),
                      fluidRow(column(width = 2, actionButton("action5", label = "Calculator", icon = icon("calculator"))))),
            hr(),
            box(title = tags$b("lollipopPlot of Selected TAAb for detection of HCC"), 
                width = 12, collapsible = T, solidHeader = T,
                status = "primary",
                plotOutput("mutation_plot")),
            hr(),
            fluidPage(width = 12, dataTableOutput("data_mutation"))
    ),
    
    ## dataset
    tabItem("public_datasets", 
            fluidPage(width = 12, 
                      HTML("<h3><strong>Public Datasets</strong></h3>"),
                      dataTableOutput("used_data"),
                      HTML("<h3><strong>Literature-based Data</strong></h3>
                           <p style='font-size: 22px; font-family:Times New Roman;'>
                           The data of the literature can be obtained by clicking on the TAAb name in the TAAb information table (first column).
                           </p>")
            )),
    tabItem("submit",
            fluidPage(width = 12, 
                      HTML("<h3><strong>Submit your Data</strong></h3>
                           <p style='font-size: 22px; font-family:Times New Roman;'>
                           hccTAAb Atlas allows for third-party data submission, which is crucial for keeping our data up-to-date. 
                           We strongly encourage users to upload yourself data by filling in the information below and submitting it. 
                           Your contribution to hccTAAb Atlas is greatly appreciated.
                           </p>")),
            fluidPage(width = 12,
              fluidRow(column(width = 3, selectInput(inputId = "data_types", label = h5(strong("Data Types")), choices = c("TAAb", "Public Data"))),
                       column(width = 3, textInput(inputId = "data_id",label = h5(strong("TAAb/Public Data")),
                                                   value = "Example: TP53/GSE144269")),
                       column(width = 6, textInput(inputId = "literature",label = h5(strong("Related Literature")),
                                                   value = "Example: A unique metastasis gene signature enables prediction of tumor relapse in early-stage hepatocellular carcinoma patients")),
                
              ),
              fluidRow(column(width = 3, textInput(inputId = "others", label = h5(strong("Other information")),
                                                   value = "Enter the other information")),
                       column(width = 3, textInput(inputId = "email",label = h5(strong("Please enter your email")),
                                                   value = "Example: xxxxxx@163.com"))),
              actionButton("submission_userdata",label = " Submit your Data",icon = icon("sign-in"),
                           style = "background-color: #3C8DBC; color: white;"),
              br(),
              br(),
              uiOutput("submission_message"),
              hr(),
            )
    ),
    
    ## About
    tabItem("about", 
            fluidPage(width = 12,
                      HTML("<h3><strong><span style='color: #3C8DBC'>Database Introduction</strong></h3>
                           <p style='font-size: 21px; font-family:Times New Roman;'>
                           hccTAAb Atlas is a comprehensive and accessible web-based resource that offers systematic knowledge on tumor-associated autoantibodies (TAAbs) 
                           in hepatocellular carcinoma (HCC) through a standardized processing pipeline. The hccTAAb Atlas provides detailed information about TAAbs 
                           and offers a range of functions, including expression analysis, survival analysis, correlation analysis, immune infiltration analysis, 
                           mutational landscape analysis, and DNA methylation analysis. This resource serves as a valuable platform for researchers and clinicians 
                           seeking to gain insights into the role of TAAbs in HCC.</p>
                           
                           <h3><strong><span style='color: #3C8DBC'>What are Tumor-Associated Autoantibodies ?</strong></h3>
                           <p style='font-size: 21px; font-family:Times New Roman;'>
                           Tumor-associated autoantibodies (TAAbs) are antibodies produced by the immune system in response to specific proteins expressed in tumor cells. 
                           They recognize and target proteins that are abnormal or modified in cancer cells. 
                           TAAbs can be detected in the blood and serve as biomarkers for cancer diagnosis, prognosis, 
                           and monitoring. They provide insights into tumor development and immune responses, offering potential for targeted therapies and personalized 
                           medicine in cancer treatment.
                           </p>

                           <h3><strong><span style='color: #3C8DBC'>Why Develop hccTAAb Atlas ?</strong></h3>
                           <p style='font-size: 21px; font-family:Times New Roman;'>
                           Numerous studies have shown the promising performance of tumor-associated autoantibodies (TAAbs) in hepatocellular carcinoma (HCC). 
                           Interestingly, during the transition from chronic liver disease to HCC, there is a stepwise increase in autoantibody levels, 
                           suggesting their significant role in HCC development. However, the current literature underreports the role of autoantibodies in immune regulation. 
                           To address these limitations, we developed hccTAAb Atlas, 
                           a comprehensive and user-friendly web-based database. This database provides a wide range of resources, 
                           including literature-based and multi-omics data, to facilitate autoantibody-related studies and promote in-depth basic 
                           and translational research in HCC.
                           </p>
                           
                           <h3><strong><span style='color: #3C8DBC'>How to Cite hccTAAb Atlas ?</strong></h3>
                           <p style='font-size: 21px; font-family:Times New Roman;'>TD Li, P Wang, GY Sun, YL Zou, YF Cheng, H Wang, Y Lu, JX Shi, KY Wang, Q Zhang, H Ye.    
                           hccTAAb Atlas: an integrated knowledge database for tumor-associated autoantibodies in hepatocellular carcinoma.</p>
                           
                           <h3><strong><span style='color: #3C8DBC'>Relevant Published Studies</strong></h3>
                           <ul>
                           <li><a href='https://pubmed.ncbi.nlm.nih.gov/32443439/' style='font-size: 21px;font-family:Times New Roman;'>Serological Biomarkers for Early Detection of Hepatocellular Carcinoma: A Focus on Autoantibodies against Tumor-Associated Antigens Encoded by Cancer Driver Genes</a></li>
                           <li><a href='https://pubmed.ncbi.nlm.nih.gov/34821436/' style='font-size: 21px;font-family:Times New Roman;'>A novel immunodiagnosis panel for hepatocellular carcinoma based on bioinformatics and the autoantibody‐antigen system</a></li>
                           <li><a href='https://pubmed.ncbi.nlm.nih.gov/36587394/' style='font-size: 21px;font-family:Times New Roman;'>Human Proteome Microarray identifies autoantibodies to tumor‐associated antigens as serological biomarkers for the diagnosis of hepatocellular carcinoma</a></li>
                           <li><a href='https://pubmed.ncbi.nlm.nih.gov/35052777/' style='font-size: 21px;font-family:Times New Roman;'>Autoantibody to GNAS in Early Detection of Hepatocellular Carcinoma: A Large-Scale Sample Study Combined with Verification in Serial Sera from HCC Patients</a></li>
                           </ul>
                           
                           <h3><strong><span style='color: #3C8DBC'>Contact</strong></h3>
                           <p style='font-size: 21px;font-family:Times New Roman;'> If you have any questions about the hccTAAb Atlas, please feel free to contact us: 
                           <a href='mailto:tiandonglee@163.com'>tiandonglee@163.com</a> (Tiandong Li);
                           <a href='mailto:tiandonglee@163.com'>yehua@zzu.edu.cn</a> (Hua Ye).
                           </p>
                           
                           <h3><strong><span style='color: #3C8DBC'>Acknowledgements</strong></h3>
                           <p style='font-size: 21px; font-family:Times New Roman;'> We express our sincere gratitude to Prof. HE Xu, GY Qin and MZ Fan for their invaluable technical support; 
                           Additionally, we extend our thanks to the <a href='http://nscc.zzu.edu.cn/' style='font-size: 20px;font-family:Times New Roman;'>National Supercomputing Center In Zhengzhou</a> for providing the platform.</p>
                           
                           <div style='margin-bottom: 20px;'></div>
                           "
                      ))),
    ## update logs
    tabItem("update", 
            fluidPage(width = 12,
                      HTML("<h3><strong><span style='color: #3C8DBC'>Update Logs</strong></h3>
                           <ul>
                           <li style='font-size: 21px;font-family:Times New Roman;'>2023.05.10: Web tool was developed, named hccTAAb Atlas (v1.0) </li>
                           <li style='font-size: 21px;font-family:Times New Roman;'>2023.08.15: Add the six aspects analysis: Expression Pattern, Survival Analysis, Immune Infiltration, 
                           Similarity Analysis, DNA Methylation, DNA Mutation (v1.1).</li>
                           <li style='font-size: 21px;font-family:Times New Roman;'>2023.10.06: Add the three proteomic data: PDC000198, PXD008373, PXD036048 (v1.2).</li>
                           <li style='font-size: 21px;font-family:Times New Roman;'>2023.10.20: Provide image download function (v1.2).</li>
                           </ul>
                           
                           <p style='font-size: 22px; color: #3C8DBC; font-family:Times New Roman;'><b>Note</b>: This update log provides a summary of recent changes and improvements made to the Shiny web page. 
                           Please feel free to explore the updated features and provide further feedback for ongoing enhancements.</p>
                           "
                           ))),
    
    ## help
    tabItem("help", 
            fluidPage(width = 12,
                      uiOutput("help"),
                      HTML("<h3><strong><span style='color: #3C8DBC'>Suggestions for Compatible Browsers</strong></h3>
                           <p style='font-size: 22px;font-family:Times New Roman;'>We recommend using the latest version of <b>Google Chrome</b>, <b>Mozilla Firefox</b>, and <b>Microsoft Edge</b>. </p>
                           </ul>
                           
                           <h3><strong><span style='color: #3C8DBC'>TAAb Information</strong></h3>
                           <img src='information.png' style='width: 1200px; height: auto;'>
                           <div style='margin-top: 1px;'></div>
                           <p style='font-size: 22px;font-family:Times New Roman;'>Click the <b>TAAb Information</b>, users can view the sensitivity and specificity of the 
                           corresponding TAAb (<b>One point represents one TAAb</b>), all the links in the table below are clickable. </p>
                           
                           <h3><strong><span style='color: #3C8DBC'>Expression Pattern</strong></h3>
                           <img src='expression.png'  style='width: 1200px; height: auto;'>
                           <div style='margin-top: 1px;'></div>
                           <p style='font-size: 22px;font-family:Times New Roman;'> 
                           Click on <b>Expression Pattern</b> and <b>Select a TAAb</b>. Then, click on the <b>Calculator</b>, users will be able to obtain figures and tables. 
                           By clicking on <b>Download</b>, users can save images locally.   | ns: <em>P</em> >0.05, *: <em>P</em><= 0.05, **: <em>P</em> <= 0.01, ***: <em>P</em> <= 0.001, ****: <em>P</em> <= 0.0001. </p>
                           
                           <h3><strong><span style='color: #3C8DBC'>Survival Analysis</strong></h3>
                           <img src='survival.png'  style='width: 1200px; height: auto;'>
                           <div style='margin-top: 1px;'></div>
                           <p style='font-size: 22px;font-family:Times New Roman;'> 
                           Click on <b>Survival Analysis</b> and <b>Select a TAAb</b>. Then, click on the <b>Calculator</b>, users will be able to obtain Kaplan-Meier plots and hazard ratios. 
                           By clicking on <b>Download</b>, users can save images locally. </p>
                           
                           
                           <h3><strong><span style='color: #3C8DBC'>Immune Infiltration</strong></h3>
                           <img src='immune.png' style='width: 1200px; height: auto;'>
                           <div style='margin-top: 1px;'></div>
                           <p style='font-size: 22px;font-family:Times New Roman;'> 
                           Click on <b>Immune Infiltration</b> and <b>Select a TAAb</b>. Then, click on the <b>Calculator</b>, users will be able to obtain 
                           the relationship between selected TAAb and immune infiltration levels (ssGSEA and CIBERSORT). By clicking on <b>Download</b>, users can save images locally. </p>
                           
                           <h3><strong><span style='color: #3C8DBC'>Similarity Analysis</strong></h3>
                           <img src='similarity.png'  style='width: 1200px; height: auto;'>
                           <div style='margin-top: 1px;'></div>
                           <p style='font-size: 22px;font-family:Times New Roman;'> 
                           Click on <b>Similarity Analysis</b> and <b>Select a TAAb</b>. Then, click on the <b>Calculator</b>, users will be able to obtain genes 
                           that exhibit both positive and negative correlations with the selected TAAb. Users can change the number of genes by clicking the <b>Select the Number to Display</b> (the default value is 100).</p>
                           
                           <h3><strong><span style='color: #3C8DBC'>DNA Methylation</strong></h3>
                           <img src='methylation.png' alt='Expression Pattern'  style='width: 1200px; height: auto;'>
                           <div style='margin-top: 1px;'></div>
                           <p style='font-size: 22px;font-family:Times New Roman;'> 
                           Click on <b>DNA Methylation</b>, users will be able to obtain the table of DNA methylation differential analysis.</b> 
                    
                           <h3><strong><span style='color: #3C8DBC'>DNA Mutation</strong></h3>
                           <img src='mutation.png' alt='Expression Pattern'  style='width: 1200px; height: auto;'>
                           <div style='margin-top: 1px;'></div>
                           <p style='font-size: 22px;font-family:Times New Roman;'> 
                           Click on <b>DNA Mutation</b> and <b>Select a TAAb</b>. Then, click on the <b>Calculator</b>, users will be able to obtain 
                           lollipop plot and detailed information table for the selected TAAb.</b>
                           
                           <div style='margin-bottom: 30px;'></div>
                           "
                      )))
    
  )
)

ui <- dashboardPage(header=header, 
                    sidebar=sidebar, 
                    body=body, 
                    controlbar = dashboardControlbar(skinSelector()),
                    title = "hccTAAb Atlas")

