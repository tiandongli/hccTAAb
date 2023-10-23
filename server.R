#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with hccTAAb here: http://nscc.zzu.edu.cn/hccTAAb/

server <- function(input, output, session) {
  ## user information
  output$user <- renderUser({
    dashboardUser(
      name = "tiandonglee@163.com", 
      image = "ltd.png", 
      title = "Tiandong Li (黎天东)",
      subtitle = "", 
      footer = p("Biostatistics | Bioinformatics | Tumor Immunology | Cancer Epidemiology", class = "text-center"),
      fluidRow(
        column(12,
               tags$div(class = "text-center",
                        socialButton(
                          href = 'mailto:tiandonglee@163.com',
                          icon = icon("envelope")
                        ),
                        socialButton(
                          href = 'mailto:litiandong@outlook.com',
                          icon = icon("envelope-open-text")
                        ),
                        socialButton(
                          href = "https://github.com/tiandongli",
                          icon = icon("github",style = "color:white")
                        )
               )
        )
      )
    )
  })
  
  ## 2.infor
  ## 2.1 plotly
  output$plotly_info <- renderPlotly({
    plot_ly(data = TAAb_infor, x = ~`100-Specificity`, y = ~Sensitivity, symbol =  ~TAAb,
            hoverinfo = "text",text = paste("TAAb:",TAAb_infor$TAAb,"<br>","Sensitivity:",
                                            TAAb_infor$Sensitivity,"%","<br>","Specificity:",
                                            TAAb_infor$Specificity,"%"))  %>% 
      add_markers(marker = list(size = 16)) %>%
      layout(xaxis = list(range = c(0,100), title = "100-Specificity(%)", titlefont = list(size = 24)),
             yaxis = list(range = c(0,100), title = "Sensitivity(%)", titlefont = list(size = 24)))
  })
  
  ## 2.2 data infor
  ## data building
  output$table_taabinfo <- DT::renderDataTable(
    TAAb_infor_web, escape = FALSE, 
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 20, 
                   scrollX = TRUE,
                   scrollCollapse = TRUE,
                   fixedColumns = TRUE,
                   columnDefs = list(list(className = 'dt-center', targets = '_all')),
                   scrollY = '60vh')
    )
  
  ## 3.analysis
  ## 3.1 expression
  observeEvent(input$action1, {
    taab_name <- isolate(input$taab_name1)
    ## gene expression data
    expr_gene_target <- merged_data[merged_data$Gene==taab_name,] #target expression matrix
    ## protein expression data
    expr_protein_target <- expr_protein[expr_protein$Symbol==taab_name,] #target expression matrix
    hpa_pathology_lctaa_filter <- hpa_pathology_lctaa[hpa_pathology_lctaa$`Gene name`== taab_name,] #HPA-pathology
    hpa_normal_lctaa_filter <- hpa_normal_lctaa[hpa_normal_lctaa$`Gene name`== taab_name,] #HPA-normal
    
    
    ## 1)gene expression boxplot
    ## define plot
    if(length(unique(expr_gene_target$Dataset))==1){
      expr_gene <- ggplot(expr_gene_target,aes(x = Group, y = Expression, fill = Group)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(size=1.5, alpha=0.8, position=position_jitter(0.25))+
          #geom_signif(comparisons = list(c("Tumor", "Non-tumor")),map_signif_level = T) + 
          stat_compare_means(label = "p.signif", method = "wilcox.test", label.x = 1.5, col = "red") +
          ylab(paste0(taab_name)) +
          xlab("") +
          theme_few() +
          facet_wrap(~Dataset,scales = "free_y",nrow = 1) + #每个面板使用独立的y轴刻度
          scale_fill_nejm() +
          scale_color_nejm() +
          theme(plot.title = element_text(size = 23, face = "bold"),
                legend.position = "right",
                text = element_text(size = 20),
                strip.text = element_text(size = 18),
                axis.text.x = element_blank())
    } else {
      expr_gene <- ggplot(expr_gene_target,aes(x = Group, y = Expression, fill = Group)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(size=1.5, alpha=0.8, position=position_jitter(0.25))+
          #geom_signif(comparisons = list(c("Tumor", "Non-tumor")),map_signif_level = T) + 
          stat_compare_means(label = "p.signif", method = "wilcox.test", label.x = 1.5, col = "red") +
          ylab(paste0(taab_name)) +
          xlab("") +
          theme_few() +
          facet_wrap(~Dataset,scales = "free_y",nrow = 1) + #每个面板使用独立的y轴刻度
          scale_fill_nejm() +
          scale_color_nejm() +
          theme(plot.title = element_text(size = 23, face = "bold"),
                legend.position = "right",
                text = element_text(size = 20),
                strip.text = element_text(size = 18),
                axis.text.x = element_blank())}
    ## Display plot
    if(length(unique(expr_gene_target$Dataset))==1){
      output$expr_gene <- renderPlot({expr_gene},
        width = 300
      )
    } else {
      output$expr_gene <- renderPlot({expr_gene},
        width = 200*length(unique(expr_gene_target$Dataset))
      )}
    ## download plot
    output$download_plot_expr_gene <- downloadHandler(
      filename = function() {
        paste0(taab_name,"_Gene_Expression_hccTAAb Atlas_",Sys.Date(),".png")
      },
      content = function(file) {
        ggsave(file,  plot = expr_gene, width = 60*length(unique(expr_gene_target$Dataset)), height = 100, dpi = 600,units = "mm", device = "png")
      }
    )
    
    ## 2)protein table (UALCAN) 
    #output$expr_protein <- renderUI({
      #url <- paste0("https://ualcan.path.uab.edu/cgi-bin/CPTAC-Result.pl?genenam=", taab_name, "&ctype=Liver")
      #tags$a(href = url, target = "_blank", HTML(paste0("<span style='font-size: 22px; font-weight: bold;'>Click to view the protein expression results of ", taab_name, "</span>")))
    #})
    
    ## 2) protein expression
    expr_protein <- ggplot(expr_protein_target,aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size=1.5, alpha=0.8, position=position_jitter(0.25))+
      #geom_signif(comparisons = list(c("Tumor", "Non-tumor")),map_signif_level = T) + 
      stat_compare_means(label = "p.signif", method = "wilcox.test", label.x = 1.5, col = "red") +
      ylab(paste0(taab_name)) +
      xlab("") +
      theme_few() +
      facet_wrap(~Datasets,scales = "free_y",nrow = 1) + #每个面板使用独立的y轴刻度
      scale_fill_jama() +
      scale_color_jama() +
      theme(plot.title = element_text(size = 23, face = "bold"),
            legend.position = "right",
            text = element_text(size = 20),
            strip.text = element_text(size = 18),
            axis.text.x = element_blank())
    
    ## display   
    if (taab_name %in% expr_protein_target$Symbol) {
      output$expr_protein <- renderPlot({
        expr_protein},
        width = 228*length(unique(expr_protein_target$Datasets)))
    } else {
      output$expr_protein <- renderPlot({
        plot(1, 1, type = "n", xlab = "", ylab = "", main = "")
        text(1, 1, paste0(taab_name," is not available in protein datasets"),cex=2, col="red")
      })
    }
     
    ## download
    output$download_plot_expr_protein <- downloadHandler(
      filename = function() {
        paste0(taab_name,"_Protein_Expression_hccTAAb Atlas_",Sys.Date(),".png")
      },
      content = function(file) {
        ggsave(file,  plot = expr_protein, width = 70*length(unique(expr_protein_target$Datasets)), height = 100, dpi = 600,units = "mm", device = "png")
      }
    )
    
    output$expr_protein_hpa_pathology <- DT::renderDataTable(
      hpa_pathology_lctaa_filter, escape = FALSE, 
      selection = "none",
      rownames = FALSE,
      options = list(pageLength = 20, 
                     scrollX = TRUE, 
                     lengthChange = FALSE, 
                     columnDefs = list(list(className = 'dt-center', targets = '_all')),
                     searching = FALSE))
    
    output$expr_protein_hpa_normal <- DT::renderDataTable(
      hpa_normal_lctaa_filter, escape = FALSE, 
      selection = "none",
      rownames = FALSE,
      options = list(pageLength = 20, 
                     scrollX = TRUE, 
                     lengthChange = FALSE, 
                     columnDefs = list(list(className = 'dt-center', targets = '_all')),
                     searching = FALSE))
  })
  
  ## 3.2 survival   
  observeEvent(input$action2, {
    taab_name <- isolate(input$taab_name2)
    
    ## 1)TCGA survival data
    survi_cutoff1 <- surv_cutpoint(data = expr_tcga_os, time = "os", event = "Event",
                                   variables = taab_name)
    sur_cat1 <- surv_categorize(survi_cutoff1)
    names(sur_cat1)[3] <- "Group"
    fit1 <- survfit(Surv(os, Event) ~Group, data = sur_cat1)
    
    output$cox_summary1 <- renderPrint({
      cox_model1 <- coxph(Surv(os, Event) ~ get(taab_name), data = expr_tcga_os)
      summary_text1 <- capture.output(summary(cox_model1))
      summary_text1 <- summary_text1[-c(1:5)]  # 去除不需要的行
      summary_text1 <- gsub("^\\s+", "", summary_text1) # 去除每行开头的空格
      summary_text1 <- gsub("get\\(taab_name\\)", taab_name, summary_text1) # 替换为对应的 taab_name
      title_text1 <- "##########################################################\n        Cox Regression Analysis(using TCGA data)\n##########################################################\n"
      cat(paste(title_text1, paste(summary_text1, collapse = "\n"), sep = "\n"), "\n")
      cat(paste("\n##########################################################\n",
                "            ","HR=", format(round(summary(cox_model1)$coef[, "exp(coef)"], 2), nsmall = 2), 
                "(95% CI: ", format(round(summary(cox_model1)$conf.int[,"lower .95"], 2), nsmall = 2), 
                "-", format(round(summary(cox_model1)$conf.int[,"upper .95"], 2), nsmall = 2),")", sep = ""), "\n##########################################################") 
    })

    ## survival_plot1
    survival_tcga <- ggsurvplot(fit1, data = sur_cat1, 
                                conf.int = T,  # 95%CI
                                size = 1,  # change line size
                                pval = T,  # p-value of log-rank test
                                risk.table = TRUE,  
                                palette = 'jco', #"npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
                                xlab = c("Months"),
                                ggtheme = theme_survminer()+ #theme_minimal() theme_gray() theme_classic() theme_bw()
                                  theme(axis.text = element_text(size = 18),  # 调整坐标轴标签字号
                                        axis.title = element_text(size = 20),  # 调整坐标轴标题字号
                                        legend.text = element_text(size = 16),  # 调整图例文本字号
                                        legend.title = element_text(size = 16)  # 调整图例标题字号
                                   ),
                                 legend.title = paste0(taab_name), 
                                 legend.labs = c("High","Low")) 
      
    
    ## display
    output$survival_tcga <- renderPlot({
      survival_tcga},
      height = 465, width = 420)
      
    ## download
    output$download_plot_sur_tcga <- downloadHandler(
      filename = function() {
        paste0(taab_name,"_Survival_TCGA_hccTAAb Atlas_",Sys.Date(),".png")
      },
      content = function(file) {
        png(file, width = 4.5, height = 5.3, units = "in", res = 600)
        print(survival_tcga)
        dev.off()
      }
    )
    
    ## 2)ICGC survival data
    if (taab_name %in% colnames(expr_icgc_os)) {
      survi_cutoff2 <- surv_cutpoint(data = expr_icgc_os, time = "os", event = "Event",
                                     variables = taab_name)
      sur_cat2 <- surv_categorize(survi_cutoff2)
      names(sur_cat2)[3] <- "Group"
      fit2 <- survfit(Surv(os, Event) ~Group, data = sur_cat2)
      
      output$cox_summary2 <- renderPrint({
        cox_model2 <- coxph(Surv(os, Event) ~ get(taab_name), data = expr_icgc_os)
        summary_text2 <- capture.output(summary(cox_model2))
        summary_text2 <- summary_text2[-c(1:5)]  # 去除不需要的行
        summary_text2 <- gsub("^\\s+", "", summary_text2) # 去除每行开头的空格
        summary_text2 <- gsub("get\\(taab_name\\)", taab_name, summary_text2) # 替换为对应的 taab_name
        title_text2 <- "##########################################################\n        Cox Regression Analysis(using ICGC data)\n##########################################################\n"
        cat(paste(title_text2, paste(summary_text2, collapse = "\n"), sep = "\n")) #输出最终的摘要信息
        cat(paste("\n##########################################################\n",
                  "            ","HR=", format(round(summary(cox_model2)$coef[, "exp(coef)"], 2), nsmall = 2), 
                  "(95% CI: ", format(round(summary(cox_model2)$conf.int[,"lower .95"], 2), nsmall = 2), 
                  "-", format(round(summary(cox_model2)$conf.int[,"upper .95"], 2), nsmall = 2),")", sep = ""), "\n##########################################################") 
      })
      
      ## survival_plot2
      survival_icgc <- ggsurvplot(fit2, data = sur_cat2, 
                                  conf.int = T,  # 95%CI
                                  size = 1,  # change line size
                                  pval = T,  # p-value of log-rank test
                                  risk.table = TRUE,  
                                  palette = 'jco', #"npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty"
                                  xlab = c("Months"),
                                  ggtheme = theme_survminer()+ #theme_minimal() theme_gray() theme_classic() theme_bw()
                                    theme(axis.text = element_text(size = 18),  # 调整坐标轴标签字号
                                          axis.title = element_text(size = 20),  # 调整坐标轴标题字号
                                          legend.text = element_text(size = 16),  # 调整图例文本字号
                                          legend.title = element_text(size = 16)  # 调整图例标题字号
                                    ),
                                  legend.title = paste0(taab_name), 
                                  legend.labs = c("High","Low"))
    ## display
    output$survival_icgc <- renderPlot({
      survival_icgc},
      height = 465,width = 420)
  } else {
    output$survival_icgc <- renderPlot({
      plot(1, 1, type = "n", xlab = "", ylab = "", main = "")
      text(1, 1, paste0(taab_name," is not available in ICGC-LIRI dataset"),cex=2, col="red")
    }, height = 465,width = 420)
  }
    
    ## download
    output$download_plot_sur_icgc <- downloadHandler(
      filename = function() {
        paste0(taab_name,"_Survival_ICGC_hccTAAb Atlas_",Sys.Date(),".png")
      },
      content = function(file) {
        png(file, width = 4.5, height = 5.3, units = "in", res = 600)
        print(survival_icgc)
        dev.off()
      }
    )  
  })
  
  ## 3.3 immune 
  observeEvent(input$action3, {
    gene <- isolate(input$taab_name3)
    ## CIBERSORT
    expr_Immune_target <- as.data.frame(expr_Immune)
    ## GSVA data
    gsva_immune_target <- gsva_taab[, c(1:29, which(names(gsva_taab) %in% gene))]
    gsva_immune_target <- pivot_longer(gsva_immune_target,
                                       cols = 1:29,
                                       names_to = "Immune_Pathway",
                                       values_to = "Immune_Score") %>% as.data.frame()
    colnames(gsva_immune_target)[1] <- "Group"
    gsva_immune_target$Group <- factor(gsva_immune_target$Group, levels=c('High','Low'))
    
    ## 1)ssGSEA
    gsva_plot <- ggplot(gsva_immune_target,aes(x=Immune_Pathway, y=Immune_Score, fill=Group))+
      #geom_violin(trim = TRUE,adjust = 1)+
      #geom_boxplot(outlier.shape = NA,color="white",width=0.2)+
      geom_boxplot(aes(fill=Group),  outlier.shape = NA)+
      geom_jitter(size=0.8, alpha=0.4, position=position_jitter(0.3))+ 
      #geom_signif(comparisons = list(c("High", "Low")),map_signif_level = T)+
      stat_compare_means(label = "p.signif", method = "wilcox.test", label.x = 1.5, 
                         label.y.npc = "top", hide.ns = TRUE, col="red")+
      theme_pubr() +
      scale_fill_nejm() +
      scale_color_nejm() +
      scale_x_discrete(name="")+ 
      labs(y="Immune score") +
      labs(fill=paste0(gene," Expression"))+
      scale_y_continuous(limits = c(0, NA))+ 
      theme(plot.title = element_text(size = 24, face = "bold"),
            legend.position = "right",
            text = element_text(size = 18),
            #axis.text.x = element_blank(), ## not x_text
            axis.title = element_text(face="bold"),
            axis.text.x=element_text(size = 16, angle = 45, hjust = 1.0, vjust = 1.0))
    
    output$gsva_plot <- renderPlot({
      gsva_plot},
      height = 570,width = 1700)
    
    output$download_plot_immune_ssGSEA <- downloadHandler(
      filename = function() {
        paste0(gene,"_immune_ssGSEA_hccTAAb Atlas_",Sys.Date(),".png")
      },
      content = function(file) {
        ggsave(file,  plot = gsva_plot, width = 450, height = 200, dpi = 600,units = "mm", device = "png")
      }
    )
    
    ## 2)CIBERSORT
    immunecells_plot <- ggplot(expr_Immune_target, aes(x=expr_Immune_target[,gene], y=infiltration)) + 
      geom_smooth(method = 'lm', se = T, color = '#1661B8')+
      stat_cor(data=expr_Immune_target, method = "spearman",color = "#E41A1C",na.rm = T)+
      geom_point() + 
      theme_few() +
      labs(x=paste0(gene)) +
      facet_wrap(~immune_cell,scales = "free_y", ncol = 6)+
      theme(text = element_text(size = 22),
            strip.text = element_text(size = 14),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 14))
 
    output$immunecells_plot <- renderPlot({
      immunecells_plot}, width = 1680
    )
    
    output$download_plot_immune_CIBERSORT <- downloadHandler(
      filename = function() {
        paste0(gene,"_immune_CIBERSORT_hccTAAb Atlas_",Sys.Date(),".png")
      },
      content = function(file) {
        ggsave(file,  plot = immunecells_plot, width = 600, height = 400, dpi = 600,units = "mm", device = "png")
      }
    )
  })
  ## 3.4 similergene
  observeEvent(input$action4, {
    taab_name <- isolate(input$taab_name4)
    ## 1)positive
    cor_positive <- corlist_web %>%
      filter(cor_coef>0) %>%
      filter(target_gene == taab_name) %>%
      dplyr::select(GeneCards,GeneNames,cor_coef) %>%
      dplyr::top_n(input$mumber_top,cor_coef)  %>%
      dplyr::arrange(desc(cor_coef)) 
    
    output$positive_table <- renderDataTable(
      cor_positive, escape = FALSE, 
      selection = "none",
      rownames = FALSE,
      options = list(pageLength = 20, scrollX = TRUE))
    
    ## 2) negtive
    cor_negative <- corlist_web %>%
      filter(cor_coef<0) %>%
      filter(target_gene == taab_name) %>%
      dplyr::select(GeneCards,GeneNames,cor_coef) %>%
      dplyr::arrange(cor_coef) %>%
      dplyr::top_n(input$mumber_top,abs(cor_coef)) 
    
    output$negative_table <- renderDataTable(
      cor_negative, escape = FALSE, 
      selection = "none",
      rownames = FALSE,
      options = list(pageLength = 20, scrollX = TRUE))
  })
  
  ## 3.5 Methylation
  output$methy_table <- DT::renderDataTable(
    lihc_DMP_taa, escape = FALSE, 
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 20, scrollX = TRUE,
                   columnDefs = list(list(className = 'dt-center', targets = '_all'))))
  
  ## 3.6 DNA mutation
  observeEvent(input$action5, {
    ## 1) target
    genename <- isolate(input$taab_name5)
    maf_target <- taab_maf[taab_maf$Hugo_Symbol==genename,]
    table_target <- maf_target[,c(7:15,36:38)]
    
    ## plot
    output$mutation_plot <- renderPlot({
      lollipopPlot(maf=taab_laml, gene=genename, 
                   labPosSize = 3, legendTxtSize = 1.3,
                   titleSize = c(1.8,0.9), domainLabelSize= 1.4,
                   cBioPortal = FALSE, 
                   AACol=NULL, 
                   #AACol="HGVSp_Short",
                   showMutationRate=TRUE)
    })
    
    output$data_mutation <- DT::renderDataTable(
      table_target,escape = FALSE, 
      selection = "none",
      rownames = FALSE,options = list(pageLength = 20, scrollX = TRUE))  
  })
  
  ## 2) mutation summary 
  output$mutation_plot1 <- renderPlot(
    plotmafSummary(maf = taab_laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE,  fs = 1.2, textSize = 1.1))
  output$mutation_plot2 <- renderPlot(
    oncoplot(maf = taab_laml, top = 20, removeNonMutated = TRUE))
  output$mutation_plot3 <- renderPlot(
    plotTiTv(res = titv(maf = taab_laml, plot = FALSE, useSyn = TRUE)))
  
  ## 4.datasets
  output$used_data <- renderDataTable(
    database, escape = FALSE, 
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 20, scrollX = TRUE,lengthChange = FALSE, searching = FALSE)
    )
}  