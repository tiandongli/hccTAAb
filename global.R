
## step 0：install packages
#cran_packages <- c("shiny", "shinydashboard", "shinydashboardPlus", "plotly", "survival", "survminer", "ggthemes", "tidyverse", "ggplot2","ggpubr", "ggsignif", "ggsci", "DT", "data.table") 
#install.packages(cran_packages)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("maftools")

## step I: libarary packages
if(T){
  library(shiny)
  library(shinydashboard)
  library(shinydashboardPlus)
  library(dashboardthemes)
  library(plotly)
  library(survival)
  library(survminer)
  library(ggthemes)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(ggsignif)
  library(ggsci)
  library(DT)
  #library(data.table)
  library(maftools)
  #library(future)
}

#plan(sequential, split = TRUE)
#no_cores <- availableCores() - 1
#plan(multicore, workers = no_cores)

## step II: load data
#data_files <- list.files("data/", pattern = "\\.Rdata$", full.names = TRUE)

#for (file in data_files) {
  #load(file)
#}

load("data/taa_ID.Rdata")




