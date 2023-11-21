####################### hccTAAb Atlas #######################################
## code for data pre-processing --------------------------------------------

## Information Links  ----------------------------------------------------------
## Information Links
## GeneCards/Symbol 
## NCBI/Entrez 
## HPA 
## UniProt
## PubMed
## neXtProt 
# createLink for GeneCards & NCBI & Esemble
createLink <- function(val,name) {
  sprintf('<a href="%s" class="btn btn-link" target="_blank" >%s</a>',val,name) ##target="_blank" 
}

hpaid <- paste0(TAAb_infor$ENSEMBL,"-",TAAb_infor$Symbol)
pubmedid <- paste0(TAAb_infor$Symbol,"+hepatocellular carcinoma")

TAAb_infor_web <- data.table::data.table(
  TAAb = paste0("<a href='https://pubmed.ncbi.nlm.nih.gov/",TAAb_infor$PMID,"'>",TAAb_infor$TAAb),
  Symbol = paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",TAAb_infor$Symbol,"'>",TAAb_infor$Symbol),
  ENSEMBL = paste0("<a href='https://asia.ensembl.org/Homo_sapiens/Gene/Summary?g=",TAAb_infor$ENSEMBL,"'>",TAAb_infor$ENSEMBL,"</a>"),
  UniProt = paste0("<a href='https://www.uniprot.org/uniprotkb/",TAAb_infor$UniProt,"'>",TAAb_infor$UniProt,"</a>"),
  NCBI = paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/",TAAb_infor$NCBI,"'>",TAAb_infor$NCBI,"</a>"),
  HPA = paste0("<a href='https://www.proteinatlas.org/",hpaid,"/pathology","'>","HPA","</a>"),
  neXtProt = paste0("<a href='https://www.nextprot.org/entry/NX_",TAAb_infor$UniProt,"/","'>","neXtProt","</a>"),
  'Related Study' = paste0("<a href='https://pubmed.ncbi.nlm.nih.gov/?term=",pubmedid,"'>","Link","</a>"),
  'Official Full Name' = TAAb_infor$GeneNames,
  'Aliases' = TAAb_infor$GeneAlias
)

save(TAAb_infor_web,file = "TAAb_infor_web.Rdata")


## GEO data download -------------------------------------------------------
## GEO data download: GSE144269, GSE14520, GSE22058, GSE25097, GSE36376, GSE63898
## example: GSE144269

# Step 0: Preparation
setwd("D:/hccTAAb/data")
GEO_ID <- c("GSE144269") #Change the GEO id ! 

if(T){
  rm(list = ls())
  library(GEOquery)
  library(tidyverse)
  if(!dir.exists(paste0("./",GEO_ID))){
    dir.create(paste0("./",GEO_ID))
    setwd(paste0("./",GEO_ID))
  } else {
    setwd(paste0("./",GEO_ID))}
}

# Step I: Download 
if(T){
  if (any(file.exists(list.files(pattern = "\\.gz$")))) {
    ## One: read local data
    print("=> Local data exists, ready to read local data")
    gset <- getGEO(filename=paste0("./",GEO_ID,"_series_matrix.txt.gz"),AnnotGPL=FALSE)
    
  } else {
    print("=> No local data exists, ready to download online data from GEO")
    ## Two: Download online data from GEO
    gset <- getGEO(GEO_ID, 
                   GSEMatrix=TRUE, 
                   AnnotGPL=FALSE,
                   destdir = "./")
  }
  
  #show(gset)
  save(gset,file = paste0("gset_",GEO_ID,".Rdata"))
  
  exprs <- as.data.frame(gset@assayData[["exprs"]])
  pheno <- as.data.frame(gset@phenoData@data)
  
  ## Check if log normalization is required
  if(T){
    qx <- as.numeric(quantile(exprs, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { exprs[which(exprs <= 0)] <- NaN
    exprs <- log2(exprs+1.01)
    print("=> log2 transform finished")}
    else{print("=> log2 transform not needed")}
  }
}

# Step II: Gene Annotation
## Option I: Using bioconductor packages
if(T){
  ## Load the required R packages here
  ## hgu133a.db: GPL96 [HG-U133A] Affymetrix Human Genome U133A Array
  ## hgu219.db: GPL13667 [HG-U219] Affymetrix Human Genome U219 Array
  ## hgu133plus2.db: 
  ## illuminaHumanv4: 
  ## GPL5188:huex10sttranscriptcluster	
  GPL_package_names <- paste0("hgu133a.db") #different GPL has different package.
  
  library(GPL_package_names, character.only = TRUE)
  library(AnnotationDbi)
  library(annotate)
  #class(exprs) 
  #head(exprs)[,1:8]
  print("=> Required R Annotation Packages have been loaded")
  exprs$symbol <- getSYMBOL(rownames(exprs), paste0(GPL_package_names))
  
  #dim(exprs)
  exprs_nona <- exprs[!is.na(exprs$symbol),]  #查看symbol列非NA的行（即完整注释的基因数）
  dim(exprs_nona)
  dim(exprs[is.na(exprs$symbol),]) #查看symbol是NA的行（即没注释上基因的列）
  symbol_dup <- exprs_nona[duplicated(exprs_nona$symbol),] #查看symbol重复的行数
  exprs_nona[exprs_nona$symbol=='',] #检查有没有空字符串
  
  exprs <- exprs_nona %>%
    aggregate(. ~symbol, max) %>%
    column_to_rownames("symbol") 
  
  ## Save data
  save(exprs, pheno, file = paste0("expr_pheno_",GEO_ID,".Rdata"))
}

## gene expression data of 170 taab 
comID <- intersect(idtaa,rownames(exprs))
comID
length(comID)

table(pheno$TYPE)
expr <- exprs[comID,]
expr <- as.data.frame(t(expr))
expr <- expr[rownames(pheno),]
expr$Group=pheno$TYPE
table(expr$Group)
expr$Group <- ifelse(expr$Group=="HCC","T","N")
library(tidyr)
expr_GSE144269 <- pivot_longer(expr,
                               cols = 1:length(comID),
                               names_to = "Gene",
                               values_to = "Expression")
expr$Dataset <- c("GSE144269")


## combined expression data
for (file_name in file_names) {
  load(file_name)
}

merged_data <- rbind(expr_gse144269, expr_gse14520,expr_gse22058,expr_gse25097,
                     expr_gse36376,expr_gse63898,expr_icgc,expr_tcga)

table(merged_data$Dataset)
save(merged_data, file = "merged_data.Rdata")



## Protein -----------------------------------------------------------------
rm(list = ls()) 
library(data.table)
library(tidyverse)
load("taa_ID.Rdata")
## pathology
hpa_pathology_tissue <- fread("pathology.tsv")
table(hpa_pathology_tissue$Cancer)
hpa_pathology_lc <- filter(hpa_pathology_tissue,Cancer %in% c("liver cancer"))

hpa_pathology_lctaa <- filter(hpa_pathology_lc,`Gene name` %in% idtaa)
setdiff(idtaa,hpa_pathology_lc$`Gene name`)
hpa_pathology_lctaa <- hpa_pathology_lctaa[,-c(1,3)]
save(hpa_pathology_lctaa,file = "hpa_pathology_lctaa.Rdata")

## normal
hpa_normal_tissue <- fread("normal_tissue.tsv")
table(hpa_normal_tissue$Tissue)
hpa_normal_lc <- filter(hpa_normal_tissue,Tissue %in% c("liver"))
hpa_normal_lctaa <- filter(hpa_normal_lc,`Gene name` %in% idtaa)
hpa_normal_lctaa <- hpa_normal_lctaa[,-c(1,3)]
save(hpa_normal_lctaa,file = "hpa_normal_lctaa.Rdata")

## Similarity Analysis -----------------------------------------------------
rm(list = ls())
library(tidyverse)
load("./TCGA-LIHC_tpm.Rdata")
comID <- intersect(idtaa,tpm$gene_name)
df <- tpm %>% filter(gene_name != comID)
comID <- intersect(idtaa,df$gene_name)
mRNA <- tpm %>% filter(gene_type=="protein_coding")
comID <- intersect(idtaa,mRNA$gene_name)
length(comID)
diffID <- setdiff(idtaa,comID)
diffID
df <- tpm %>% filter(gene_name %in% diffID)

expr_mRNA <- rbind(mRNA,df)
table(expr_mRNA$gene_type)

dim(expr_mRNA[is.na(expr_mRNA$gene_name),]) 
dim(expr_mRNA[duplicated(expr_mRNA$gene_name),]) 
dim(expr_mRNA[expr_mRNA$gene_name=='',]) 

is.data.frame(expr_mRNA)
class(expr_mRNA$gene_name)
expr_mRNA <- expr_mRNA %>%
  select(-id,-gene_type) %>%
  aggregate(. ~gene_name, max) %>%
  column_to_rownames("gene_name") %>% 
  select(starts_with("TCGA"))

pheno <- data.frame(sample=colnames(expr_mRNA),
                    group=ifelse(str_sub(colnames(expr_mRNA),14,15)=="01","T","N"))

table(pheno$group)
expr_mRNA <- expr_mRNA[,pheno$group=="T"]
expr <- as.data.frame(t(expr_mRNA))

gene_correlation <- function(data, gene) {
  gene_data <- data[,gene]
  other_data <- data[,!colnames(data) == gene, drop=FALSE]
  cor_results <- apply(other_data, 2, function(x) cor.test(gene_data, x, method = "pearson",adjust = "fdr"))
  cor_coef <- sapply(cor_results, function(x) x$estimate)
  p_value <- sapply(cor_results, function(x) x$p.value)
  cor_result_df <- data.frame(
    target_gene = gene,
    gene_all = colnames(other_data),
    Correlation = cor_coef,
    FDR = p_value,
    #adjusted_pvalue = adjusted_pvalue,
    stringsAsFactors = FALSE,
    row.names = NULL
  ) %>%
    na.omit() %>%
    #filter(p_value < 0.05, abs(cor_coef) > 0.5)
    filter(FDR < 0.05) %>%
    return(cor_result_df)
}


load("./idtaa.Rdata")
idtaa
corlist <- lapply(idtaa, gene_correlation, data=expr)
#corlist <- mclapply(idtaa, gene_correlation, data=expr, mc.cores = 4)
corlist_web <- do.call(rbind, corlist)
table(corlist_web$target_gene)
length(table(corlist_web$target_gene))
save(corlist_web, file = "corlist_web.Rdata")

## DNA mutation ------------------------------------------------------------
library(TCGAbiolinks)
library(tidyverse)
library(maftools)

query <- GDCquery(
  project = "TCGA-LIHC", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)
maf_LIHC <- GDCprepare(query)
id <- str_sub(maf_LIHC$Matched_Norm_Sample_Barcode,1,15)
id <- unique(id)
id <- str_sub(id,14,15)
taab_maf <- filter(maf_LIHC,Hugo_Symbol %in% geneid)
taab_maf_laml = read.maf(maf = taab_maf)
save(taab_maf_laml,file = "taab_maf_laml.Rdata")

## DNA Methylation ---------------------------------------------------------
library(TCGAbiolinks)
library(ChAMP)
library(doParallel)
library(foreach)
cl.cores = detectCores()
cl <- makeCluster(cl.cores)
registerDoparallel(cl-2)  
load("taa_ID.Rdata")

query_met.hg38 <- GDCquery(
  project = "TCGA-LIHC", 
  legacy = FALSE,
  data.category = "DNA Methylation", 
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)
GDCdownload(query_met.hg38)
lihc_met_hg38 <- GDCprepare(query_met.hg38)
object.size(lihc_met_hg38)
#save(lihc_met_hg38,file = "lihc_met_hg38.Rdata")


lihc_met_hg38 <- assay(lihc_met_hg38)
#lihc_met_hg38 <- as.data.frame(lihc_met_hg38) 
myLoad <- lihc_met_hg38
colnames(myLoad) <- substr(colnames(myLoad),1,15)
myLoad <- as.data.frame(myLoad)
group <- data.frame(Sample=colnames(myLoad), 
                    Group=ifelse(substr(colnames(myLoad), 14,15)=="01","Tumor","Normal"))
table(group$Group)

myfilter <- champ.filter(as.matrix(myLoad), pd=group) #data filtering 
myfilter$beta <- na.omit(myfilter$beta) 
CpG.GUI(CpG=rownames(myfilter$beta), arraytype="450K") #CpG distribution

lihc_met_Norm <- champ.norm(beta=myfilter$beta, arraytype="450K", cores=8) 
save(lihc_met_Norm, file = "lihc_met_Norm.Rdata") 
lihc_met_Norm[1:5,1:5]
QC.GUI(beta=lihc_met_Norm)

## Finding differentiated CpG 
lihc_DMP <- champ.DMP(beta = lihc_met_Norm, pheno=myfilter$pd$Group) 
lihc_DMP <- lihc_DMP[["Tumor_to_Normal"]] 

lihc_DMP_taa <- lihc_DMP[lihc_DMP$gene %in% idtaa,]
lihc_DMP_taa <- lihc_DMP_taa[,c(14,1:8,15:18)]
lihc_DMP_taa$logFC <- round(lihc_DMP_taa$logFC, 3)
lihc_DMP_taa$Tumor_AVG <- round(lihc_DMP_taa$Tumor_AVG, 3)
lihc_DMP_taa$Normal_AVG <- round(lihc_DMP_taa$Normal_AVG, 3)
lihc_DMP_taa$AveExpr <- round(lihc_DMP_taa$AveExpr, 3)
lihc_DMP_taa$B <- round(lihc_DMP_taa$B, 3)
lihc_DMP_taa$t <- round(lihc_DMP_taa$t, 3)
lihc_DMP_taa$P.Value <- round(lihc_DMP_taa$P.Value, 6)
lihc_DMP_taa$adj.P.Val <- round(lihc_DMP_taa$adj.P.Val, 6)

lihc_DMP_taa$Probs <- rownames(lihc_DMP_taa)
lihc_DMP_taa <- lihc_DMP_taa[,c(1,14,2:13)]
length(unique(lihc_DMP_taa$gene))
lihc_DMP_taa <- lihc_DMP_taa[order(lihc_DMP_taa$gene),]
names(lihc_DMP_taa)[1] <- "Symbol"

save(lihc_DMP_taa, file = "lihc_DMP_taa.Rdata")





