---
title: "RBIF111 - Homework 8"
author: "Neshita Kipuw"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 4
    fig_caption: yes
    code_folding: hide
    number_sections: true
  pdf_document:
    toc: true
vignette: >
    %\VignetteIndexEntry{Text}
    %\usepackage[utf8]{inputenc}
    %\VignetteEngine{knitr::rmarkdown}
fontsize: 15pt
editor_options: 
  chunk_output_type: console
---

<style>
pre code, pre, code {
  white-space: pre !important;
  overflow-x: auto !important;
  word-break: keep-all !important;
  word-wrap: initial !important;
}
body {
text-align: justify}
</style>
---
  

```{r library, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Load necessary packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(stats)
library(DESeq2)
library(reshape2)
library(data.table)
library(survival)
```
  
  
  
# Question 1  
## Logistic regression analysis on gene expression levels and the categorical variable.  
#### From normalized data set  
  
The data set obtained from GEO is a study on gene expression profile characterization and gene-of-interest identification in normal, stenotic (AS), and regurgitant (AI) human aortic valves using RNA sequencing. The metadata contain several categorical variables such as gender, age, BMI, status of disease, and left-ventricle ejection fraction (to classify heart failure). The logistic regression perform analysis between gene expression and gender to identify which genes have the most significant effect on patient gender.  
  
The code and results (hidden) below will download the normalized data set from GEO database. *May collapse the code and results*  
  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='hide', cache=FALSE}
#GEO download
GEODataDownload <- function(DS, gpl, gsm, PlateAnnotInfo, GenerateMetaData, Technology){
  library(DESeq2); library(data.table)
  if(Technology == "Array"){
    gset <- getGEO(DS)
    if(length(gset) > 1) idx <- grep(gpl, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    fvarLabels(gset) <- make.names(fvarLabels(gset))
    comp <- gsub(" ", "", gsm)
    comp <- gsub(",", "", comp)
    gsms <- paste0(comp)
    #### Set up raw names ####
    sml <- c()
    for(i in 1:nchar(gsms)){ sml[i] <- substr(gsms,i,i)}
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] < 2)
    if(LogC){ ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) }
    sml <- paste("G", sml, sep="")
    f1 <- as.factor(sml)
    gset$description2 <- f1
    design <- model.matrix(~description2 + 0, gset)
    colnames(design) <- levels(f1)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(G1-G0, levels = design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust="fdr", sort.by = "B", number = 25000000000)
    #### subset ####
    ex2 <- data.table(subset(tT, select=c("ID", "logFC", "P.Value", "adj.P.Val")))
    ex2$ID <- as.character(ex2$ID)
    #### annotate with gene names ####
    plat <- PlateAnnotInfo[GPLID == gpl,][,!"GPLID", with = FALSE]
    if(nrow(plat) == 0){ print(paste("There is no annotation information available for", gpl)) }
    plat$ID <- as.character(plat$ID)
    plat <- plat[!duplicated(plat$ID),]
    ex2 <- merge(plat, ex2, by = "ID")
    ex2$ID <- as.character(ex2$ID)
    exraw <- data.table(ex)
    exraw$ID <- as.character(rownames(ex))
    #### annotate raw data with gene names ####
    ex2 <- merge(ex2, exraw, by = "ID")
    #### Generate Meta Data ####
    Pdat <- pData(gset)
    #### Add Meta data ####
    Pdat <- as.data.table(Pdat)
    return(list(Data = ex2, MetaData = Pdat))
  }
  
  if(Technology == "RNAseq"){
    ACC <- paste("acc=", DS, sep = "")
    file <- paste("file=", DS, "_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep = "")
    comp <- gsub(" ", "", gsm)
    comp <- gsub(",", "", comp)
    gsms <- paste0(comp)
    #### Set up DEG names ####
    urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
    path <- paste(urld, ACC, file, sep="&");
    tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
    exraw <- tbl 
    apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
    annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
    rownames(annot) <- annot$GeneID
    sml <- strsplit(gsms, split="")[[1]]
    sel <- which(sml != "X")
    sml <- sml[sel]
    tbl <- tbl[ ,sel]
    gs <- factor(sml)
    groups <- make.names(c("Ctrl", "Tx"))
    levels(gs) <- groups
    sample_info <- data.frame(Group = gs, row.names = colnames(tbl))
    keep <- rowSums( tbl >= 10 ) >= min(table(gs))
    tbl <- tbl[keep, ]
    ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)
    
    # Estimate size factors and get normalized counts
    ds <- estimateSizeFactors(ds)
    # Continue with DESeq analysis
    ds <- DESeq(ds, test="Wald", sfType="poscount")
    normalized_counts <- counts(ds, normalized=TRUE)

    r <- results(ds, contrast=c("Group", groups[2], groups[1]), alpha=0.05, pAdjustMethod ="fdr")
    tT <- r[order(r$padj)[1:length(r$padj)],]
    tT <- merge(as.data.frame(tT), annot, by.x="row.names", by.y="GeneID", sort=F)
    tT <- subset(tT, select=c("Row.names","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean","Symbol","Description"))
    
    #### subset ####
    ex2 <- data.table(subset(tT, select=c("Row.names", "Symbol", "Description", "log2FoldChange", "pvalue", "padj")))
    setnames(ex2, c("Row.names", "Symbol", "Description"), c("ENTREZID", "SYMBOL", "GENENAME"))
    
    #### Get Raw data ####
    exraw <- as.data.table(exraw)
    exraw$ENTREZID <- rownames(exraw)
    
    #### merge FC and raw data together ####
    mer <- merge(ex2, exraw, by = "ENTREZID")
    
    # Add normalized counts to the merged data
    normalized_counts_dt <- as.data.table(normalized_counts, keep.rownames = "ENTREZID")
    final_table <- merge(mer, normalized_counts_dt, by = "ENTREZID")
    
    return(final_table)
  }
}

# Execute the GEODataDownload function and obtain the RNAseq data set
RNAseqData <- GEODataDownload(DS = "GSE153555",
                              gpl = "GPL16791",
                              gsm = "000000000111111111111111111",
                              PlateAnnotInfo = "Only needed if the data was generated using MicroArrays",
                              GenerateMetaData = "Only needed if the data was generated using MicroArrays",
                              Technology = "RNAseq")

# Save the RNAseq data set to the working directory
saveRDS(RNAseqData, file = "~/RNAseq_GSE153555_norm.rds")

# Reload the file into R environment
GSE153555 <- readRDS("~/RNAseq_GSE153555_norm.rds")

# Edit column names that end with ".y"
colnames(GSE153555) <- gsub(".y$", "", colnames(GSE153555))

# Clean up column names
remove_cols <- grepl(".x$", colnames(GSE153555))
GSE153555_2 <- subset(GSE153555, select = !remove_cols)
head(GSE153555_2)
```
  
  
The column below will add columns for the continuous variables associated with data set, obtained from the metadata.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Add in metadata from folder
metadata <- read.csv("./SRA_GSE153555.csv")
metadata <- metadata %>% mutate(degree_of_disease = as.factor(degree_of_disease),
                              Gender = as.factor(Gender))
head(metadata)

# Convert data set to long format to merge with sample annotation
GSE153555_3 <- reshape2::melt(GSE153555_2, id.vars = c("ENTREZID", "SYMBOL", "GENENAME", "log2FoldChange", "pvalue", "padj"), variable.name = "Sample_ID", value.name = "Expression")

# Merge the data frames by Sample_ID
GSE153555_4 <- merge(GSE153555_3, metadata, by = "Sample_ID", all.x = TRUE)

# Check the resulting converted data set
head(GSE153555_4)

```
  
  
The code below will perform the logistic regression analysis on gene expression and the effect on patient gender, and the results will be visualized.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
# Logistic Regression on gene expression and gender
log_reg <- GSE153555_4 %>%
  group_by(SYMBOL) %>%
  summarise(
    pvalue = {
        model <- glm(Gender ~ Expression,
                     family = binomial(link = "logit"),
                     data = pick(Expression, Gender))
        summary(model)$coefficients[2,4]
    }
  )

# Sort by p-value
top_genes <- log_reg %>% arrange(pvalue) %>% slice_head(n = 5)
print("Top 5 genes with the lowest p-values:")
print(top_genes)

# Plot p-values
ggplot(log_reg, aes(x = reorder(SYMBOL, pvalue), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  coord_flip() +
  labs(title = "Gene significance by gender", hjust = 0.5, face = "bold", size = 15,
       x = "Genes", size = 12, face = "bold",
       y = "-log10(p-value)", size = 12, face = "bold") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

```
  
  
### Multiple hypothesis testing and leave-n-out resampling from normalized data set  
  
  
The code below will perform multiple hypothesis testing by n = 100 resampling and leaving out n random samples each time. The top 5 most significant genes by p-value will be identified. The average p-value and variance from resampling will be calculated and reported. The gene with low average p-value and low variance across resampling is considered to be robustly associated with gender. Although not found in this top 5 list, genes with high variance may indicate unstable or sample-specific effects.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
n_resampling <- 100
n_leaveout <- 1

# Split data by gene
GSE153555_5 <- split(GSE153555_4, GSE153555_4$SYMBOL)

## Function for multiple hypothesis test & leave-n-out resampling
resample_pvalue <- function(GSE153555_4, n_resampling, n_leaveout) {
    replicate(n_resampling, {
        # leave out 5 random samples
        validation_sample <- sample(1:nrow(GSE153555_4), size = n_leaveout)
        training_data <- GSE153555_4[-validation_sample, ] # training set
        validation_data <- GSE153555_4[validation_sample, ] # validation set
        
        # fit logistic regression model
        model <- glm(Gender ~ Expression, family = binomial(link = "logit"), data = training_data)
        # get the p-value of the Expression coefficient
        summary(model)$coefficients[2,4]
    })
}

# Execute the function
my_resampling <- sapply(GSE153555_5, function(GSE153555_4) {
  p_values <- resample_pvalue(GSE153555_4, n_resampling, n_leaveout)
  c(avg_pvalue = mean(p_values, na.rm = TRUE), var_pvalue = var(p_values, na.rm = TRUE))
})

# Convert to data frame
my_summary <- as.data.frame(t(my_resampling))
colnames(my_summary) <- c("avg_pvalue", "var_pvalue")
my_summary$SYMBOL <- rownames(my_summary)
head(my_summary)

# Get top 5 results based on average p-values
top_resampling <- my_summary %>%
    arrange(avg_pvalue) %>%
    head(n = 5)

print("Top 5 genes from resampling with the lowest average p-values:")
print(top_resampling)
```
  
  
**Discussion**  
The logistic regression show that the genes with the smallest p-values (< 0.05) indicate a significant effect on gender. The result of resampling help to validate these results by assessing consistency across multiple samples, which showed the same top 5 genes as the original logistic regression. 
  
  
## Repeat logistic regression analysis and multiple hypothesis testing on non-normalized data set  
  
The (hidden) code below will download the data set as-is from GEO, i.e., not normalized.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='hide', cache=FALSE}
#GEO download
GEODataDownload <- function(DS, gpl, gsm, PlateAnnotInfo, GenerateMetaData, Technology){
  library(DESeq2); library(data.table)
  if(Technology == "Array"){
    gset <- getGEO(DS)
    if(length(gset) > 1) idx <- grep(gpl, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    fvarLabels(gset) <- make.names(fvarLabels(gset))
    comp <- gsub(" ", "", gsm)
    comp <- gsub(",", "", comp)
    gsms <- paste0(comp)
    #### Set up raw names ####
    sml <- c()
    for(i in 1:nchar(gsms)){ sml[i] <- substr(gsms,i,i)}
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] < 2)
    if(LogC){ ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) }
    sml <- paste("G", sml, sep="")
    f1 <- as.factor(sml)
    gset$description2 <- f1
    design <- model.matrix(~description2 + 0, gset)
    colnames(design) <- levels(f1)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(G1-G0, levels = design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust="fdr", sort.by = "B", number = 25000000000)
    #### subset ####
    ex2 <- data.table(subset(tT, select=c("ID", "logFC", "P.Value", "adj.P.Val")))
    ex2$ID <- as.character(ex2$ID)
    #### annotate with gene names ####
    plat <- PlateAnnotInfo[GPLID == gpl,][,!"GPLID", with = FALSE]
    if(nrow(plat) == 0){ print(paste("There is no annotation information available for", gpl)) }
    plat$ID <- as.character(plat$ID)
    plat <- plat[!duplicated(plat$ID),]
    ex2 <- merge(plat, ex2, by = "ID")
    ex2$ID <- as.character(ex2$ID)
    exraw <- data.table(ex)
    exraw$ID <- as.character(rownames(ex))
    #### annotate raw data with gene names ####
    ex2 <- merge(ex2, exraw, by = "ID")
    #### Generate Meta Data ####
    Pdat <- pData(gset)
    #### Add Meta data ####
    Pdat <- as.data.table(Pdat)
    return(list(Data = ex2, MetaData = Pdat))
  }
  
  if(Technology == "RNAseq"){
    ACC <- paste("acc=", DS, sep = "")
    file <- paste("file=", DS, "_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep = "")
    comp <- gsub(" ", "", gsm)
    comp <- gsub(",", "", comp)
    gsms <- paste0(comp)
    #### Set up DEG names ####
    urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
    path <- paste(urld, ACC, file, sep="&");
    tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
    exraw <- tbl 
    apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
    annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
    rownames(annot) <- annot$GeneID
    sml <- strsplit(gsms, split="")[[1]]
    sel <- which(sml != "X")
    sml <- sml[sel]
    tbl <- tbl[ ,sel]
    gs <- factor(sml)
    groups <- make.names(c("Ctrl", "Tx"))
    levels(gs) <- groups
    sample_info <- data.frame(Group = gs, row.names = colnames(tbl))
    keep <- rowSums( tbl >= 10 ) >= min(table(gs))
    tbl <- tbl[keep, ]
    ds <- DESeqDataSetFromMatrix(countData=tbl, colData=sample_info, design= ~Group)
    ds <- DESeq(ds, test="Wald", sfType="poscount")
    r <- results(ds, contrast=c("Group", groups[2], groups[1]), alpha=0.05, pAdjustMethod ="fdr")
    tT <- r[order(r$padj)[1:length(r$padj)],]
    tT <- merge(as.data.frame(tT), annot, by=0, sort=F)
    tT <- subset(tT, select=c("GeneID","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean","Symbol","Description"))
    #### subset ####
    ex2 <- data.table(subset(tT, select=c("GeneID", "Symbol", "Description", "log2FoldChange", "pvalue", "padj")))
    #### Adjust column names ####
    setnames(ex2, c("GeneID", "Symbol", "Description"), c("ENTREZID", "SYMBOL", "GENENAME"))
    #### Get Raw data ####
    GeneID <- as.integer(rownames(exraw))
    exraw <- as.data.table(exraw)
    #### Update column names ####
    exraw$ENTREZID <- GeneID
    #### merge FC and raw data together ####
    mer <- merge(ex2, exraw, by = "ENTREZID")
    return(mer)  }
}

# Execute the GEODataDownload function and obtain the RNAseq data set
RNAseqData2 <- GEODataDownload(DS = "GSE153555",
                              gpl = "GPL16791",
                              gsm = "000000000111111111111111111",
                              PlateAnnotInfo = "Only needed if the data was generated using MicroArrays",
                              GenerateMetaData = "Only needed if the data was generated using MicroArrays",
                              Technology = "RNAseq")

# Save the RNAseq data set to the working directory
saveRDS(RNAseqData2, file = "~/RNAseq_GSE153555.rds")

# Reload the file into R environment
GSE153555v2 <- readRDS("~/RNAseq_GSE153555.rds")

# Edit column names that end with ".y"
#colnames(GSE153555v2) <- gsub(".y$", "", colnames(GSE153555v2))

# Clean up column names
#remove_cols <- grepl(".x$", colnames(GSE153555v2))
#GSE153555v2_2 <- subset(GSE153555v2, select = !remove_cols)
head(GSE153555v2)
```

The column below will add columns for the continuous variables associated with data set, obtained from the metadata.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Grab previously downloaded metadata
head(metadata)

# Convert data set to long format to merge with sample annotation
GSE153555v2_2 <- reshape2::melt(GSE153555v2, id.vars = c("ENTREZID", "SYMBOL", "GENENAME", "log2FoldChange", "pvalue", "padj"), variable.name = "Sample_ID", value.name = "Expression")

# Merge the data frames by Sample_ID
GSE153555v2_3 <- merge(GSE153555v2_2, metadata, by = "Sample_ID", all.x = TRUE)

# Check the resulting converted data set
head(GSE153555v2_3)

```
  
      
The code below will perform the logistic regression analysis on gene expression and the effect on patient gender, and the results will be visualized.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
# Logistic Regression on gene expression and gender
log_reg2 <- GSE153555v2_3 %>%
  group_by(SYMBOL) %>%
  summarise(
    pvalue = {
        model <- glm(Gender ~ Expression,
                     family = binomial(link = "logit"),
                     data = pick(Expression, Gender))
        summary(model)$coefficients[2,4]
    }
  )

# Sort by p-value
top_genes2 <- log_reg2 %>% arrange(pvalue) %>% slice_head(n = 5)
print("Top 5 genes with the lowest p-values:")
print(top_genes2)

# Plot p-values
ggplot(log_reg2, aes(x = reorder(SYMBOL, pvalue), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  coord_flip() +
  labs(title = "Gene significance by gender", hjust = 0.5, face = "bold", size = 15,
       x = "Genes", size = 12, face = "bold",
       y = "-log10(p-value)", size = 12, face = "bold") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

```
  
  
### Multiple hypothesis testing and leave-one-out resampling from non-normalized data set  
  
The code below will perform multiple hypothesis testing by n = 100 resampling and leaving out n random sample each time. The top 5 most significant genes by p-value will be identified. The average p-value and variance from resampling will be calculated and reported.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
## Recall the previous variables
# n_resampling <- 100
# n_leaveout <- 1

# Split data by gene
GSE153555v2_4 <- split(GSE153555v2_3, GSE153555v2_3$SYMBOL)

## Function for multiple hypothesis test & leave-n-out resampling
resample_pvalue2 <- function(GSE153555v2_3, n_resampling, n_leaveout) {
    replicate(n_resampling, {
        # leave out 5 random samples
        validation_sample2 <- sample(1:nrow(GSE153555v2_3), size = n_leaveout)
        training_data2 <- GSE153555v2_3[-validation_sample2, ] # training set
        validation_data2 <- GSE153555v2_3[validation_sample2, ] # validation set
        
        # fit logistic regression model
        model2 <- glm(Gender ~ Expression, family = binomial(link = "logit"), data = training_data2)
        # get the p-value of the Expression coefficient
        summary(model2)$coefficients[2,4]
    })
}

# Execute the function
my_resampling2 <- sapply(GSE153555v2_4, function(GSE153555v2_3) {
  p_values2 <- resample_pvalue2(GSE153555v2_3, n_resampling, n_leaveout)
  c(avg_pvalue = mean(p_values2, na.rm = TRUE), var_pvalue = var(p_values2, na.rm = TRUE))
})

# Convert to data frame
my_summary2 <- as.data.frame(t(my_resampling2))
colnames(my_summary2) <- c("avg_pvalue", "var_pvalue")
my_summary2$SYMBOL <- rownames(my_summary2)
head(my_summary2)

# Get top 5 results based on average p-values
top_resampling2 <- my_summary2 %>%
    arrange(avg_pvalue) %>%
    head(n = 5)

print("Top 5 genes from resampling with the lowest average p-values:")
print(top_resampling2)

```
  
  
**Discussion**  
The logistic regression for non-normalized data set also show small p-values < 0.05, indicating significant effect on gender. The resampling helps to validate the top 5 genes results from logistic regression, showing a similar list. Despite the low average p-values and small variance of p-values from resampling, only 2 genes (EIF1AY and XIST) out of the top 5 most significant genes were identified in both non- and normalized data sets. Therefore, careful examination of the original data should be considered to determine whether normalization would be required, as there is quite a notable difference between non- and normalized data, which impact the results greatly. Non-normalized data will have large differences in scale, leading to certain genes dominating the model due to extreme values, even if these differences are potentially not meaningful. Normalization then reduces this effect, allowing for more reliable comparisons and reducing the variability. Additionally, non-normalized data may display skewness or non-uniform distribution across samples, while normalized data is closer to a Gaussian distribution thereby the significant genes are more likely to reflect genuine biological results.  
  
  
  
# Question 2  
## Time-to-event survival analysis  
  
The data was downloaded from the SurvSet Library, and contains information on the survival of subjects who are habitual smokers for a time period measured in years. The data contains information such as age, gender, the number of years the subject has been smoking, etc. The code below will load the data set into the R environment to perform the survival analysis. THen the data will be further analyzed to take into account the gender of the subjects.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Load the data set
mydata <- read.csv("./pharmacoSmoking.csv")
head(mydata)
summary(mydata)

# Analyse the data
survival_data <- Surv(time = mydata$num_yearsSmoking, event = mydata$event)
head(survival_data)

# Fit Kaplan-Meier survival curves
km_fit <- survfit(survival_data ~ mydata$fac_gender)
summary(km_fit)
```
  
  
The (hidden) code below will plot the survival curves for the survival probability vs. years of smoking, classified by gender.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Plot survival curves
plot(km_fit, col = c("blue", "red"), 
     xlab = "Years of Smoking", 
     ylab = "Survival Probability",
     main = "Survival Curves for Smoking Status",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 1.8,
     lwd = 2)
legend("topright", 
       legend = levels(as.factor(mydata$fac_gender)), 
       col = c("blue", "red"),
       cex = 1.5,
       lwd = 3)
```