---
title: "RBIF111 - Homework 9"
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
library(tibble)
library(plotly)
library(scatterplot3d)

```
  
  
  
# Question 1  
## Description of PCA analysis  
  
Principal Component Analysis (PCA) is a powerful statistical technique used to simplify complex, high-dimensional data sets by identifying the most important patterns and reducing data complexity, while preserving as much of the original information as possible.
How it works:  
1. Data Transformation: at its core, PCA is about transforming the original data set into a new coordinate system. Imagine there is a cloud of data points in a multidimensional space. PCA finds the directions (principal components) along which the data varies the most.  
2. Variance Maximization: the first principal component is the direction that explains the largest amount of variance in the data. This is like finding the "best-fit line" that captures the most spread-out variation, but in multiple dimensions. The second principal component is orthogonal (i.e., perpendicular) to the first and captures the next most important pattern of variation.  
3. Dimensionality Reduction: by ranking these principal components by the amount of variance they explain, it can decide how many components to keep. For example, it is possible that just two or three components capture 80-90% of the total variation in the original data set.  

Mathematically, PCA uses *eigenvalue* decomposition of the data's covariance matrix. Therefore, each principal component is an *eigenvector*. The corresponding *eigenvalue* represents how much variance that component explains. Some of the advantages are reduced data complexity, helps to identify hidden patterns, useful for visualization, and it is the pre-processing step for many machine learning algorithms.  

References:
[1] Jolliffe, I. T. (2002). Principal Component Analysis (2nd ed.). Springer.
[2] Lever, J., Krzywinski, M., & Altman, N. (2017). Principal component analysis. Nature Methods, 14(7), 641-642.
  
  
# Question 2  
## PCA analysis  
### Non-normalized data set  
  
The data set obtained from GEO is a study to investigate how lemafulin (as part of a combination therapy strategy) overcomes the sorafenib resistance in hepatocellular carcinoma (HCC). The data set is a RNA-seq on the HepG2 cells treated with DMSO (Control) or lemafulin (Tx).  
  
The (hidden) code below will download the data set as-is from GEO, i.e., not normalized.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
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
RNAseqData <- GEODataDownload(DS = "GSE252988",
                              gpl = "GPL24676",
                              gsm = "000111",
                              PlateAnnotInfo = "Only needed if the data was generated using MicroArrays",
                              GenerateMetaData = "Only needed if the data was generated using MicroArrays",
                              Technology = "RNAseq")

# Save the RNAseq data set to the working directory
saveRDS(RNAseqData, file = "./RNAseq_GSE252988.rds")

# Reload the file into R environment
GSE252988 <- readRDS("./RNAseq_GSE252988.rds")
head(GSE252988)
```
  
  
The code below will convert the data into long format to add the sample annotation of "Ctrl" and "Tx." Then the variance for each gene is calculated and sorted in descending order to obtain the top 100 and top 1000 genes. Finally, the data will be converted back into wide format in preparation for PCA analysis.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Get gene expressions from a control and a treatment sample
GSE252988_2 <- as.data.frame(GSE252988[, c(2, 7:12)])
head(GSE252988_2)

# Convert to long format to add annotation
new_df <- data.frame(
  Sample_ID = c("GSM8011934","GSM8011935", "GSM8011936", "GSM8011937", "GSM8011938", "GSM8011939"), 
  group = c("Ctrl", "Ctrl", "Ctrl", "Tx", "Tx", "Tx")
)

GSE252988_3 <- reshape2::melt(GSE252988_2, id.vars = c("SYMBOL"), variable.name = "Sample_ID", value.name = "Expression")

# merge with sample annotation
GSE252988_4 <- as.data.frame(merge(GSE252988_3, new_df, by = "Sample_ID"))
head(GSE252988_4)

# Calculate variance for each gene
GSE252988_5 <- GSE252988_4 %>%
  group_by(SYMBOL) %>%
  summarize(Variance = var(Expression), .groups = "drop")
head(GSE252988_5)

# Get top 100 & 1000 genes with greatest variance
top100 <- GSE252988_5 %>% arrange(desc(Variance)) %>% head(100)  # get top 100 by variance
top1000 <- GSE252988_5 %>% arrange(desc(Variance)) %>% head(1000)  # get top 1000 by variance

# Get data for top 100 & 1000 genes
top100_data <- GSE252988_4 %>% filter(SYMBOL %in% top100$SYMBOL) # get expression values
head(top100_data)

top1000_data <- GSE252988_4 %>% filter(SYMBOL %in% top1000$SYMBOL) # get expression values
tail(top1000_data)

# Function to convert back to wide format
data_wide <- function(data) {
  data %>%
    select(-group) %>%  # Remove group column before pivoting
    pivot_wider(names_from = SYMBOL, values_from = Expression) %>%
    column_to_rownames(var = "Sample_ID") %>%  # Convert sample IDs to row names
    as.data.frame()  # Output as a data frame
}

# Execute the wide format function
top100_df <- data_wide(top100_data)
top1000_df <- data_wide(top1000_data)

```
  
  
The code below will perform the PCA analysis.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Get the group annotation to add into PCA analysis
group_annotate <- GSE252988_4 %>% distinct(Sample_ID, group) %>% pull(group)

# Function for PCA analysis
my_pca <- function(data, annotation) {
  pca_result <- prcomp(data, scale. = TRUE)

  # explain variance
  var_explained <- summary(pca_result)$importance[2, ]
  cumulative_var <- cumsum(var_explained)

  # get the number of PCs that explain 75% of the variance
  num_pcs <- which(cumulative_var >= 0.75)[1]

  # get scores and annotation
  pca_scores <- as.data.frame(pca_result$x)
  pca_scores$Group <- annotation

  return(list(
    pca_result = pca_result,
    pca_scores = pca_scores,
    num_pcs = num_pcs,
    cumulative_variance = cumulative_var
  ))
}

# Execute PCA analysis function
top100_pca <- my_pca(top100_df, group_annotate)
top1000_pca <- my_pca(top1000_df, group_annotate)

```
  
  
The (hidden) code below will visualize the results of PCA analysis on the top 100 and top 1000 genes using a 3D graph.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
## Visualize PCA results
# Function to create 3D PCA plots
plot_pca_3d <- function(pca_scores, title, colors) {
  # Create the 3D scatter plot
  s3d <- scatterplot3d(
    pca_scores$PC1, 
    pca_scores$PC2, 
    pca_scores$PC3,
    color = colors[as.factor(pca_scores$Group)],
    pch = 16,
    main = title,
    xlab = "PC1",
    ylab = "PC2",
    zlab = "PC3",
    angle = 40,  # Rotation angle
    box = TRUE
  )
  
  # Add legend
  legend("right", 
         legend = levels(as.factor(pca_scores$Group)),
         col = colors,
         pch = 16,
         title = "Groups")
}

# Set up plotting layout
par(mfrow = c(2, 1), mar = c(4, 4, 4, 4))

# Plot top 100 PCA
plot_pca_3d(
  top100_pca$pca_scores,
  "PCA of Top 100 Genes by Variance",
  c("black", "magenta")
)

# Plot top 1000 PCA
plot_pca_3d(
  top1000_pca$pca_scores,
  "PCA of Top 1000 Genes by Variance",
  c("blue", "darkorange")
)

# Reset plotting parameters
par(mfrow = c(1, 1))
```
  
  
**Discussion**  
The PCA analysis on the top 100 genes with the highest variance shows a clear separation of the two groups (Ctrl and Tx) along the primary principal components. The smaller number of genes seems sufficient to capture the major variance in the data, as evidenced by the distinct clustering. Grouping is visible, but the variance may be less distributed.  
Including the top 1000 genes provides a more comprehensive analysis of the data set. The separation between Ctrl and Tx is still apparent, but distributed across a larger number of principal components. The larger gene subset introduces additional dimensions of variance, which might result in more nuanced grouping of samples. This plot suggests that while more variance is captured by increasing the number of genes, much of the separation between the groups can still be observed with a smaller subset.  
The number of principal components explaining most of the variance increases with the inclusion of more genes, from 100 to 1000. The 3D scatter plots for both analyses show clustering patterns aligned with the "Group" annotation, with the treatment (Tx) and control (Ctrl) groups occupying different/separate regions in the 3D space. Using only the top 100 genes seems to be sufficient to capture the major biological variation associated with the grouping.  
  
  
### Normalized data set  
  
The (hidden) code below will download the normalized data set from GEO.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Download normalized data set from GEO
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
RNAseqData2 <- GEODataDownload(DS = "GSE252988",
                              gpl = "GPL24676",
                              gsm = "000111",
                              PlateAnnotInfo = "Only needed if the data was generated using MicroArrays",
                              GenerateMetaData = "Only needed if the data was generated using MicroArrays",
                              Technology = "RNAseq")

# Save the RNAseq data set to the working directory
saveRDS(RNAseqData2, file = "./RNAseq_GSE252988_norm.rds")

# Reload the file into R environment
GSE252988_v2 <- readRDS("./RNAseq_GSE252988_norm.rds")

# Clean up the data set
colnames(GSE252988_v2) <- gsub(".y$", "", colnames(GSE252988_v2))
remove_cols <- grepl(".x$", colnames(GSE252988_v2))
GSE252988_v2.2 <- subset(GSE252988_v2, select = !remove_cols)
head(GSE252988_v2.2)

```
  
The code below will convert the data into long format to add the sample annotation of "Ctrl" and "Tx." Then the variance for each gene is calculated and sorted in descending order to obtain the top 100 and top 1000 genes. Finally, the data will be converted back into wide format in preparation for PCA analysis.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Get gene expressions from a control and a treatment sample
GSE252988_v2.3 <- as.data.frame(GSE252988_v2.2[, c(2, 7:12)])
head(GSE252988_v2.3)

# Convert to long format to add annotation
new_df2 <- data.frame(
  Sample_ID = c("GSM8011934","GSM8011935", "GSM8011936", "GSM8011937", "GSM8011938", "GSM8011939"), 
  group = c("Ctrl", "Ctrl", "Ctrl", "Tx", "Tx", "Tx")
)

GSE252988_v2.4 <- reshape2::melt(GSE252988_v2.3, id.vars = c("SYMBOL"), variable.name = "Sample_ID", value.name = "Expression")

# merge with sample annotation
GSE252988_v2.5 <- as.data.frame(merge(GSE252988_v2.4, new_df2, by = "Sample_ID"))
head(GSE252988_v2.5)

# Calculate variance for each gene
GSE252988_v2.6 <- GSE252988_v2.5 %>%
  group_by(SYMBOL) %>%
  summarize(Variance = var(Expression), .groups = "drop")
head(GSE252988_v2.6)

# Get top 100 & 1000 genes with greatest variance
top100_v2 <- GSE252988_v2.6 %>% arrange(desc(Variance)) %>% head(100)  # get top 100 by variance
top1000_v2 <- GSE252988_v2.6 %>% arrange(desc(Variance)) %>% head(1000)  # get top 1000 by variance

# Get data for top 100 & 1000 genes
top100v2_data <- GSE252988_v2.5 %>% filter(SYMBOL %in% top100_v2$SYMBOL) # get expression values
head(top100v2_data)

top1000v2_data <- GSE252988_v2.5 %>% filter(SYMBOL %in% top1000_v2$SYMBOL) # get expression values
tail(top1000v2_data)

# Function to convert back to wide format
data_wide <- function(data) {
  data %>%
    select(-group) %>%  # Remove group column before pivoting
    pivot_wider(names_from = SYMBOL, values_from = Expression) %>%
    column_to_rownames(var = "Sample_ID") %>%  # Convert sample IDs to row names
    as.data.frame()  # Output as a data frame
}

# Execute the wide format function
top100v2_df <- data_wide(top100v2_data)
top1000v2_df <- data_wide(top1000v2_data)

```
  
The (hidden) code below will perform the PCA analysis.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Get the group annotation to add into PCA analysis
group_annotate <- GSE252988_4 %>% distinct(Sample_ID, group) %>% pull(group)

# Function for PCA analysis
my_pca <- function(data, annotation) {
  pca_result <- prcomp(data, scale. = TRUE)

  # explain variance
  var_explained <- summary(pca_result)$importance[2, ]
  cumulative_var <- cumsum(var_explained)

  # get the number of PCs that explain 75% of the variance
  num_pcs <- which(cumulative_var >= 0.75)[1]

  # get scores and annotation
  pca_scores <- as.data.frame(pca_result$x)
  pca_scores$Group <- annotation

  return(list(
    pca_result = pca_result,
    pca_scores = pca_scores,
    num_pcs = num_pcs,
    cumulative_variance = cumulative_var
  ))
}

# Execute PCA analysis function
top100v2_pca <- my_pca(top100v2_df, group_annotate)
top1000v2_pca <- my_pca(top1000v2_df, group_annotate)

```
  
The (hidden) code below will visualize the results of the PCA analysis.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
## Visualize PCA results
# Function to create 3D PCA plots
plot_pca_3d <- function(pca_scores, title, colors) {
  # Create the 3D scatter plot
  s3d <- scatterplot3d(
    pca_scores$PC1, 
    pca_scores$PC2, 
    pca_scores$PC3,
    color = colors[as.factor(pca_scores$Group)],
    pch = 16,
    main = title,
    xlab = "PC1",
    ylab = "PC2",
    zlab = "PC3",
    angle = 40,  # Rotation angle
    box = TRUE
  )
  
  # Add legend
  legend("right", 
         legend = levels(as.factor(pca_scores$Group)),
         col = colors,
         pch = 16,
         title = "Groups")
}

# Set up plotting layout
par(mfrow = c(2, 1), mar = c(4, 4, 4, 4))

# Plot top 100 PCA
plot_pca_3d(
  top100v2_pca$pca_scores,
  "PCA of Top 100 Genes by Variance (normalized data)",
  c("black", "magenta")
)

# Plot top 1000 PCA
plot_pca_3d(
  top1000v2_pca$pca_scores,
  "PCA of Top 1000 Genes by Variance (normalized data)",
  c("blue", "darkorange")
)

# Reset plotting parameters
par(mfrow = c(1, 1))
```
  
  
**Discussion**  
The top 100 genes plot show that the groups tend to cluster separately, indicating that the top 100 variance can effectively differentiate the two groups, Ctrl and Tx. The separation is most obvious along the PC1 axis, indicating that the PC1 captures the primary source of variation that separates the control and treatment samples.  
The 1000 genes plot show a similar pattern to the top 100 gene plot, but it does have a larger separation distances between the Ctrl and Tx samples. This suggests that expanding the gene set from 100 to 1000 has increased the ability to discriminate between the two groups based on the gene expression. However, due to similar overall pattern between the two plots, either would be effective in separating the control and treatment samples.  
  
  
  
# Question 3  
## Hierchical clustering (default) analysis and repeat with other clustering metrics  
  
  
The code below are the functions for perform the default hierchical clustering analysis and additional analysis using several clustering metrics. The tested selected clustering metrics are: ward, average, and centroid.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Clustering methods
cluster_methods <- c("ward.D", "average", "centroid")

# Function to perform clustering
my_cluster <- function(data, method, metric) {
  # Drop group column
  if ("group" %in% colnames(data)) {
    data <- data[, -which(colnames(data) == "group")]
  }
    # Convert data to numeric matrix
  data_matrix <- as.matrix(data)
  
  if (metric == "spearman") {
    # Correlate genes (columns) for Spearman
    dist_matrix <- as.dist(1 - cor(data_matrix, method = "spearman"))
  } else {
    # Cluster genes (columns) for Euclidean
    dist_matrix <- dist(t(data_matrix), method = metric)
  }  
  return(hclust(dist_matrix, method = method))
}

# Function to perform clustering and create plots for each method
perform_clustering <- function(data, gene_set_name) {
  # Create a list to store all plots
  plot_list <- vector("list", length(cluster_methods) + 1)
  names(plot_list) <- c("spearman", cluster_methods)

  # Spearman correlation clustering
  spearman_cluster <- my_cluster(data, "ward.D2", "spearman")
  plot_list[["spearman"]] <- plot(spearman_cluster,
                                 main = paste("Spearman Correlation -", gene_set_name),
                                 xlab = "", ylab = "Height")
  
  # Euclidean distance-based clustering methods
  for (i in cluster_methods) {
    cluster_result <- my_cluster(data, i, "euclidean")
    plot_list[[i]] <- plot(cluster_result,
                          main = paste(i, "-", gene_set_name),
                          xlab = "", ylab = "Height")
  }
  return(plot_list)
}

# Execute function for top 100 
cluster_100 <- perform_clustering(top100_df, "Top 100 Genes")

# Set up plotting layout for multiple plots
# Adjust numbers based on how many plots you want to show at once
par(mfrow = c(5, 1))

# Display plots for top 100 genes
for (i in seq_along(cluster_100)) {
  print(cluster_100[[i]])
}

```
  
The code below will continue the clustering analysis for the top 1000 genes.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Execute function for top 1000 
cluster_1000 <- perform_clustering(top1000_df, "Top 1000 Genes")

# Set up plotting layout for multiple plots
# Adjust numbers based on how many plots you want to show at once
par(mfrow = c(5, 1), mar = c(4,4,3,1))

# Display plots for top 1000 genes
for (i in seq_along(cluster_1000)) {
  print(cluster_1000[[i]])
}
```
  
The code below will continue the clustering analysis for all of the genes in the data set.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Get all genes
all_gene <- GSE252988_2 %>%
  pull(SYMBOL)

# Function to filter & convert data to wide format
filter_convert <- function(selected_genes) {
  GSE252988_4 %>%
    filter(SYMBOL %in% selected_genes) %>%
    pivot_wider(names_from = SYMBOL, values_from = Expression) %>%
    column_to_rownames("Sample_ID")
}

# Execute the function
data_all <- filter_convert(all_gene)

# Execute clustering function for all genes 
cluster_all <- perform_clustering(data_all, "All Genes")

# Set up plotting layout for multiple plots
# Adjust numbers based on how many plots you want to show at once
par(mfrow = c(5, 1), mar = c(4,4,3,1))

# Display plots for top 1000 genes
for (i in seq_along(cluster_all)) {
  print(cluster_all[[i]])
}

```
  
  
**Discussion**
Since there are several clustering results, the discussion will focus on the comparison of *ward.D* clustering for the top 100 and top 1000 genes. The top 100 genes show distinct high-level clusters, with a large cluster on the left and smaller sub-clusters branching off. Then there are more isolated clusters on the right side. Compared to the top 100 genes, the top 1000 genes plot show a more complex dendrogram with a larger number of high-level clusters. The main difference observed here is that the top 1000 genes plot has a more diverse set of clusters, with a greater number of separate groups branching off at higher levels. This can suggest that increasing the gene set from 100 to 1000 presented more sources of variation. In the top 100 genes the main clusters looks more consolidated, while the top 1000 genes show more distinct separation, which help to show the increased number of genes helped to further differentiate the samples.  
  
Brief description of each clustering method used:  
1. Ward's method: a hierarchical clustering algorithm that aims to minimize the variance within each cluster. It starts with each point as a separate cluster and then iteratively merges the two clusters that result in the smallest increase in total within-cluster variance. This continues until all points are in the same cluster. Ward's method is sensitive to the scaling and metric used for measuring the distances between the points. **Ref**: Großwendt, A., Heiko Röglin and Schmidt, M. (2019). Analysis of Ward’s Method. *Society for Industrial and Applied Mathematics eBooks*, pp.2939–2957. doi:https://doi.org/10.1137/1.9781611975482.182.  
  
2. Average-linkage: a hierarchical clustering method that calculates the distance between two clusters as the average distance between all pairs of observations in the clusters. The algorithm then merges the two clusters with the smallest average distance. Average linkage method is one of the several methods of hierarchical clustering that researchers use to detect the patterns in the data. This method tends to produce compact, well-separated clusters and is less sensitive to outliers compared to single linkage. **Ref**: Saraçli, S., Doğan, N. and Doğan, İ., 2013. Comparison of hierarchical cluster analysis methods by cophenetic correlation. *Journal of inequalities and Applications*, 2013, pp.1-8.  
  
3. Centroid-linkage: a hierarchical clustering method that calculates the distance between two clusters by computing the difference between their centroids (the average value of each feature in the cluster). The algorithm then merges the two clusters with the smallest centroid distance. This method is also known as the "UPGMA" (unweighted pair-group method with arithmetic mean) clustering method. Like other hierarchical clustering methods, centroid linkage identifies similarities between observations and groups them into hierarchical structures or clusters.**Ref**: Saraçli, S., Doğan, N. and Doğan, İ., 2013. Comparison of hierarchical cluster analysis methods by cophenetic correlation. *Journal of inequalities and Applications*, 2013, pp.1-8.  
  
  
  
  
# Question 4  
### Identification of most stable dendogram across exponentially increasing slices of the data set sorted by row variance.  
  
The data set used for the code below is the data frame of the top 100 genes (top100_df), in wide format, generated in question #2 (line 248). In the first lines of the code below, the corresponding sample IDs will be added to the data frame since it only contains gene names and gene expression values. Then the variance will be calculated and the genes will be sorted in descending variance order.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Add sample IDs to the selected data
top100_df2 <- top100_df
top100_df2$Sample_ID <- rownames(top100_df)
top100_df2 <- top100_df2 %>% relocate(Sample_ID, .before = "A2M")

set.seed(42)
# Step 1: Calculate variance for each gene & order by descending variance
gene_var <- apply(top100_df2[, -1], 2, var)
sorted_genes <- top100_df2[, c(1, order(gene_var, decreasing = TRUE) + 1)]
head(gene_var)
```
  
The code below will generate exponentially increasing slices of the data ordered by row variance.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Step 2: Make exponentially increasing slices
# define slice sizes
total_genes <- ncol(sorted_genes)
base_numbers <- c(2, 4, 16, 256)
slice_sizes <- base_numbers^2
slice_sizes <- slice_sizes[slice_sizes <= total_genes]  # don't let sizes > total genes
slice_sizes <- c(slice_sizes, total_genes - 1) # reduce by 1 for the sample ID column
print(slice_sizes)

```
  
The code below will perform the clustering analysis by *ward.D* for each slice generated above.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Step 3: Perform clustering for each slice
# empty list to store clustering results
dendrograms <- list()
# for loop to perform clustering for each slice size
for (size in slice_sizes) {
  # subset data to include only gene values
  subset_data <- t(sorted_genes[, 2:(size + 1)])
  # calculate distance matrix
  dist_matrix <- dist(subset_data)
  # perform hierarchical clustering using "ward.D"
  dendrograms[[as.character(size)]] <- hclust(dist_matrix, method = "ward.D")
}

head(dendrograms)
```
  
Next, for each slice the tree is cut at k=3, and the number of times individual tree leaves (samples) are identified in the same branch across all slices are counted.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Step 4: Cut trees at specific height
# empty list to store clustering memberships
membership_counts <- list()
# for loop to cut trees at specific height
for (size in names(dendrograms)) {
  # cut the tree at a specific height
  memberships <- cutree(dendrograms[[size]], k = 3)  # change k to adjust clustering
  # store memberships
  membership_counts[[size]] <- memberships
}
print(membership_counts)

# Compare leaves (samples) and their groupings across slices
# create matrix to store common leaf counts
common_leaf_counts <- matrix(0, nrow = nrow(sorted_genes), ncol = nrow(sorted_genes))
# set the row and column names
rownames(common_leaf_counts) <- colnames(common_leaf_counts) <- sorted_genes$Sample_ID
# for loop to compare leaves (samples)
for (size in names(membership_counts)) {
  current_membership <- membership_counts[[size]]
  # for each pair of samples, check if they're in the same cluster
  for (i in 1:(nrow(sorted_genes) - 1)) {
    for (j in (i + 1):nrow(sorted_genes)) {
      # check if both values exist and are equal
      if (!is.na(current_membership[i]) &&
          !is.na(current_membership[j]) &&
          current_membership[i] == current_membership[j]) {
        common_leaf_counts[i, j] <- common_leaf_counts[i, j] + 1
        common_leaf_counts[j, i] <- common_leaf_counts[j, i] + 1
      }
    }
  }
}

print(common_leaf_counts)
```
  
Finally, the code below will find the dendrogram with the most stable leaves.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Find the dendrogram with the most stable leaves
# get the scores
stability_scores <- rowSums(common_leaf_counts)
print(stability_scores)

# get the most stable dendrogram
most_stable <- which.max(stability_scores)
print(most_stable)
```
  
The dendrogram and the histograms for each slice are shown below.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
par(mfrow = c(2, 2))
# Plot the most stable dendrogram
most_stable_dendrogram <- dendrograms[[names(dendrograms)[most_stable]]]
plot(most_stable_dendrogram, main = "Most Stable Dendrogram")

# Make histograms for each slice
for (size in names(membership_counts)) {
  hist(membership_counts[[size]], main = paste("Cluster Distribution -", size, "Genes"),
       xlab = "Cluster", col = "skyblue", breaks = 10)
}
```
