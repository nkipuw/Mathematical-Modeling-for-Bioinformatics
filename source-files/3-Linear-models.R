---
title: "RBIF111 - Homework 3"
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
library(cowplot)
```

The following RNAseq data set is obtained from GEO accession **GSE248274**. It is a study of the mechanisms involved in colonic epithelial differentiation that is key to unraveling the alterations causing inflammatory conditions and colon cancer. The researchers performed global transcriptomic analyses (RNA-seq) in colon normal organoids derived from 6 patients that were incubated in proliferation medium or differentiation medium, containing Bone Morphogenetic Protein 4 and the Notch inhibitor dibenzazepine, for 48 h. Additionally, organoids cultured in differentiation medium were treated with calcitriol or ethanol (control) to investigate the effect of vitamin D on the differentiation of colon normal organoids. *Note*: code is hidden due to extensive lines of code.

```{r GEO, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
GEODataDownload <- function(DS, gpl, gsm, PlateAnnotInfo, GenerateMetaData, Technology){
  library(DESeq2); library(limma); library(data.table)
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
    return(mer)
  }
}
#### Execute the GEODataDownload function and obtain the RNAseq data set ####
RNAseqData <- GEODataDownload(DS = "GSE248274",
                              gpl = "GPL21697",
                              gsm = "01001001010101",
                              PlateAnnotInfo = "Only needed if the data was generated using MicroArrays",
                              GenerateMetaData = "Only needed if the data was generated using MicroArrays",
                              Technology = "RNAseq")
data.table(head(as.data.frame(RNAseqData),10), filter = 'top', options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px", autoWidth = TRUE))

```

File was downloaded and saved into the current working directory under a new name, and reloaded into R environment.

```{r saveData, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
####Define file path to save into working directory####
#getwd()
#file_path <- "~/RBIF111"

####Save the RNAseq data set to the working directory####
saveRDS(RNAseqData, file = "RNAseq_GSE248274.rds")

####Reload the file into R environment####
RNAseq_GSE248274 <- readRDS("RNAseq_GSE248274.rds")
print(head(RNAseq_GSE248274))
```

Added new columns called age, sex, group, and Subset_Comparison for sample annotation

```{r longData, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
####Create a new data frame to group the samples by ctrl vs. Tx ####
new_df <- data.frame(
  Sample_ID = c("GSM7910769", "GSM7910771", "GSM7910772", "GSM7910773", "GSM7910774", "GSM7910775", "GSM7910776", "GSM7910777", "GSM7910778", "GSM7910780", "GSM7910782", "GSM7910783", "GSM7910785", "GSM7910786"),
  group = c("Ctrl", "Tx", "Ctrl", "Ctrl", "Tx", "Ctrl", "Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx"),
  age = c(58,58,55,55,55,64,64,64,71,71,88,88,68,68),
  sex = c("M","M","M","M","M","M","M","M","F","F","F","F","M","M")
)

#### Convert data set to long format to merge the sample group####
RNAseq_GSE248274_2 <- melt(RNAseq_GSE248274, id.vars = c("ENTREZID", "SYMBOL", "GENENAME", "log2FoldChange", "pvalue", "padj"), 
                           variable.name = "Sample_ID", value.name = "Expression")

####Merge the sample group to the long format data set and append the new column #### 
RNAseq_GSE248274_2 <- merge(RNAseq_GSE248274_2, new_df, by = "Sample_ID", all.x = TRUE)

####Append the new column and use a placeholder value####
RNAseq_GSE248274_2$Subset_Comparison <- "Unknown"

####Annotate the samples with correct sample groups####
RNAseq_GSE248274_2$Subset_Comparison <- ifelse(RNAseq_GSE248274_2$group == "Tx", "Treatment", "Control")

####Check the resulting converted data set####
print(head(RNAseq_GSE248274_2))
print(tail(RNAseq_GSE248274_2))
```

# Question 1

## Finding genes with the most and least significant effect on the continuous variable (dependent variable).
The code below takes *age* as the continuous variable and *gene expression* as the independent variable. Then a linear model is fitted for each gene, and the ANOVA test was done to extract the p-values.
```{r lm+ANOVA, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
#Create placeholder for p-values
gene_pvalues <- data.frame(SYMBOL = unique(RNAseq_GSE248274_2$SYMBOL), p_value = NA)

# Loop each gene then fit linear model and perform ANOVA test
gene_pvalues$p_value <- sapply(gene_pvalues$SYMBOL, function(gene) {
  # Subset the data for each gene
  gene_subset <- RNAseq_GSE248274_2[RNAseq_GSE248274_2$SYMBOL == gene, ]
  
  # Fit a linear model and perform ANOVA test
  anova_test <- anova(lm(age ~ Expression, data = gene_subset))
  
  # Return the p-value
  anova_test$`Pr(>F)`[1]
})

#Find genes with most and least significant effect
most_signif <- gene_pvalues[order(gene_pvalues$p_value), ]
least_signif <- gene_pvalues[order(gene_pvalues$p_value, decreasing = TRUE), ]

print(paste("Most significant gene based on ANOVA test is:", most_signif[1, "SYMBOL"], ", with p-value =", format(most_signif[1, "p_value"], scientific = TRUE)))
print(paste("Least significant gene based on ANOVA test is:", least_signif[1, "SYMBOL"], ", with p-value =", format(least_signif[1, "p_value"], scientific = TRUE)))
```

The code below generates the p-value distribution for gene effect on age for all of the genes.
```{r histogram, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

ggplot(gene_pvalues, aes(x = p_value)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of p-values for all genes", x = "p-value", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10))
```

The code below performs the linear model of the gene with the most and least significant effects.
```{r models, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
#Create models for most and least significant genes
most_signif_model <- lm(age ~ Expression, data = RNAseq_GSE248274_2[RNAseq_GSE248274_2$SYMBOL == most_signif[1, "SYMBOL"], ])
least_signif_model <- lm(age ~ Expression, data = RNAseq_GSE248274_2[RNAseq_GSE248274_2$SYMBOL == least_signif[1, "SYMBOL"], ])
```

This is the summary and ANOVA results for the linear models of the most significant gene
```{r summ1, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
summary(most_signif_model)
anova(most_signif_model)
```

This is the summary and ANOVA results for the linear models of the least significant gene
```{r summ2, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
summary(least_signif_model)
anova(least_signif_model)
```


The code below generates the diagnostic plots for the linear models of the gene with the most significant effects.
```{r plots1-1, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

par(mfrow = c(2, 2))
plot(most_signif_model, main = paste("Diagnostic Plots for", most_signif[1, "SYMBOL"], "Model"))
```

The code below generates the diagnostic plots for the linear models of the gene with the least significant effects.
```{r plots1-2, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
par(mfrow = c(2, 2))
plot(least_signif_model, main = paste("Diagnostic Plots for", least_signif[1, "SYMBOL"], "Model"))
```


# Question 2

## Batch correction prior to fitting the models.

Batch correction was performed by modifying the DEseq2 code to export the normalized values with the following code (hidden).
```{r normBatch, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
#Get data from GEO
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
RNAseqData2 <- GEODataDownload(DS = "GSE248274",
                              gpl = "GPL21697",
                              gsm = "01001001010101",
                              PlateAnnotInfo = "Only needed if the data was generated using MicroArrays",
                              GenerateMetaData = "Only needed if the data was generated using MicroArrays",
                              Technology = "RNAseq")
data.table(head(as.data.frame(RNAseqData2),10), filter = 'top', options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px", autoWidth = TRUE))

####Save the RNAseq data set to the working directory####
saveRDS(RNAseqData2, file = "~/RNAseq_GSE248274_norm.rds")

####Reload the file into R environment####
RNAseq_GSE248274_norm <- readRDS("~/RNAseq_GSE248274_norm.rds")
```

Add new columns called age, sex, group, and Subset_Comparison for sample annotation

```{r longData2, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
####Create a new data frame to group the samples by ctrl vs. Tx ####
new_df2 <- data.frame(
  Sample_ID = c("GSM7910769.y", "GSM7910771.y", "GSM7910772.y", "GSM7910773.y", "GSM7910774.y", "GSM7910775.y", "GSM7910776.y", "GSM7910777.y", "GSM7910778.y", "GSM7910780.y", "GSM7910782.y", "GSM7910783.y", "GSM7910785.y", "GSM7910786.y"),
  group = c("Ctrl", "Tx", "Ctrl", "Ctrl", "Tx", "Ctrl", "Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx"),
  age = c(58,58,55,55,55,64,64,64,71,71,88,88,68,68),
  sex = c("M","M","M","M","M","M","M","M","F","F","F","F","M","M")
)

#### Convert data set to long format to merge the sample group####
RNAseq_GSE248274_norm2 <- melt(RNAseq_GSE248274_norm, id.vars = c("ENTREZID", "SYMBOL", "GENENAME", "log2FoldChange", "pvalue", "padj"), 
                           variable.name = "Sample_ID", value.name = "Expression")

####Merge the sample group to the long format data set and append the new column #### 
RNAseq_GSE248274_norm2 <- merge(RNAseq_GSE248274_norm2, new_df2, by = "Sample_ID", all.x = TRUE)

####Append the new column and use a placeholder value####
RNAseq_GSE248274_norm2$Subset_Comparison <- "Unknown"

####Annotate the samples with correct sample groups####
RNAseq_GSE248274_norm2$Subset_Comparison <- ifelse(RNAseq_GSE248274_norm2$group == "Tx", "Treatment", "Control")

####Check the resulting converted data set####
print(head(RNAseq_GSE248274_norm2))
print(tail(RNAseq_GSE248274_norm2))
```


## Finding genes with the most and least significant effect on the continuous variable (dependent variable) within the normalized dataset.

The code below takes *age* as the continuous variable and *gene expression* as the independent variable. Then a linear model is fitted for each gene, and the ANOVA test was done to extract the p-values.
```{r lm+ANOVA2, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
#Create placeholder for p-values
gene_pvalues2 <- data.frame(SYMBOL = unique(RNAseq_GSE248274_norm2$SYMBOL), p_value = NA)

# Loop each gene then fit linear model and perform ANOVA test
gene_pvalues2$p_value <- sapply(gene_pvalues2$SYMBOL, function(gene) {
  # Subset the data for each gene
  gene_subset <- RNAseq_GSE248274_norm2[RNAseq_GSE248274_norm2$SYMBOL == gene, ]
  
  # Fit a linear model and perform ANOVA test
  anova_test <- anova(lm(age ~ Expression, data = gene_subset))
  
  # Return the p-value
  anova_test$`Pr(>F)`[1]
})

#Find genes with most and least significant effect
most_signif2 <- gene_pvalues2[order(gene_pvalues2$p_value), ]
least_signif2 <- gene_pvalues2[order(gene_pvalues2$p_value, decreasing = TRUE), ]

print(paste("Most significant gene in the normalized dataset based on ANOVA test is:", most_signif2[1, "SYMBOL"], ", with p-value =", format(most_signif2[1, "p_value"], scientific = TRUE)))
print(paste("Least significant gene in the normalized dataset based on ANOVA test is:", least_signif2[1, "SYMBOL"], ", with p-value =", format(least_signif2[1, "p_value"], scientific = TRUE)))
```

The code below generates the p-value distribution for gene effect on age for all of the genes in the normalized dataset.
```{r histogram2, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

ggplot(gene_pvalues2, aes(x = p_value)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of p-values for all genes in normalized data", x = "p-value", y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10))
```

The code below performs the linear models of the gene with the most and least significant effects in the normalized dataset.
```{r models2, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
#Create models for most and least significant genes
most_signif_model2 <- lm(age ~ Expression, data = RNAseq_GSE248274_norm2[RNAseq_GSE248274_norm2$SYMBOL == most_signif2[1, "SYMBOL"], ])
least_signif_model2 <- lm(age ~ Expression, data = RNAseq_GSE248274_norm2[RNAseq_GSE248274_norm2$SYMBOL == least_signif2[1, "SYMBOL"], ])
```

The code below present the model summary and ANOVA results for the linear models of the gene with the most significant effects in the normalized dataset.
```{r summ3, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
summary(most_signif_model2)
anova(most_signif_model2)
```

The code below present the model summary and ANOVA results for the linear models of the gene with the least significant effects in the normalized dataset.
```{r summ4, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
summary(least_signif_model2)
anova(least_signif_model2)
```


The code below generates the diagnostic plots for the linear models of the gene with the most significant effects.
```{r plots2-1, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

par(mfrow = c(2, 2))
plot(most_signif_model2, main = paste("Diagnostic Plots for", most_signif2[1, "SYMBOL"], "(Normalized)"))
```

The code below generates the diagnostic plots for the linear models of the gene with the least significant effects.
```{r plots2-2, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
par(mfrow = c(2, 2))
plot(least_signif_model2, main = paste("Diagnostic Plots for", least_signif2[1, "SYMBOL"], "(Normalized)"))
```

**Explanation of Results**

Histogram of PMSC4 (raw) data set clearly show a large variance in the frequency scale, ranging from about 50 to 600, thereby skewing the distribution to the left. Whereas PMS2P5 (normalized) data show a much closer frequency scale between 25 to 200, showing mostly equal distribution with a slight skew to the left.

Discussion on the most significant genes, PMSC4 (raw) and PMS2P5 (normalized) data:

Based on the residuals vs. fitted plots of PSMC4 (raw data) and PMS2P5 (normalized), both plots doesn't show a clear pattern which indicates non-linear relationships between the dependent and independent variables. In this case, the model was able to capture the linear relationship. However, data points #2 & 12 are pointed out in the PSMC4 9 (raw) data set, and data point #20 is shown to be the farthest in PMS2P5 (normalized) data set.

The Q-Q plots show that both PMSC4 (raw) and PMS2P5 (normalized) data look quite normally distributed. Although one can see that PMS2P5 (normalized) data have a couple data points (#20 & 26) that are a bit further off the straight line, which could be a potential problem. There are also outliers, #2 & 12, away from the straight line in the PMSC4 (raw) data set.

The Scale-LOcation plot checks for the assumption of equal variance (homoscedasticity). Both plots for PMSC4 (raw) and PMS2P5 (normalized) data show an equally random spread of data points indicating good, equal variance (for the most part). Although, we see again data point #20 in PMS2P5 (normalized) data has the widest spread. Whereas, #2 & 12 are seen again in this plot to have a wider spread in PMSC4 (raw) data set.

Seeing the outliers that have been consistent across the plots within each data set, the residuals vs. leverage plots helps to determine whether these data points are influential and does not follow the trend. It seems that data points #2 & 12, in PMSC4 (raw) data set, are influential to the regression results because they have a high Cook's distance score as they touch the dashed line. Within PMS2P5 (normalized) data set, data point #20 reach close to the Cook's distance line, also indicating that it might be influential to the regression. 

Therefore, for the most significant gene: PMSC4 (raw) data set, the four plots show potential problematic cases with two data points that is influential to the regression analysis and probably needs a closer look individually whether there is really a linear relationship between the dependent and independent variables. If so, perhaps *age* may not be the best variable to represent the model of the differential gene expression. In comparison, PMS2P5 (normalized) data set is the *better* data set because it improved upon the previous by having only one problematic data point that could be influential to the regression; However, the fact that four plots show the consistent, problematic data point may warrant further look into the data.    

Discussion on the least significant genes, PCSK9 (raw) and TSPAN1 (normalized) data:

In both cases, raw and normalized data does not show any distinctive pattern in the residual vs. fitted values plots. The Q-Q plots show close enough to a normal distribution in both data sets, however both also show 3 outliers each, away from the straight line that could be problematic. The Scale-Location for PCSK9 (raw) data does not show homoscedasticity, while TSPAN1 (normalized) data show better equal (randomly) spread points. Finally, PCSK9 (raw) plot of residuals vs. leverage show more data points with a high Cook's distance indicating the problematic outliers seen in previous plots could be influential to the linear regression. Whereas, TSPAN1 (normalized) plot does not show data points that is really close to Cook's distance, indicating any outliers observed in previous plots may not be as influential/problematic to the linear regression. 
Based on these observations, the normalized data again proves to be the *better* dataset to create models to study the relationship between age and gene expression in this study.

# Question 3

## The relationship between the fitted and residual values of the most significant gene to the dependent variable.

Using the gene with the greatest association to the dependent variable (in this case, age), the code below obtains the associated data and creates a new data frame for plotting.
```{r getData, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Get the most significant gene's data
most_sig_gene <- most_signif2[1, "SYMBOL"]
gene_data <- RNAseq_GSE248274_norm2[RNAseq_GSE248274_norm2$SYMBOL == most_sig_gene, ]
# Clean up the data
gene_data <- gene_data[!is.na(gene_data$age), ]

# Fit the model
fit_model <- lm(age ~ Expression, data = gene_data)

# Get predicted (fitted) values
fitted_values <- predict(fit_model)
# Get observed values
obsvd_values <- gene_data$age
# Get residuals
resid_values <- resid(fit_model)

# Create a data frame for plotting
new_df3 <- data.frame(
  Expression = gene_data$Expression,
  Observed = obsvd_values,
  Fitted = fitted_values,
  Residuals = resid_values
)
head(new_df3)
```

The code below will create plots between the fitted data and the observed values, and between the fitted data and the residual values.

```{r plots3, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

# Plot 1: Observed vs Fitted values with regression line
plot1 <- ggplot(new_df3, aes(x = Expression)) +
  geom_point(aes(y = Observed), color = "blue", size = 3, alpha = 0.6) +
  geom_line(aes(y = Fitted), color = "magenta", linewidth = 1) +
  labs(title = paste("Observed vs Fitted Values for", most_sig_gene),
       x = "Expression Level",
       y = "Age") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10, face = "bold"))

# Plot 2: Residuals vs Fitted values
plot2 <- ggplot(new_df3, aes(x = Fitted, y = Residuals)) +
  geom_point(color = "black", size = 3, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "magenta") +
  labs(title = paste("Residuals vs Fitted Values for", most_sig_gene),
       x = "Fitted Values",
       y = "Residuals") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10, face = "bold"))

# Arrange plots side by side
#dev.off() #optional in case of invalid graphics state
plot_grid(plot1, plot2, ncol = 2)

```

The following code will print the summary statistics of the Observed vs. Expression plot.
```{r summary3, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Print summary statistics
cat("Summary Statistics for most significant gene,", most_sig_gene, ":\n", "Correlation between Observed and Fitted values =", cor(obsvd_values, fitted_values), "\n", "Mean Absolute Error =", mean(abs(resid_values)), "\n", "Root Mean Square Error =", sqrt(mean(resid_values^2)), "\n")

```