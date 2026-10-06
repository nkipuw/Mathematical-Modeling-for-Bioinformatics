---
title: "RBIF111 - Homework 6"
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
```


```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}

```
  
  
# Download Data from GEO  
  
The downloaded data set is a study to evaluate the changes in molecular expression of circulating tumor cells (CTC) in patients with hepatocellular carcinoma (HCC) treated with Atezolizumab & Bevacizumab. Any changes in molecular expression and CTC numbers were analyzed to identify effective biomarkers. The RNAseq data set contain 234 genes, 10 samples for control, and 10 samples for treatment.  
  
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
RNAseqData <- GEODataDownload(DS = "GSE261186",
                              gpl = "GPL15520",
                              gsm = "01010101010101010101",
                              PlateAnnotInfo = "Only needed if the data was generated using MicroArrays",
                              GenerateMetaData = "Only needed if the data was generated using MicroArrays",
                              Technology = "RNAseq")

####Save the RNAseq data set to the working directory####
saveRDS(RNAseqData, file = "~/RNAseq_GSE261186_norm.rds")

####Reload the file into R environment####
GSE261186 <- readRDS("~/RNAseq_GSE261186_norm.rds")

##Edit column names that end with ".y"
colnames(GSE261186) <- gsub(".y$", "", colnames(GSE261186))

##clean up column names
remove_cols2 <- grepl(".x$", colnames(GSE261186))
GSE261186_2 <- subset(GSE261186, select = !remove_cols2)

head(GSE261186_2)
```

  
The code below will add new columns: age and group, as sample annotations.
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Create a new data frame to group the samples by ctrl vs. Diabetic
new_df2 <- data.frame(
  Sample_ID = colnames(GSE261186_2)[7:ncol(GSE261186_2)],
  age = c(79.9, 79.9, 74.7, 74.7, 43.8, 43.8, 80.7, 80.7, 70.2, 70.2, 87.7, 87.7, 64.1, 64.1, 62, 62, 69.4, 69.4, 40.6, 40.6),
  group = c("Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx", "Ctrl", "Tx"),
  stringsAsFactors = FALSE
)

# Convert data set to long format to merge the sample group
GSE261186_3 <- reshape2::melt(GSE261186_2, id.vars = c("ENTREZID", "SYMBOL", "GENENAME", "log2FoldChange", "pvalue", "padj"), variable.name = "Sample_ID", value.name = "Expression")

# Merge the sample group to the long format data set and append the new column 
GSE261186_4 <- merge(GSE261186_3, new_df2, by = "Sample_ID", all.x = TRUE)

# Check the resulting converted data set
print(head(GSE261186_4))
print(tail(GSE261186_4))

```
  
  
# Question 1  
  
## Build linear model, calculate p-value with ANOVA, and perform Bonferroni & FDR correction  
  
  
The code below will build a linear model using ANOVA and obtain p-value for each gene. 

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
## fit linear model with ANOVA
lin_mod <- GSE261186_4 %>%
    group_by(SYMBOL) %>%  ## group by gene
    summarise(
        ANOVA_pvalue = {
            fit <- aov(Expression ~ age, data = pick(everything()))  ## fit ANOVA model
            summary(fit)[[1]]$`Pr(>F)`[1]   ## get ANOVA p-value
        }
    )
head(lin_mod)

## perform Bonferroni & FDR corrections
lin_mod2 <- lin_mod %>%
    mutate(
        Bonferroni_adj_pval = p.adjust(ANOVA_pvalue, method = "bonferroni"),
        FDR_adj_pval = p.adjust(ANOVA_pvalue, method = "fdr")
    )
head(lin_mod2)

## merge all statistical results
GSE261186_5 <- merge(GSE261186_4, lin_mod2, by = "SYMBOL", all.x = TRUE)
head(GSE261186_5)

```
  
  
The code below will generate a scatter plot of significant genes for each p-value measurement and the -log10(p-value) on the y-axis.  
  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Add -log10(p-value) columns for graphs
GSE261186_6 <- GSE261186_5 %>%
  mutate(
    log_ANOVApval = -log10(ANOVA_pvalue),
    log_Bonferroni = -log10(Bonferroni_adj_pval),
    log_FDR = -log10(FDR_adj_pval)
  )

## Create scatter plots for each p-value type
# ANOVA p-value plot
ggplot(GSE261186_6, aes(x = reorder(SYMBOL, -log_ANOVApval), y = log_ANOVApval)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Genes", y = "-log10(ANOVA p-value)", 
       title = "Gene Significance by ANOVA p-value", hjust = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "magenta")

# Bonferroni adjusted p-value plot
ggplot(GSE261186_6, aes(x = reorder(SYMBOL, -log_Bonferroni), y = log_Bonferroni)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Genes", y = "-log10(Bonferroni adjusted p-value)",
       title = "Gene Significance by Bonferroni Adjusted p-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange")

# FDR adjusted p-value plot
ggplot(GSE261186_6, aes(x = reorder(SYMBOL, -log_FDR), y = log_FDR)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Genes", y = "-log10(FDR adjusted p-value)",
       title = "Gene Significance by FDR Adjusted p-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

```
  
  
**Discussion**  
Bonferroni correction is conservative, especially with a large number of tests, therefore many genes that were significant under the original p-values will lose significance. The FDR correction is less stringent than Bonferroni and allows more genes to remain significant. The gene order under FDR is less likely to change dramatically compared to Bonferroni, as FDR focuses on controlling the false discovery rate rather than individual test significance. Under Bonferroni correction, genes with moderate (original) p-values may become insignificant, leading to a reordered ranking with a smaller set of significant genes. The FDR correction maintains more significant genes, and the overall rank of genes is similar to the original p-values, but a few genes with borderline original p-values may shift.  
  
  
  
# Question 2  
  
  
The code below will calculate the average level of expression for each gene for the two groups: "Tx" for treatment and "Ctrl" for control.  
  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# calculate the average expression of each gene for the Tx group
avgExp_Tx <- GSE261186_4 %>%
  filter(group == "Tx") %>%
  group_by(SYMBOL) %>%
  summarise(Tx_avgExpr = mean(Expression))

# calculate the average expression of each gene for the Ctrl group
avgExp_Ctrl <- GSE261186_4 %>%
  filter(group == "Ctrl") %>%
  group_by(SYMBOL) %>%
  summarise(Ctrl_avgExpr = mean(Expression))

# merge the two data frames
avgExp <- merge(avgExp_Tx, avgExp_Ctrl, by = "SYMBOL")

# calculate the ad-hocfold change
avgExp$adhoc_FC <- avgExp$Tx_avgExpr / avgExp$Ctrl_avgExpr

head(avgExp)
```
  
  
The code below will perform non-parametric Wilcoxon test for each gene. Then the Bonferroni and FDR correction on the p-values will be calculated.
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# perform nonparametric significance test for each gene
wilcox_test <- GSE261186_4 %>%
  group_by(SYMBOL) %>%
  summarise(
    adhoc_pvalue = wilcox.test(Expression[group == "Tx"], Expression[group == "Ctrl"])$p.value, exact = FALSE)

# calculate FDR & Bonferroni correction
wilcox_test$adhoc_FDR <- p.adjust(wilcox_test$adhoc_pvalue, method = "fdr")
wilcox_test$adhoc_Bonferroni <- p.adjust(wilcox_test$adhoc_pvalue, method = "bonferroni")

# clean up data frame
wilcox_test <- subset(wilcox_test, select = c("SYMBOL", "adhoc_pvalue", "adhoc_FDR", "adhoc_Bonferroni"))

head(wilcox_test)
```
  
  
The code below will merge the generated statistical results to the original downloaded data.
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# merge the two data frames
merge_stats <- merge(avgExp, wilcox_test, by = "SYMBOL")

# merge with original downloaded data
GSE261186_7 <- merge(GSE261186_2, merge_stats, by = "SYMBOL")
GSE261186_7 <- GSE261186_7 %>%
  relocate(c(Tx_avgExpr, Ctrl_avgExpr, adhoc_FC), .after = padj) %>%
  relocate(c(adhoc_pvalue, adhoc_FDR, adhoc_Bonferroni), .after = adhoc_FC)

head(GSE261186_7)

```
  
  
**Discussion**  
The results from Bonferroni correction shows that all of the adjusted p-values = 1, while the adjusted p-values from FDR vary. The Bonferroni correction adjusts p-values by multiplying them by the total number of tests. It is highly conservative, especially when performing many tests because it aims to strictly control the family-wise error rate. The FDR method is less strict and focuses on controlling the expected proportion of false positives; therefore, it allows for more significant findings in datasets with numerous tests compared to Bonferroni. The raw p-values show that the majority have values >0.05. With 234 genes in the data set multiplied by the raw p-values, the Bonferroni correction pushes all adjusted p-value above 1, which gets truncated to 1. Having a sample size of 10 each for control and treatment group might be too small and the Wilcoxon test may not have enough power to detect significant differences, leading to higher raw p-values. The FDR adjusted p-values vary because the method is less strict and adapts based on the distribution of raw p-values and their ranks. It is better suited for cases such as this, gene expression studies, where some false positives are tolerable.  
  
  
  
# Question 3  
  
  
The code below will find the most & least significant genes, as well as obtain the expression values for each of them.  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
dev.off()
# find the most & least significant genes
most_sig <- GSE261186_2 %>%
  filter(padj == min(padj))
print(most_sig)

least_sig <- GSE261186_2 %>%
  filter(padj == max(padj))
print(least_sig)

# get expression values for most & least significant gene
most_sig_expr <- GSE261186_4 %>%
  filter(SYMBOL == most_sig$SYMBOL) %>%
  select(Expression, age)
head(most_sig_expr)

least_sig_expr <- GSE261186_4 %>%
  filter(SYMBOL == least_sig$SYMBOL) %>%
  select(Expression, age)
head(least_sig_expr)

```
  
  
The code below will perform the permutation test to obtain the permuted p-value slopes and sum of squares.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# get original linear model fit
orig_fit <- lm(Expression ~ age, data = GSE261186_4)
orig_slope <- coef(orig_fit)[2]  ##get slope
orig_ANOVA <- anova(orig_fit) # nolint
orig_sumsq <- orig_ANOVA$"Sum Sq"[1]  ##get sum sq for age

## Function to perform permutation test
perm_test <- function(expr_data, n_perm = 1000){
  # store permuted values
  perm_slopes <- numeric(n_perm)
  perm_sumsq <- numeric(n_perm)

  # perform permutation test
  for (i in 1:n_perm) {
    # permutation of the age variable
    perm_data <- expr_data
    perm_data$age <- sample(expr_data$age)

    # fit linear model to permuted data
    perm_fit <- lm(Expression ~ age, data = perm_data)
    perm_slopes[i] <- coef(perm_fit)[2]  ##get slope
    perm_sumsq[i] <- anova(perm_fit)$"Sum Sq"[1]  ##get sum sq for age
  }

  # calculate proportions
  pval_slope <- mean(abs(perm_slopes) >= abs(orig_slope))
  pval_sumsq <- mean(perm_sumsq <= orig_sumsq)

  # return results
  return(list(
    pval_slope = pval_slope,
    pval_sumsq = pval_sumsq,
    permuted_slopes = perm_slopes,
    permuted_sumsq = perm_sumsq))
}

set.seed(123)
# execute function
most_sig_perm <- perm_test(most_sig_expr)
least_sig_perm <- perm_test(least_sig_expr)

```
  
  
The code below will visualize the distribution of the simulated slopes and sum of squares generated with permutation using histograms.
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
## visualize the results
par(mfrow = c(2,2))
# histogram for most significant gene
hist(most_sig_perm$permuted_slopes,
     main = paste("Permuted Slopes Distribution for\n most significant gene: CD28"),
     xlab = "Slopes",
     col = "skyblue",
     xlim = c(-200, 30))
abline(v = orig_slope, col = "red", lwd = 1)

hist(most_sig_perm$permuted_sumsq,
     main = paste("Permuted Sum of Squares Distribution for\n most significant gene: CD28"),
     xlab = "Sum of Squares",
     col = "#f0a00c")

hist(least_sig_perm$permuted_slopes,
     main = paste("Permuted Slopes Distribution for\n least significant gene: NOTCH1"),
     xlab = "Slopes",
     col = "#e4bff5",
     xlim = c(-200, 100))
abline(v = orig_slope, col = "red", lwd = 1)

hist(least_sig_perm$permuted_sumsq,
     main = paste("Permuted Sum of Squares Distribution for\n least significant gene: NOTCH1"),
     xlab = "Sum of Squares",
     col = "#e2c8af")

```
  
  
**Discussion**  
Based on the original linear model of the whole data set, the most significant genes are at the top of the graph showing largest -log10(p-values), which are the smallest raw p-values, and a significant decline after the top few genes indicating a subset of genes has a strong association between age and gene expression. Permutation testing provides confirmation of the ANOVA results that the observed significance of CD28 gene is unlikely due to chance, and that NOTCH1 gene remain consistent with the null hypothesis.  
  
The slope distribution for the most significant gene, CD28, show that it is centered around zero. The observed slope from the original linear model is far outside the range of this distribution, indicating that it is unlikely to occur by chance. Thus, it supports the significance of CD28 gene's relationship with age as identified by ANOVA.  
The sum of squares distribution show that it is heavily skewed to the left towards small values. The observed sum of squares is much larger than most of the permuted values, thereby confirming that the variance explained by the model for CD28 is highly significant.  
  
The slope distribution for the least significant gene, NOTCH1, is also centered around zero, but the observed slope is closer near the center suggesting that it is not significantly different from random noise. This aligns with the result from ANOVA that identified NOTCH1 as having the least significant association with age.  
The sum of squares distribution also show heavy skew to the left towards small values.  
  
For the most significant gene, CD28, both ANOVA and permutation show strong relationship between gene expression and age. For the least significant gene, NOTCH1, both methods show no meaningful association. The consistent results suggest that the assumptions of the linear model are likely valid for this data set.