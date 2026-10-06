---
title: "RBIF111 - Homework 7"
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
library(limma)
library(pscl)
```
  
  
  
# Question 1  
## Cross validation  
  
Cross-validation is a technique used in statistical analysis and machine learning to evaluate how well a model will estimate the average prediction error of models fit on other unseen training sets drawn from the same population.[1-3] The method involves partitioning the original data set randomly into subsets, and then using one of the subsets for testing and the remaining subsets for training the model.[2,3] The process is repeated *x* times, with each of the subsets used once for testing and the rest for training, and the resulting estimates of the model's performance are averaged to give an overall measure of the model's quality.[2,3] The goal of cross-validation is to estimate the expected level of performance of a model on a new, independent data, and to avoid overfitting -- where a model performs well on the training data but poorly on new data.[2,3]  
  
### Partitioning data for cross validation  
In cross-validation, the original data set is partitioned randomly into *k* subsets, or "folds", where *k* is a predetermined number usually chosen between 5 and 10.[2,3] Each fold is used once for testing the model performance, and the remaining *k-1* folds are used for training. This process is repeated *k* times, with each fold used once for testing and the other *k-1* folds used for training. The results of each iteration are then averaged to give an overall estimate of the model's performance.[2,3]  
  
### Potential issues  
1. If the training or test data is not representative of the overall data set, it can introduce bias into the model.[2]  
2. If the model is too complex, it may fit the training data too closely and not generalize well to new data, i.e., overfitting.[2,3]  
3. If the model is too simple, it may not capture the relevant features of the data and fail to generalize well to new data, i.e., underfitting.[2,3]  
4. If the data set is very large, cross-validation can be computationally expensive and time-consuming.[2]  
5. A larger number of folds typically provides a more reliable estimate of the model's performance, but it is also computationally expensive.[3]  
6. In time series or sequential data, it may not be appropriate to randomly partition the data. Instead, a window or sliding approach may be more appropriate to capture the relevant time aspects in the data.[2,3]  
  
  
**References**  
[1] Hawkins, D.M., Basak, S.C. and Mills, D., 2003. Assessing model fit by cross-validation. *Journal of chemical information and computer sciences*, 43(2), pp.579-586.  
[2] Krstajic, D., Buturovic, L.J., Leahy, D.E. and Thomas, S., 2014. Cross-validation pitfalls when selecting and assessing regression and classification models. *Journal of cheminformatics*, 6, pp.1-15.  
[3] Bates, S., Hastie, T. and Tibshirani, R., 2024. Cross-validation: what does it estimate and how well does it do it?. *Journal of the American Statistical Association*, 119(546), pp.1434-1445.  
  
  
  
# Question 2  
  
## Download data set from GEO  
The data set obtained from GEO is a study on decreased carnitine status with aging that may contribute to frailty and fatigue. It is a cross-sectional study comparing *fit* and *frail* older people, in which the authors hypothesized that muscle carnitine deficiency is associated with pre-frailty, diminished physical performance, and altered mitochondrial function; therefore, muscle mitochondrial gene expression was analyzed. There are 80 samples total, with 26 young and healthy individuals serving as the control.  
  
The (hidden) code below will download the GEO data set.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='hide', cache=FALSE}
GEODataDownload <- function(DS, gpl, gsm, PlateAnnotInfo, GenerateMetaData, Technology){
  library(DESeq2); library(limma)
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
```
  
  
The code below will save and load the data set into the R environment.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# get RNAseq data
RNAseqData <- GEODataDownload(DS = "GSE144304",
                              gpl = "GPL18573",
                              gsm = "11111111111111111111111111111111111111111111111111111100000000000000000000000000",
                              PlateAnnotInfo = "Only needed if the data was generated using MicroArrays",
                              GenerateMetaData = "Only needed if the data was generated using MicroArrays",
                              Technology = "RNAseq")
#data.table(head(as.data.frame(RNAseqData),10), filter = 'top', options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px", autoWidth = TRUE))

####Save the RNAseq data set to the working directory####
saveRDS(RNAseqData, file = "~/RNAseq_GSE144304_norm.rds")

####Reload the file into R environment####
GSE144304 <- readRDS("~/RNAseq_GSE144304_norm.rds")

```
  
  
The column below will add columns for the continuous variables associated with data set, obtained from the metadata.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
##Edit column names that end with ".y"
colnames(GSE144304) <- gsub(".y$", "", colnames(GSE144304))

##clean up column names
remove_cols <- grepl(".x$", colnames(GSE144304))
GSE144304_2 <- subset(GSE144304, select = !remove_cols)

head(GSE144304_2)

# Convert data set to long format to merge with sample annotation
GSE144304_3 <- reshape2::melt(GSE144304_2, id.vars = c("ENTREZID", "SYMBOL", "GENENAME", "log2FoldChange", "pvalue", "padj"), variable.name = "Sample_ID", value.name = "Expression")

# Add in metadata from folder
metadata <- read.csv("./SRA_GSE144304.csv")
metadata <- metadata %>% mutate(outcome = as.factor(outcome))

head(metadata)

# Merge the data frames by Sample_ID
GSE144304_4 <- merge(GSE144304_3, metadata, by = "Sample_ID", all.x = TRUE)

# Check the resulting converted data set
print(head(GSE144304_4))
print(tail(GSE144304_4))

```
  
  
# Question 2  
## Determine the feature with the greatest prediction value for the outcome  
  
The codes below will perform logistic regression analysis to determine which continuous feature has the greatest predictive value for the outcome. This feature will be analyzed further in question #3.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Figure out the feature with the greatest predictive value for the outcome
set.seed(123)

# set up the model
my_model <- glm(outcome ~ Expression + age + weight + BMI, data = GSE144304_4, family = binomial(link = "logit"))
# check the model
summary(my_model)

# analyze the model
my_anova <- anova(my_model, test = "Chisq")
print(my_anova)

# assess the model fit
my_fit <- pR2(my_model)
print(my_fit)
```
  
  
**Discussion**  
The result of the model show that Expression is not a significant feature, but features: age, weight, and BMI, have the lowest p-values and therefore more significant. The result of anova of the model show that features: age and BMI, have the highest deviance; which means both these features significantly reduces the residual deviance, with the weight feature less so. Because the control group in this data set are young and healthy people, age is obviously a major factor and will not be selected as the most significant feature here. Instead, *BMI* will be selected as the feature with the greatest predictive value for the outcome (whether the subject will be fit or frail in old age).
The McFadden R-squared index may be used as an equivalent to the R-squared of linear regression, to assess the model fit, which resulted in a value of ~0.312.  
  
  
  
# Question 3  
  
## Cross-validation and bootstrap  
  
  
The code below is function to perform the cross-validation, with folds = 5.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Cross-validation function
xvalidation <- function(train_data, test_data, n_folds = 5) {
  # Extract unique Sample_IDs and genes
  train_ids <- unique(train_data$Sample_ID)
  unique_genes <- unique(train_data$SYMBOL)
  
  # Create cross-validation groups
  xval_groups <- sample(1:n_folds, size = length(train_ids), replace = TRUE)
  names(xval_groups) <- train_ids
  
  # Function to evaluate genes for each fold
  run_fold <- function(fold) {
    # Split into training and test sets for cross-validation
    cv_test_ids <- names(xval_groups[xval_groups == fold])
    cv_train_ids <- names(xval_groups[xval_groups != fold])
    
    train_cv <- train_data %>% filter(Sample_ID %in% cv_train_ids)
    test_cv <- train_data %>% filter(Sample_ID %in% cv_test_ids)
    
    # Function to evaluate each gene
    evaluate_gene <- function(gene) {
      train_gene <- train_cv %>% filter(SYMBOL == gene)
      test_gene <- test_cv %>% filter(SYMBOL == gene)
      
      if(nrow(train_gene) < 2 || nrow(test_gene) < 1) {
        return(Inf)
      }
      
      tryCatch({
        model <- lm(BMI ~ Expression, data = train_gene)
        pred <- predict(model, newdata = test_gene)
        mean((pred - test_gene$BMI)^2, na.rm = TRUE)
      }, error = function(e) Inf)
    }
    
    # Evaluate all genes
    gene_mse <- sapply(unique_genes, evaluate_gene)
    names(which.min(gene_mse))
  }
  
  # Run cross-validation for all folds
  fold_best_genes <- sapply(1:n_folds, run_fold)
  
  # Get best gene and compute its MSE
  best_gene <- names(sort(table(fold_best_genes), decreasing = TRUE)[1])
  
  final_model <- lm(BMI ~ Expression, 
                   data = train_data %>% filter(SYMBOL == best_gene))
  final_pred <- predict(final_model, 
                       newdata = test_data %>% filter(SYMBOL == best_gene))
  mse <- mean((final_pred - test_data$BMI)^2, na.rm = TRUE)
  
  list(best_genes = fold_best_genes, best_gene = best_gene, mse = mse)
}

```
  
  
The code below is a function that will perform the bootstrap, at least n = 50, and up to n = 250.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Bootstrap function
bootstrap_genes <- function(long_data, n_boots, n_folds) {
  sample_ids <- unique(long_data$Sample_ID)
  
  # Function to run one bootstrap iteration
  run_bootstrap <- function(boot_iter) {
    # Bootstrap: sample with replacement
    train_ids <- sample(sample_ids, size = length(sample_ids), replace = TRUE)
    test_ids <- setdiff(sample_ids, train_ids)
    
    train_data <- long_data %>% filter(Sample_ID %in% train_ids)
    test_data <- long_data %>% filter(Sample_ID %in% test_ids)
    
    # Execute cross-validation function on this bootstrap sample
    xvalidation(train_data, test_data, n_folds)
  }
  
  # Run all bootstrap iterations
  results <- lapply(1:n_boots, run_bootstrap)
  
  # Extract results
  all_best_genes <- unlist(lapply(results, function(x) x$best_genes))
  mse_list <- unlist(lapply(results, function(x) x$mse))
  
  # Get top 5 genes
  top_genes <- names(sort(table(all_best_genes), decreasing = TRUE)[1:5])
  
  list(top_genes = top_genes, mse = mse_list)
}

```
  
  
The code below will execute the functions starting with bootstrapping n = 50, then higher at n = 250.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Execute the functions
set.seed(123)  # For reproducibility
result_10 <- bootstrap_genes(GSE144304_4, n_boots = 10, n_folds = 5)
result_50 <- bootstrap_genes(GSE144304_4, n_boots = 50, n_folds = 5)
#result_250 <- bootstrap_genes(GSE144304_4, n_boots = 250)

# Print results
cat("Top 5 genes from n=10 bootstraps:\n")
print(result_10$top_genes)
cat("Top 5 genes from n=50 bootstraps:\n")
print(result_50$top_genes)
```
  
  
The code below will create a boxplot to visualize the comparison of the top most frequently identified genes, and the mean squared error estimates.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Prepare data for plotting
mse_data <- data.frame(
  MSE = c(result_10$mse, result_10$mse),
  Simulation = factor(c(rep("10 boots", length(result_10$mse)),
                       rep("50 boots", length(result_50$mse)))))

# Create boxplot
ggplot(mse_data, aes(x = Simulation, y = MSE)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  labs(title = "Comparison of MSE Distribution",
       x = "Number of Bootstrap Simulations",
       y = "Mean Squared Error")

```
  
  
**Discussion**  
In the boxplot for both 10 and 50 bootstrap simulations, there is a clear outlier above the upper whisker (around MSE of 25). Outliers can distort the mean and standard deviation, leading to skewed results. In this case, the outlier raises the overall spread of the MSE distribution and may indicate a specific scenario where the model performed poorly. If the outliers are rare but valid cases then they should be investigated; but if they are due to data error or incomplete data pre-processing, then they could be removed to improve the reliability.  
Sample normalization adjusts the data to ensure comparability and uniform scaling, especially when dealing with gene expression data that can have high variability. Fortunately, the data set has been normalized from the beginning of the analysis. If the data set was not normalized, differences in sample scaling can lead to inconsistent results, in which the model may weight certain samples or genes incorrectly. The mean squared error can also be overestimated in some simulations, leading to more variability in the results.  
The boxplot show that both 10 and 50 simulations have similar median mean squared error values, showing a consistent model performance. The spread (interquartile range) and whiskers are slightly wider for the n=10 bootstrap simulation, which could result from higher variability with the fewer simulation. Increasing the bootstrap simulation n>50 can help reduce the variability and provide more robust estimates of the mean squared error distributions.  
