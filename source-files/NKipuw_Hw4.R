---
title: "RBIF111 - Homework 4"
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

# Question 1

The data set from UCI machine learning repository consists of a study on the donor database of Blood Transfusion Service Center in Hsin-Chu City in Taiwan. The center passes their blood transfusion service bus to one university in Hsin-Chu City to gather blood donated about every three months. The data contains 748 random donors from the donor database. The data features include continuous variables of: months since last donation, frequency (total number of donation event), total blood donated (in c.c.), months since first donation, and a classification variable representing whether donor donated blood in March 2007 (1 = donating blood; 0 = not donating blood).  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
getwd()

# Read in the data as CSV and name all the columns 
blood_dataset <- read.csv("transfusion.data", header = TRUE)
names(blood_dataset) <- c("last_donation", "frequency", "total_donation", "first_donation", "outcome")
head(blood_dataset)

##Partition data set based on hospital outcome
partition_data <- function(data){
    #Partition the data set into two based on outcome
    partition_data <- split(data, data$outcome)
    
    # Convert the list to a data frame
    df_0 <- partition_data[[1]]
    df_1 <- partition_data[[2]]
    
    # Add a column to identify the group
    df_0$group <- "No blood donation"
    df_1$group <- "Donate blood"
    
    # Combine the two data frames
    result_df <- rbind(df_0, df_1)
    
    return(result_df)
}
# Execute the function
blood_data <- partition_data(blood_dataset)
head(blood_data)

```

Check whether the continuous data variables are normally distributed
```{r, echo=FALSE}
par(mfcol = c(2,2))
hist(blood_data$last_donation)
hist(blood_data$frequency)
hist(blood_data$total_donation)
hist(blood_data$first_donation)
#dev.off()

```

It seems none of the continuous data variables are normally distributed, therefore the *Wilcoxon test* will be used instead of the *t-test* for non-parametric data.  
  
  
## Random sampling and Wilcoxon test
The code below takes a random sample of N=30 measurement from continuous variable columns from the subset data that were partitioned by outcome. Then the mean difference is calculated and a 95% confidence interval is constructed for the difference using Wilcoxon test. This is done 1,000 times and all the confidence intervals are recorded.  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source="fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Function to randomly sample and do non-parametric test
rsample_wcox <- function(data, column_name){
    #split data by group
    nondonor_grp <- data[data$group == "No blood donation", ]
    donor_grp <- data[data$group == "Donate blood", ]

    #take random sample
    nondonor_sample <- sample(nondonor_grp[[column_name]], size = 30, replace = FALSE)
    donor_sample <- sample(donor_grp[[column_name]], size = 30, replace = FALSE)

    #perform wilcoxon test
    wcox_result <- wilcox.test(nondonor_sample, donor_sample, conf.int = TRUE, conf.level = 0.95, exact = FALSE) #setting 'exact=FALSE' will eliminate ties#

    #Return mean difference and CI
    return(c(wcox_result$conf.int, mean_diff = mean(nondonor_sample) - mean(donor_sample)))
}
set.seed(123)
# Execute the function
blood_data1 <- replicate(1000, rsample_wcox(blood_data, "last_donation"))
blood_data2 <- replicate(1000, rsample_wcox(blood_data, "frequency"))
blood_data3 <- replicate(1000, rsample_wcox(blood_data, "total_donation"))
blood_data4 <- replicate(1000, rsample_wcox(blood_data, "first_donation"))

# Get CI & mean results for each variable
lower_CI1 <- blood_data1[1, ]
upper_CI1 <- blood_data1[2, ]
mean_diff1 <- blood_data1[3, ]

lower_CI2 <- blood_data2[1, ]
upper_CI2 <- blood_data2[2, ]
mean_diff2 <- blood_data2[3, ]

lower_CI3 <- blood_data3[1, ]
upper_CI3 <- blood_data3[2, ]
mean_diff3 <- blood_data3[3, ]

lower_CI4 <- blood_data4[1, ]
upper_CI4 <- blood_data4[2, ]
mean_diff4 <- blood_data4[3, ]

```

## Results from random sampling and Wilcoxon test

### a) The average lengths of the confidence interval for each variable:

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
avg_CI4_len <- mean(upper_CI4 - lower_CI4)
avg_CI3_len <- mean(upper_CI3 - lower_CI3)
avg_CI2_len <- mean(upper_CI2 - lower_CI2)
avg_CI1_len <- mean(upper_CI1 - lower_CI1)

print(paste("The average length of the confidence interval for the last donation is =", format(round(avg_CI1_len, 3), nsmall = 3)))

print(paste("The average length of the confidence interval for donation frequency is =", format(round(avg_CI2_len, 3), nsmall = 3)))

print(paste("The average length of the confidence interval for total blood donation is =", format(round(avg_CI3_len, 3), nsmall = 3)))

print(paste("The average length of the confidence interval for the time since first donation is =", format(round(avg_CI4_len, 3), nsmall = 3)))

```

### b) The population level difference (i.e., the difference between the means) for each variable:

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
pop_diff4 <- mean(blood_data$first_donation[blood_data$group == "Donate blood"]) - mean(blood_data$first_donation[blood_data$group == "No blood donation"])
pop_diff3 <- mean(blood_data$total_donation[blood_data$group == "Donate blood"]) - mean(blood_data$total_donation[blood_data$group == "No blood donation"])
pop_diff2 <- mean(blood_data$frequency[blood_data$group == "Donate blood"]) - mean(blood_data$frequency[blood_data$group == "No blood donation"])
pop_diff1 <- mean(blood_data$last_donation[blood_data$group == "Donate blood"]) - mean(blood_data$last_donation[blood_data$group == "No blood donation"])

print(paste("The population level difference for the last donation is =", format(round(pop_diff1, 3), nsmall = 3)))

print(paste("The population level difference for donation frequency is =", format(round(pop_diff2, 3), nsmall = 3)))

print(paste("The population level difference for total blood donation is =", format(round(pop_diff3, 3), nsmall = 3)))

print(paste("The population level difference for the time since first donation is =", format(round(pop_diff4, 3), nsmall = 3)))

```

### c) How often the confidence intervals contain the population level difference for each variable:

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
freq_poplevel4 <- mean(lower_CI4 < pop_diff4 & upper_CI4 > pop_diff4)

freq_poplevel3 <- mean(lower_CI3 < pop_diff3 & upper_CI3 > pop_diff3)

freq_poplevel2 <- mean(lower_CI2 < pop_diff2 & upper_CI2 > pop_diff2)

freq_poplevel1 <- mean(lower_CI1 < pop_diff1 & upper_CI1 > pop_diff1)

print(paste("The frequency that the CI contains the population level difference is =", format(round(freq_poplevel1, 3), nsmall = 3)))

print(paste("The frequency that the CI contains the population level difference is =", format(round(freq_poplevel2, 3), nsmall = 3)))

print(paste("The frequency that the CI contains the population level difference is =", format(round(freq_poplevel3, 3), nsmall = 3)))

print(paste("The frequency that the CI contains the population level difference is =", format(round(freq_poplevel4, 3), nsmall = 3)))
```


## Non-random sampling and Wilcoxon test analysis for the entire data for each variable

The code below will calculate the confidence interval for the full data set of each variable. Wilcoxon test will be used due to skewed distribution of the variables. The result is the lower and upper confidence interval values for each variable, along with the differences in mean.   

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
rsample_wcox2 <- function(data, column_name){
    #split data by group
    nondonor_grp <- data[data$group == "No blood donation", ]
    donor_grp <- data[data$group == "Donate blood", ]
    
    # get all the data per group
    nondonor_sample <- nondonor_grp[[column_name]]
    donor_sample <- donor_grp[[column_name]]

    #perform Wilcoxon test
    wcox_result <- wilcox.test(nondonor_sample, donor_sample, conf.int = TRUE, conf.level = 0.95, exact = FALSE) #setting 'exact = FALSE' will eliminate ties#

    #Return mean difference and CI
    return(c(wcox_result$conf.int, mean_diff = mean(nondonor_sample) - mean(donor_sample)))
}

set.seed(123)
# Execute the function for each variable
blood_data5 <- rsample_wcox2(blood_data, "last_donation")
blood_data6 <- rsample_wcox2(blood_data, "frequency")
blood_data7 <- rsample_wcox2(blood_data, "total_donation")
blood_data8 <- rsample_wcox2(blood_data, "first_donation")

# Get CI & mean results for full dataset
lower_CI5 <- blood_data5[1]
upper_CI5 <- blood_data5[2]
mean_diff5 <- blood_data5[3]

lower_CI6 <- blood_data6[1]
upper_CI6 <- blood_data6[2]
mean_diff6 <- blood_data6[3]

lower_CI7 <- blood_data7[1]
upper_CI7 <- blood_data7[2]
mean_diff7 <- blood_data7[3]

lower_CI8 <- blood_data8[1]
upper_CI8 <- blood_data8[2]
mean_diff8 <- blood_data8[3]

```
  
  
The code (hidden) below will visualize the comparison between the random sampling and the full data set for the age variable. *Note: code hidden due to extensive lines and may be collapsed for viewing*  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, fig.height=10, fig.width=8, class.source = "fold-hide", results='show', cache=FALSE}
# Visualize data
par(mfcol = c(2,2))
par(xpd=FALSE)

## Plot last blood donation
plot(1:1000, blood_data1[3,], ylim = range(c(blood_data1[1,], blood_data1[2,])),
     type = "n", xlab = "Iteration", ylab = "Confidence Interval",
     main = "Random Sample vs Full Dataset CI for Last Blood Donation")
# Add random sample CI
segments(1:1000, blood_data1[1,], 1:1000, blood_data1[2,], col = "skyblue", alpha = 0.1)
# Add full dataset CI as horizontal lines
abline(h = lower_CI5, col = "purple", lwd = 2, lty = 2)
abline(h = upper_CI5, col = "purple", lwd = 2, lty = 2)
# Add legend
par(xpd=TRUE)
legend("topright", inset = c(0, -0.05), legend = c("Random Sample CI", "Full Dataset CI"),
       col = c("skyblue", "purple"), lty = c(1, 2), lwd = c(1, 2))

## Plot frequency of blood donation
par(xpd=FALSE)
plot(1:1000, blood_data2[3,], ylim = range(c(blood_data2[1,], blood_data2[2,])),
     type = "n", xlab = "Iteration", ylab = "Confidence Interval",
     main = "Random Sample vs Full Dataset CI for Donation Frequency")
# Add random sample CI
segments(1:1000, blood_data2[1,], 1:1000, blood_data2[2,], col = "skyblue", alpha = 0.1)
# Add full dataset CI as horizontal lines
abline(h = lower_CI6, col = "magenta", lwd = 2, lty = 2)
abline(h = upper_CI6, col = "magenta", lwd = 2, lty = 2)
# Add legend
par(xpd=TRUE)
legend("topright", inset = c(0, -0.05), legend = c("Random Sample CI", "Full Dataset CI"),
       col = c("skyblue", "magenta"), lty = c(1, 2), lwd = c(1, 2))

## Plot total blood donation
par(xpd=FALSE)
plot(1:1000, blood_data3[3,], ylim = range(c(blood_data3[1,], blood_data3[2,])),
     type = "n", xlab = "Iteration", ylab = "Confidence Interval",
     main = "Random Sample vs Full Dataset CI for Total Donation")
# Add random sample CI
segments(1:1000, blood_data3[1,], 1:1000, blood_data3[2,], col = "skyblue", alpha = 0.1)
# Add full dataset CI as horizontal lines
abline(h = lower_CI7, col = "darkorange", lwd = 2, lty = 2)
abline(h = upper_CI7, col = "darkorange", lwd = 2, lty = 2)
# Add legend
par(xpd=TRUE)
legend("bottomright", inset = c(0,0), legend = c("Random Sample CI", "Full Dataset CI"),
       col = c("skyblue", "darkorange"), lty = c(1, 2), lwd = c(1, 2))

## Plot time since first blood donation
par(xpd=FALSE)
plot(1:1000, blood_data4[3,], ylim = range(c(blood_data4[1,], blood_data4[2,])),
     type = "n", xlab = "Iteration", ylab = "Confidence Interval",
     main = "Random Sample vs Full Dataset CI for First Donation")
# Add random sample CI
segments(1:1000, blood_data4[1,], 1:1000, blood_data4[2,], col = "skyblue", alpha = 0.1)
# Add full dataset CI as horizontal lines
abline(h = lower_CI8, col = "darkred", lwd = 2, lty = 2)
abline(h = upper_CI8, col = "darkred", lwd = 2, lty = 2)
# Add legend
par(xpd=TRUE)
legend("topright", inset = c(0, -0.05), legend = c("Random Sample CI", "Full Dataset CI"),
       col = c("skyblue", "darkred"), lty = c(1, 2), lwd = c(1, 2))

```
  
  
  
# Question 2  

The data obtained from GEO database is a study to identify the potential regulators in the progression of diabetic nephropathy (DN) using RNAseq of blood samples from healthy volunteers, patients with type-2 diabetes, and patients with DN. The data contains 14 total samples.

The code and results (hidden) below will download the data set from GEO database. *May collapse the code and results*
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
RNAseqData <- GEODataDownload(DS = "GSE154881",
                              gpl = "GPL24676",
                              gsm = "00000111111111",
                              PlateAnnotInfo = "Only needed if the data was generated using MicroArrays",
                              GenerateMetaData = "Only needed if the data was generated using MicroArrays",
                              Technology = "RNAseq")
data.table(head(as.data.frame(RNAseqData),10), filter = 'top', options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px", autoWidth = TRUE))

####Save the RNAseq data set to the working directory####
saveRDS(RNAseqData, file = "~/RNAseq_GSE154881_norm.rds")

####Reload the file into R environment####
GSE154881_norm <- readRDS("~/RNAseq_GSE154881_norm.rds")
```

The code below add new columns: group and subset comparison, for sample annotation.
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Create a new data frame to group the samples by ctrl vs. Diabetic
new_df <- data.frame(
  Sample_ID = c("GSM4681806.y", "GSM4681807.y", "GSM4681808.y", "GSM4681809.y", "GSM4681810.y", "GSM4681811.y", "GSM4681812.y", "GSM4681813.y", "GSM4681814.y", "GSM4681815.y", "GSM4681816.y", "GSM4681817.y", "GSM4681818.y", "GSM4681820.y"),
  age = c(28,35,49,62,57,54,48,65,68,47,49,54,69,52),
  group = c("Ctrl", "Ctrl", "Ctrl", "Ctrl", "Ctrl", "Diabetic2", "Diabetic2", "Diabetic2", "Diabetic2", "Diabetic2", "Diabetic1", "Diabetic1", "Diabetic1", "Diabetic1")
)

# Convert data set to long format to merge the sample group
GSE154881_2 <- melt(GSE154881_norm, id.vars = c("ENTREZID", "SYMBOL", "GENENAME", "log2FoldChange", "pvalue", "padj"), 
                           variable.name = "Sample_ID", value.name = "Expression")

# Merge the sample group to the long format data set and append the new column 
GSE154881_2 <- merge(GSE154881_2, new_df, by = "Sample_ID", all.x = TRUE)

# Append the new column and use a placeholder value
GSE154881_2$Subset_Comparison <- "Unknown"

# Annotate the samples with correct sample groups
GSE154881_2$Subset_Comparison <- with(GSE154881_2, 
    case_when(
        Sample_ID %in% c("GSM4681806.y", "GSM4681807.y", "GSM4681808.y", "GSM4681809.y", "GSM4681810.y") ~ "Healthy",
        Sample_ID %in% c("GSM4681811.y", "GSM4681812.y", "GSM4681813.y", "GSM4681814.y", "GSM4681815.y") ~ "Type2Diabetes",
        Sample_ID %in% c("GSM4681816.y", "GSM4681817.y", "GSM4681818.y", "GSM4681820.y") ~ "DiabeticNephropathy",
        TRUE ~ ""  # Default value for all other rows
    )
)

# Check the resulting converted data set
print(head(GSE154881_2))
print(tail(GSE154881_2))

```
  
  
## Build linear models  

The code below will build linear models and obtain p-value and confidence intervals of categorical variable effect on gene expression for each gene, using a for loop. 

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Build linear model data frame
lin_model <- data.frame(
    gene = unique(GSE154881_2$SYMBOL),
    lower_CI = NA,
    upper_CI = NA,
    pvalue = NA,
    estimate = NA
)
set.seed(123)
# Loop through each gene
# Linear model built on Expression as independent variable and age as dependent (continuous) variable
for (gene in lin_model$gene) {
    # get the data for the specific gene
    gene_data <- GSE154881_2[GSE154881_2$SYMBOL == gene, ]
    # fit the linear model
    lin_mod <- lm(Expression ~ age, data = gene_data)
    # get CI
    ci <- confint(lin_mod, level = 0.95)
    # get summary stats
    summary <- summary(lin_mod)
    # store results
    lin_model$lower_CI[lin_model$gene == gene] <- ci[2,1] # second row, first column
    lin_model$upper_CI[lin_model$gene == gene] <- ci[2,2] # second row, second column
    lin_model$pvalue[lin_model$gene == gene] <- summary$coefficients[2,4] # put pvalue in second row, fourth column
    lin_model$estimate[lin_model$gene == gene] <- summary$coefficients[2,1] # put size in second row, first column
}

head(lin_model)
```
  
The genes will be sorted in ascending order according to the p-value, and the confidence intervals will be visualized.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", results='show', cache=FALSE}
# Order genes by p-value
lm_ordered <- lin_model[order(lin_model$pvalue), ]

# Get top genes for visualization
num_genes <- 50  # Number of top genes to plot
top_genes <- head(lm_ordered, num_genes)
head(top_genes)

# Visualization using ggplot2
library(ggplot2)

ggplot(top_genes, aes(x = reorder(gene, -pvalue), y = estimate)) +
    geom_point(aes(color = pvalue < 0.05), size = 2) +
    geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "magenta") +
    coord_flip() +  # Flip coordinates for horizontal plot
    theme_minimal() +
    labs(title = "Top 50 Genes: Effect Sizes and Confidence Intervals",
         subtitle = "In ascending order by p-value",
         x = "Gene",
         y = "Effect Size",
         color = "Significant") +
    theme(
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "bottom"
    ) +
    scale_color_manual(values = c("black", "magenta"))

```
  
**Explanation of Plot**  
  
The linear model reveal whether gene expression (independent variable) differs among age (dependent variable). The resulting confidence interval shows the uncertainty in the effect size estimates, i.e., wider intervals means more uncertainty. If the interval crosses zero, then it is not significant. The error bars show the confidence intervals in the plot. The smaller the p-value, the more evidence against the null hypothesis; therefore, it was the variable used for sorting the (ascending) order of the genes. Of course this plot only covers the top 50 genes with p-vales < 0.05, when there are 138 genes that fit that criteria.  
  
  
## The least and most significant genes  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Print summary of significant genes
sig_genes <- sum(lin_model$pvalue < 0.05)
cat("The total number of significant genes (p < 0.05) =", sig_genes, ", out of 6888 unique genes.\n\n")


# Get the most and least significant genes
most_signif <- lin_model[which.min(lin_model$pvalue), ]
least_signif <- lin_model[which.max(lin_model$pvalue), ]

# Print the most and least significant genes
cat("The least significant gene is:", least_signif$gene, "with a p-value of =", round(least_signif$pvalue, 5), ",\n and the 95% confidence interval is =", round(least_signif$lower_CI, 3), "to", round(least_signif$upper_CI, 3), "\n\n")
cat("The most significant gene is:", most_signif$gene, "with a p-value of =", most_signif$pvalue,",\n and the 95% confidence interval is =", round(most_signif$lower_CI, 3), "to", round(most_signif$upper_CI, 3), "\n\n")

```
  
  
*Explanation of Results*  
  
The gene at the top, FLVCR1 is significant, because it has the smallest p-value (0.00158) and a tight confidence interval. The small p-value could be due to the very small variance in the confidence interval. It also shows a negative effect size but not a great magnitude, thereby indicating a negative (but small) relationship with the dependent variable.  If we look towards the middle of the plot, the gene AZU1 is also significant, because it also has a small p-value = 0.00935 (< $\alpha$ = 0.05) and tight confidence interval. However, it shows a positive effect size albeit a small magnitude. This indicates that this gene has a positive (but also small) relationship with the dependent variable.  
The gene at the bottom of the list, TBKBP1 has a p-value = 0.999 and a wide confidence interval. This indicates that the gene is not significant as the p-value is > 0.05 and the null hypothesis cannot be rejected.  

  

