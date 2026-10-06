---
title: "RBIF111 - Final Exam"
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
  

```{r library, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='hide', cache=FALSE}
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
library(survival)

```
  
    
      
        
# Question 1  
## Logistic regression  
### Input is ER IHC status (categorical variable) and output is the overall survival months (continuous variable)
Logistic regression is used to construct a linear model using "ER-IHC" as the input and "overall survival months" as the output, then the 95% confidence interval is calculated.  
  
The (hidden) code below will load the provided clinical data set and the mRNA expression data.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Load the data sets
clin_data <- read.csv("./TCGA_BRCA_clinical_data.csv")
mrna_data <- read.csv("./TCGA_BRCA_mRNA_expression.csv")
clin_data[1:3, 1:5]
mrna_data[1:3, 1:5]
```
  
  
```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
# Get relevant columns & omit NA values
clin_data2 <- clin_data %>%
  select(ER_Status_By_IHC, Overall_Survival_Months) %>%
  na.omit() %>%
  filter(ER_Status_By_IHC != "Indeterminate")

# Make sure ER_Status_By_IHC is a factor
clin_data2$ER_Status_By_IHC <- as.factor(clin_data2$ER_Status_By_IHC)

head(clin_data2)

# Logistic Regression on survival months and ER IHC status
log_reg <- glm(ER_Status_By_IHC ~ Overall_Survival_Months,
                family = binomial(link = "logit"),
                data = clin_data2)
# Get coefficients
summary(log_reg)
# Exponentiate the coefficient to get odds ratio
exp(coefficients(log_reg))
# Get confidence interval
confint(log_reg, level = 0.95)
# Exponentiate the CI
exp(confint(log_reg, level = 0.95))
```  
  
**Discussion**  
The logistic regression coefficients give the change in the log odds of the outcome for a one unit increase in the predictor variable. In the model summary results, a negative coefficient (-0.009) indicates that as the "Overall Survival Months" increases, the odds of ER IHC status being "Positive" decrease slightly, but the relationship is pretty weak since the coefficient is close to zero. Further, the odds ratio = 0.991 implies that for every 1-unit increase in the "Overall Survival Months", the odds of a positive ER IHC status decrease by (1 - 0.991 =) 0.009. The 95% confidence interval show a very small range for the effect of "Overall Survival Months," and since the odds ratio is so close to 1 it suggests that the continuous variable does not strongly influence the outcome. Finally, since the p-value = 0.3 (> 0.05), the null hypothesis is not rejected and the observed relationship between "overall Survival Months" and "ER IHC Status" is not statistically significant.  
  
  
  
# Question 2  
## Cox regression model (Survival Analysis)  
### Input is ER IHC status (categorical variable), outputs are the overall survival months (continuous variable) and overall survival status (categorical outcome)  
  
The code below will create a new data frame for analysis by Cox proportional hazard model.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Get relevant columns & omit NA values
clin_data3 <- clin_data %>%
  select(ER_Status_By_IHC, Overall_Survival_Months, Overall_Survival_Status) %>%
  na.omit() %>%
  filter(ER_Status_By_IHC != "Indeterminate")

# Convert survival status to binary (0 = alive, 1 = dead)
clin_data3$Overall_Survival_Status <- ifelse(clin_data3$Overall_Survival_Status == "1:DECEASED", 1, 0)

# Convert ER IHC status to a factor
clin_data3$ER_Status_By_IHC <- as.factor(clin_data3$ER_Status_By_IHC)

head(clin_data3)

```
  
    
The code below will first create a "Surv" object within the data frame, consisting of the overall survival months as the time variable and the overall survival outcome as the event. This object represents the survival information. Then the Cox proportional hazard model will be performed using the "Surv" object and the ER IHC status.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
# Create Surv object for model analysis
clin_data3.Surv <- transform(clin_data3,
                            surv = Surv(
                              time = clin_data3$Overall_Survival_Months,
                              event = clin_data3$Overall_Survival_Status,
                              type = "right"
                            )
)
# show the first 5 observations
clin_data3.Surv$surv[1:5]

## Cox proportional hazard model
# use the 'surv' column for the formula, not the whole data frame
cox_reg <- coxph(formula = surv ~ ER_Status_By_IHC, data = clin_data3.Surv) 
summary(cox_reg)

```
  
  
**Discussion**  
The p-value tests the null hypothesis that the coefficient is zero, meaning there is no association between ER IHC status and the outcome (i.e., death). The resulting p-value = 0.929, which is greater than the typical $\alpha$ = 0.05; Therefore, it fails to reject the null hypothesis and there is no significant association between ER IHC status and survival time.  
The calculated lower 95% = 0.463 and the upper 95% = 2.326. The range is quite wide, as it reaches below **1** (which indicates reduced risk) and above **1** (which indicates increased risk) -- where **1** indicates no effect. This wide range indicates that there is a lot of uncertainty in the estimate, and suggests that the association between ER IHC status and survival outcome is not statistically significant.  
The resulting coefficient (or log hazard ratio) = 0.0369, which is a positive value and indicates that having *Positive* ER IHC status is associated with an increased risk of survival outcome (i.e., death) compared to the *Negative* status; However, since the value is so close to zero, this suggests that there is no meaningful difference in risk.  
The calculated exponentiated coefficient (or hazard ratio) = 1.0376, implies that *Positive* ER IHC status have 3.76% increased risk of death compared with *Negative* status; However, since the hazard ratio is close to **1** the difference in risk is minimal.  
  
  
```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Fit Kaplan-Meier survival curve
km_fit <- survfit(formula = surv ~ ER_Status_By_IHC, data = clin_data3.Surv)
#summary(km_fit)

# Plot the Kaplan-Meier curve
plot(km_fit, col = c("blue", "red"), 
     xlab = "Overall Survival Months", 
     ylab = "Survival Probability",
     main = "Survival Curves for ER IHC Status",
     cex.lab = 1.5,
     cex.axis = 1.5,
     cex.main = 1.8,
     lwd = 2)
legend("topright", 
       legend = levels(as.factor(clin_data3$ER_Status_By_IHC)), 
       col = c("blue", "red"),
       cex = 1.5,
       lwd = 3)
```
  
  
  
# Question 3  
## Comparison of logistic regression (Q1) and Cox proportional hazard model (Q2)  
  
#### Which model demonstrates higher accuracy?  
Both the logistic regression model and the Cox proportional hazard model failed to demonstrate a statistically significant association between the input and output. So I suppose I could argue that both are just as accurate in this case. However, since the context of the analysis is concerned with looking at the amount of time that elapses before an event occurs, the Cox proportional hazard model would be more suitable for survival analysis compared to the logistic regression. Logistic regression ignores the time component and only analyzes a binary outcome. Additionally, the Cox model provides a *concordance index* that is a measure of predictive accuracy, and in this case the result of 0.539 indicates a weak ability (a value of 0.5 equals random guessing). Based on this, I now could argue that the Cox proportional hazard model demonstrates higher accuracy.  
  
#### Why is “Overall_Survival_Status” necessary when “overall survival months” already serves as an output?  
In survival analysis, Cox model, both variables are needed. The overall survival months captures when the event happens (i.e., time component), while the overall survival status captures whether the event occurred (i.e., alive or deceased). Without the overall survival status, the model won't be able to differentiate between subjects who experienced the event/outcome. Furthermore, it would be assumed that all the subjects experienced the event (i.e., death) which can lead to biased results and/or incorrect survival estimates.  
  
  
  
      
# Question 4  
## Identifying genes most significantly associated with "overall survival months" and "overall survival status" in tumor samples  
  
The code below will first subset the mRNA expression data set to include only the tumor sample IDs and the corresponding gene expressions. Then the median for each gene is calculated and stored in a new column appended in the data set. Finally, a new variable is created to store gene expressions greater-than **1**, which will be used for subsequent analysis. (This is done because having gene expression values as zero creates errors in subsequent codes.)  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Get the tumor data set
tumor_data <- mrna_data[, 1:104]
# rename the first column to gene
colnames(tumor_data)[1] <- "gene"

# Calculate the median for each gene and append as a new column
median_expr <- tumor_data %>%
  rowwise() %>%
  mutate(MedianExpression = median(c_across(starts_with("TCGA")))) %>%
  ungroup()

# Sort data & exclude median expression < 1
tumor_data2 <- median_expr %>%
  arrange(desc(MedianExpression)) %>% # sort in descending order
  filter(MedianExpression >= 1) # exclude median expression < 1
# show the generated data
tumor_data2[1:3, 103:105]
```
  
The code below will prepare the necessary data from both the clinical and mRNA expression data, and merge them into one data frame. Both *overall survival months* and *overall survival status* along with their corresponding sample IDs from the clinical data are needed together with the expression data.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
## Prepare mRNA data
# Transpose the expression data (exclude gene names)
mrna_tumor <- t(tumor_data2[, -c(1, 105)])
colnames(mrna_tumor) <- tumor_data2$gene # assign gene names to columns
mrna_tumor_df <- as.data.frame(mrna_tumor) # convert to data frame
# Add patient IDs & move as the first column
mrna_tumor_df$Sample_ID <- rownames(mrna_tumor_df)
mrna_tumor_df <- mrna_tumor_df %>%
  select(Sample_ID, everything()) # move column for last to first
# show the first row of the first 5 columns
mrna_tumor_df[1:5, 1:5]

## Prepare clinical data
# Subset clinical data for tumor samples
clin_data_tumor <- clin_data %>%
  select(Sample_ID = Sample_ID_Tumor, 
         survival_time = Overall_Survival_Months,
         survival_status = Overall_Survival_Status) %>%
  na.omit()
# Convert survival status to binary
clin_data_tumor$survival_status <- ifelse(clin_data_tumor$survival_status == "1:DECEASED", 1, 0)

# Merge mRNA and clinical data
merged_data <- merge(clin_data_tumor, mrna_tumor_df, by = "Sample_ID")
# show the first 5 rows & columns of the merged data
merged_data[1:5, 1:5]

```
  
The code below will perform the hold-out cross-validation method, using 60% of the data as training set and 40% as the test set. The median expression of each gene will be calculated again (within the code function) to classify the expression as high or low for input into the model. Then, the Cox proportional hazard model will be performed, the p-values will be recorded for each gene, and each p-value will be adjusted with multiple hypothesis testing.  
  
```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
# Split data into 60% training and 40% testing
train_idx <- sample(1:nrow(merged_data), 0.6 * nrow(merged_data))
train_data <- merged_data[train_idx, ]
test_data <- merged_data[-train_idx, ]

# Function to perform Cox regression for each gene
cox_analysis <- function(data) {
  genes <- colnames(data)[!(colnames(data) %in% c("Sample_ID", "survival_time", "survival_status"))]
  results <- data.frame(gene=character(), p_value=numeric(), hazard_ratio=numeric())
  
  for(gene in genes) {
    # calculate median expression for the gene in training data
    median_expr <- median(data[[gene]])
    
    # create high/low expression groups
    expr_group <- ifelse(data[[gene]] > median_expr, "High", "Low")
    
    # create Surv object
    surv_obj <- Surv(data$survival_time, data$survival_status)
    
    # perform cox regression
    cox_model <- tryCatch({
      coxph(surv_obj ~ expr_group)
    }, error = function(e) NULL)
    
    if(!is.null(cox_model)) {
      summary <- summary(cox_model)
      new_row <- data.frame(
        gene = gene,
        p_value = summary$waldtest["pvalue"],
        hazard_ratio = exp(coef(cox_model)),
        stringsAsFactors = FALSE
      )
      results <- dplyr::bind_rows(results, new_row)
    }
  }
  return(results)
}

# Execute function with training data
train_analysis <- cox_analysis(train_data)
# Adjust p-values for multiple testing
train_analysis$p_adj <- p.adjust(train_analysis$p_value, method = "BH")
# Sort by adjusted p-values and get top 5 genes
train_top5_genes <- train_analysis %>%
  arrange(p_adj) %>%
  head(5)
print(train_top5_genes)
```
  
The next code will validate the model on the test data set generated previously. The p-values of the training set's top 5 genes, resulting from both the training and testing set, will be reported. Finally, the top 5 genes obtained from the testing set will be reported as well.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
## Validate on test data
testSet_analysis <- cox_analysis(test_data)
# Adjust p-values for multiple testing
testSet_analysis$p_adj <- p.adjust(testSet_analysis$p_value, method = "BH")

# Compare results between training and test sets for top 5 genes
validation <- testSet_analysis %>%
  filter(gene %in% train_top5_genes$gene) %>%
  arrange(match(gene, train_top5_genes$gene))

print(validation)

# Get top 5 genes from test set
test_top5_genes <- testSet_analysis %>%
  arrange(p_adj) %>%
  head(5)
print(test_top5_genes)
```
  
The corresponding Kaplan-Meier plots for the top 5 genes resulting from the testing set are displayed below.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=10, fig.width=8, results='show', cache=FALSE}

## Plot Kaplan-Meier Curve for testing set's top 5 genes
# Function to plot KM curves
plot_km_curves <- function(test_data, tumor_data2, test_top5_genes) {
  # plot parameters
  par(mfrow = c(3, 2))
  # for loop to plot each gene
  for(gene in test_top5_genes$gene) {
    # make surv object
    surv_obj <- Surv(test_data$survival_time, test_data$survival_status)
    
    # get median expression
    median_expr <- tumor_data2$MedianExpression[tumor_data2$gene == gene]
    
    # classify into group based on median
    expr_group <- ifelse(test_data[[gene]] > median_expr, "High", "Low")
    
    # make KM curve using expression group
    km_curve <- survfit(surv_obj ~ expr_group)
    
    # plot the KM curve
    plot(km_curve, 
         main = paste("Kaplan-Meier Curve for", gene),
         xlab = "Time (Months)", 
         ylab = "Survival Probability",
         col = c("blue", "magenta"))
    
    # add legend
    legend("bottomleft", 
           legend = c("Low Expression", "High Expression"),
           col = c("blue", "magenta"),
           lty = 1)
  }
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
}

# Execute the function
plot_km_curves(test_data, tumor_data2, test_top5_genes)

```
  
    
      
      
# Question 5  
## External validation with outside dataset and comparison with top 5 genes generated previously  
  
The data set obtained from cbioportal.com is part of The Cancer Genome Atlas (TCGA) program which collected clinicopathologic annotation data with multi-platform molecular profiles of >11,000 human tumors of 33 different cancer types. The gene expression data contains 20,531 rows for genes and 1,084 columns for sample IDs.  

The code below will load the clinical and mRNA expression data sets into the R environment. The following codes will subset both data sets to clean-up and gather all relevant data for further analysis.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Load data sets into R environment
brca_clin_data <- read.csv("./brca_tcga_pan_can_atlas_clinical-data.csv")
brca_mrna_data <- read.csv("./brca_tcga_pan_can_atlas_mRNA.csv")
brca_mrna_data <- as.data.frame(brca_mrna_data)

## Prepare clinical data ##
# Subset relevant columns from clinical data
brca_clinical <- brca_clin_data %>%
  select(Sample_ID = Sample_ID, 
         survival_time = Overall_Survival_Months,
         survival_status = Overall_Survival_Status)
# Convert survival status to binary
brca_clinical$survival_status <- ifelse(brca_clinical$survival_status == "1:DECEASED", 1, 0)

head(brca_clinical)


## Prepare mRNA expression data ##
# Subset mrna data set: exclude missing values & second column
brca_seq_data <- brca_mrna_data %>%
    select(-c(2)) %>% # exclude second column
    filter(!is.na(Hugo_Symbol)) %>% # exclude rows with missing values
    rename(gene = Hugo_Symbol) # rename the first column to gene

# Calculate median for each gene and append as a new column
median_genes <- brca_seq_data %>%
  rowwise() %>%
  mutate(MedianExpression = median(c_across(starts_with("TCGA")))) %>%
  ungroup()

# Sort data & exclude median expression < 1
brca_seq_data2 <- median_genes %>%
  arrange(desc(MedianExpression)) %>% # sort in descending order
  filter(MedianExpression >= 1) # exclude median expression < 1

# Transpose the expression data (exclude gene names)
mrna_seq <- t(brca_seq_data2[, -c(1, 1084)])
colnames(mrna_seq) <- brca_seq_data2$gene # assign gene names to columns
mrna_seq_df <- as.data.frame(mrna_seq) # convert to data frame

# Add patient IDs & rearrange columns
mrna_seq_df$Sample_ID <- rownames(mrna_seq_df)
# Substitute "." with "-" in Sample_ID
mrna_seq_df$Sample_ID <- gsub("\\.", "-", mrna_seq_df$Sample_ID)
# Remove duplicate columns keeping first replicate
mrna_seq_df <- mrna_seq_df[, !duplicated(colnames(mrna_seq_df))] %>%
  select(Sample_ID, everything())
# show the first 5 row & columns
mrna_seq_df[1:5, 1:5]


## Merge mRNA and clinical data ##
data_merged <- merge(brca_clinical, mrna_seq_df, by = "Sample_ID")
# show the first 5 rows & columns of the merged data
data_merged[1:5, 1:5]

```
  
The code below will perform Cox proportional hazard model on the entire data, instead of splitting the data into training & testing sets. The median expression of each gene will be calculated again (within the code function) to classify the expression as high or low for input into the model. Then, the Cox proportional hazard model will be performed, the p-values will be recorded for each gene, and each p-value will be adjusted with multiple hypothesis testing.
  
```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
## Cox Regression ##
set.seed(123)
# Function to perform Cox regression for each gene
cox_analysis2 <- function(data) {
  genes <- colnames(data)[!(colnames(data) %in% c("Sample_ID", "survival_time", "survival_status"))]
  results <- data.frame(gene=character(), p_value=numeric(), hazard_ratio=numeric())
  
  for(gene in genes) {
    # calculate median expression for the gene in training data
    median_genes <- median(data[[gene]])
    
    # create high/low expression groups
    expr_group <- ifelse(data[[gene]] > median_genes, "High", "Low")
    
    # create Surv object
    surv_obj <- Surv(data$survival_time, data$survival_status)
    
    # perform cox regression
    cox_model <- tryCatch({
      coxph(surv_obj ~ expr_group)
    }, error = function(e) NULL)
    
    if(!is.null(cox_model)) {
      summary <- summary(cox_model)
      new_row <- data.frame(
        gene = gene,
        p_value = summary$waldtest["pvalue"],
        hazard_ratio = exp(coef(cox_model)),
        stringsAsFactors = FALSE
      )
      results <- dplyr::bind_rows(results, new_row)
    }
  }
  return(results)
}

# Execute the function
brca_analysis <- cox_analysis2(data_merged)
# Adjust p-values for multiple testing
brca_analysis$p_adj <- p.adjust(brca_analysis$p_value, method = "BH")

head(brca_analysis)
```
  
The next code will compare the results of the model on the external data set with the top 5 genes generated from test data set in Question #4. The results of the model from the previous testing set's top 5 genes, as well as the external data set, will be reported for comparison.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Compare results of top 5 genes from Q4 with external data
validation2 <- brca_analysis %>%
  filter(gene %in% test_top5_genes$gene) %>%
  arrange(match(gene, test_top5_genes$gene))

print(validation2)
print(test_top5_genes)
```
  
The corresponding Kaplan-Meier plots for the top 5 genes resulting from the external data set are displayed below.  

```{r , eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=10, fig.width=8, results='show', cache=FALSE}
## Plot Kaplan-Meier Curve for top 5 genes using external data results
# Function to plot KM curves
plot_km_curves2 <- function(data_merged, brca_seq_data2, test_top5_genes) {
  # plot parameters
  par(mfrow = c(3, 2))
  # for loop to plot each gene
  for(gene in test_top5_genes$gene) {
    # make surv object
    surv_obj <- Surv(data_merged$survival_time, data_merged$survival_status)
    
    # get median expression
    median_genes <- brca_seq_data2$MedianExpression[brca_seq_data2$gene == gene]
    
    # classify into group based on median
    expr_group <- ifelse(data_merged[[gene]] > median_genes, "High", "Low")
    
    # make KM curve using expression group
    km_curve <- survfit(surv_obj ~ expr_group)
    
    # plot the KM curve
    plot(km_curve, 
         main = paste("Kaplan-Meier Curve for", gene),
         xlab = "Time (Months)", 
         ylab = "Survival Probability",
         col = c("#00b7ff", "#ff5e00"))
    
    # add legend
    legend("bottomleft", 
           legend = c("Low Expression", "High Expression"),
           col = c("#00b7ff", "#ff5e00"),
           lty = 1)
  }
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
}

# Call the function with your data
plot_km_curves2(data_merged, brca_seq_data2, test_top5_genes)

```
  
**Discussion**  
There appears to be some consistent trends across both original and external data sets, where **COL1A1** and **EEF1A1** show high expression is associated with worse overall survival. While **GAPDH** shows consistent trend where high expression is associated with slightly better overall survival. On the other hand, **ACTG1** and **AHNAK** genes show inconsistent results between the data sets, indicating uncertainty in its association with survival.  
The observed inconsistencies suggest that the Cox model may have nicely performed for the original data set, as it displayed robust association results for some genes, however it did not generalize well to the external data set. This may reflect the differences between the two data sets, in terms of patient population or biological variability, for example.
