---
title: "RBIF111 - HW5 - Midterm Exam"
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
  
  
The following (hidden) code will load the relevant files into the R environment for analysis.  
  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, cache=FALSE}
mRNA_data <- read.csv("TCGA_BRCA_mRNA_expression.csv")
protexpr_data <- read.csv("TCGA_BRCA_RPPA_protein_expression.csv")
clin_data <- read.csv("TCGA_BRCA_clinical_data.csv")
```
  
  

# Question 1  
## Perform un-/paired t-tests  
  
### The code below will perform the paired t-test on the mRNA expression data.  
The **Null Hypothesis** is that the population mean of the differences between the normal and tumor tissues equals zero. The alternative hypothesis is that the population mean of the differences between the normal and tumor tissues does not equal to zero.  
  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# log10 transform the data
mRNA_data_2 <- mRNA_data
mRNA_data_2[, 2:ncol(mRNA_data_2)] <- log10(mRNA_data_2[, 2:ncol(mRNA_data_2)] + 1)

# get normal and tumor ID from column names
types <- substr(colnames(mRNA_data_2), 14, 15)

# identify normal and tumor samples
normal_samples <- which(types == "11")
tumor_samples <- which(types == "01")

# Function to loop through all genes
ttest_results1 <- c()
mean_diff1 <- c()
for (gene in rownames(mRNA_data_2)) {
  normal_expr <- as.numeric(mRNA_data_2[gene, normal_samples])
  tumor_expr <- as.numeric(mRNA_data_2[gene, tumor_samples])

  # perform paired t-test
  paired_ttest <- t.test(tumor_expr, normal_expr, paired = TRUE)

  # store p-value
  ttest_results1[gene] <- paired_ttest$p.value

  # store mean difference
  mean_diff1[gene] <- paired_ttest$estimate
}

# new data frame of the generated p-values & mean differences
plot_df <- data.frame(mRNA_data_2[1], ttest_results1, mean_diff1)
head(plot_df)

```
  

```{r, eval=TRUE, echo=FALSE}
# Visualize with histograms
par(mfcol = c(1, 2))
# plots
hist(plot_df$ttest_results1,
  main = "Histogram of p-values (paired)",
  xlab = "log10 p-value",
  ylab = "Frequency")

hist(plot_df$mean_diff1,
  main = "Histogram of mean difference (paired)",
  xlab = "mean difference",
  ylab = "Frequency")

```
  
The result of the paired t-test show that a large distribution of p-values are close to zero in the p-value histogram. This suggests that most of the genes do not have statistically significant differences in expression between normal and breast cancer tissue. The expression data was transformed with log10, therefore p-values towards the right would indicate lower p-values and imply to be significant. Given the large p-value distribution near zero, there are likely only a few genes with very small p-values indicating significant difference between the groups. If there were more significant differences we would expect to see a higher distribution towards the right, at higher values of log10 p-value.  
The mean differences are centered around zero, showing normal distribution. This indicates that most genes have a small difference in expression between the two groups, where only a few genes are showing a large negative or positive mean difference. This distribution suggests that there is no widespread increase or decrease in expression across all genes when comparing normal to breast cancer tissue.  
  
  
### The code below will perform the unpaired t-test on the mRNA expression data.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Function to loop through all genes
ttest_results2 <- c()
mean_diff2 <- c()
for (gene in rownames(mRNA_data_2)) {
  normal_expr <- as.numeric(mRNA_data_2[gene, normal_samples])
  tumor_expr <- as.numeric(mRNA_data_2[gene, tumor_samples])

  # perform paired t-test
  paired_ttest2 <- t.test(tumor_expr, normal_expr, paired = FALSE)

  # store p-value
  ttest_results2[gene] <- paired_ttest2$p.value

  # store mean difference
  mean_diff2[gene] <- paired_ttest2$estimate
}

# new data frame of the generated p-values & mean differences
plot_df2 <- data.frame(mRNA_data_2[1], ttest_results2, mean_diff2)

## Histogram
par(mfcol = c(1, 2))
# plots
hist(plot_df2$ttest_results2,
  main = "Histogram of p-values - unpaired",
  xlab = "log10 p-value",
  ylab = "Frequency")
hist(plot_df2$mean_diff2,
  main = "Histogram of mean difference - unpaired",
  xlab = "mean difference",
  ylab = "Frequency")
```

```{r, eval=TRUE, echo=FALSE, class.source="fold-hide", results='hide'}
dev.off() ##reset graphics##
```  
  
**Discussion:**  
Similar to the unpaired results, the p-value histogram of the unpaired t-test show the majority of p-values are close to zero. This indicates that most genes don't show statistically significant differences in expression between the two groups.  
Unlike the paired result, the mean differences histogram of the unpaired t-test shows a distribution skewing to the right. This suggests that there may be an increase in expression in the breast cancer tissue compared to the normal tissue; However, as shown with the p-value histogram, these differences may not be statistically significant.
  
Based on the data, the most appropriate test is the **paired t-test** to determine whether the ERBB (HER2) mRNA expression is significantly different between normal and tumor tissues. Because the two sets of observations are uniquely paired so that an mRNA expression observation in one set matches the observation in the other. Most genes showed non-significant differences, but a few genes have small p-values. And many genes showed positive mean differences in expression, suggesting that cancer tissues might express certain genes at higher levels than normal tissues. From these results, the genes with very small p-values and large mean differences could be selected as important genes for further study.  
  
  

# Question 2  
## Correlation between mRNA expression and protein levels in HER2-IHC negative patients.  
  
The code below will first check the normality of the protein expression data as a whole.
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# check protein data for normality using Shapiro-Wilk test
prot_normtest <- as.data.frame(sapply(protexpr_data[, 2:107], shapiro.test))
# get just the statistic and p-value rows
prot_normtest <- prot_normtest %>% filter(row_number() <= n()-2)
# unlist the data values
prot_normtest <- do.call(data.frame, lapply(prot_normtest, function(x) unlist(x)))
# transpose the data into long format
prot_normlong <- t(prot_normtest)
# visualize the generated p-values
hist(prot_normlong[, "p.value"], breaks = 10)

```
  
The results of the Shapiro-Wilk test show that the majority of the p-values are very small, therefore indicating that the data is not normally distributed.  
  
  
  
The code below will subset the expression data from the mRNA and protein data sets by matching the patient ID with IHCHER2 *negative* status from the clinical data set. The subset data sets are then merged into one data frame to perform the correlation calculations.  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# get patient ID from clinical data
patient_ID <- clin_data %>% filter(IHCHER2 == "Negative") %>% select("Patient_ID")
# rename the patient ID to match the expression data
patient_ID$Patient_ID <- paste(substr(patient_ID$Patient_ID, 1, 13), "01", sep = "_")
# compile the patient ID to a vector
id_cols <- c(patient_ID$Patient_ID)

# get data from protein expression by patient ID
colnames(protexpr_data)[1] <- "gene" ##rename first column
prot_IHCHER2_neg <- subset(protexpr_data, select = c("gene", id_cols))
# arrange protein data by gene
prot_IHCHER2_neg <- prot_IHCHER2_neg %>% arrange(gene)


# get data from mRNA expression by patient ID
colnames(mRNA_data)[1] <- "gene" ##rename first column##
mrna_subset <- subset(mRNA_data, select = c("gene", id_cols))
## log10 transform the data to be in uniform scale with protein data
mrna_subset[, 2:ncol(mrna_subset)] <- log10(mrna_subset[, 2:ncol(mrna_subset)] + 1)


# merge data sets
merge_data <- inner_join(prot_IHCHER2_neg, mrna_subset, by = "gene")
merge_data <- merge_data %>% drop_na() ##drop any invalid values##

# Calculate pearson correlation & p-value
corr_function <- lapply(1:nrow(merge_data), function(x) {
    prot_values <- as.numeric(merge_data[x, 2:51])
    mrna_values <- as.numeric(merge_data[x, 52:101])

    pearson <- cor.test(mrna_values, prot_values, method = "pearson")
    spearman <- cor.test(mrna_values, prot_values, method = "spearman", exact = FALSE) ##exact = False helps eliminate ties##

    data.frame(
        gene = merge_data$gene[x],
        pearson_corr = pearson$estimate,
        pearson_pvalue = pearson$p.value,
        spearman_corr = spearman$estimate,
        spearman_pvalue = spearman$p.value
    )
})

# create data frame for results
corr_results <- do.call(rbind, corr_function)
head(corr_results)
```
  
  
The plots below visualizes the correlation results as histograms.
```{r, eval=TRUE, echo=FALSE, message=FALSE, warning=FALSE, class.source = "fold-hide", results='show', cache=FALSE}
# plot correlations
par(mfcol = c(2, 2))
hist(corr_results$pearson_corr,
  main = "Histogram of Pearson correlation",
  xlab = "Pearson correlation coefficient",
  ylab = "Frequency")
abline(v = 0.0, col = "blue", lty = "dashed")
hist(corr_results$pearson_pvalue, breaks = 20,
  main = "Histogram of Pearson p-value",
  xlab = "p-value",
  ylab = "Frequency")
abline(v = 0.05, col = "magenta", lty = "dashed")

hist(corr_results$spearman_corr,
  main = "Histogram of Spearman correlation",
  xlab = "Spearman correlation coefficient",
  ylab = "Frequency")
abline(v = 0.0, col = "blue", lty = "dashed")
hist(corr_results$spearman_pvalue, breaks = 20,
  main = "Histogram of Spearman p-value",
  xlab = "log10 p-value",
  ylab = "Frequency")
abline(v = 0.05, col = "magenta", lty = "dashed")

```
  
  
**Discussion:**
The Pearson and Spearman correlations measure the strength and direction of association between two continuous variables. The Pearson correlation measures the linear relationship between two variables, while the Spearman correlation measures the increasing or decreasing relationship between two variables which may not be necessarily linear.  
The Pearson correlation (r) assumes the data is normally distributed and the relationship between variables should be linear. The Null hypothesis here is that there is no linear relationship between the variables (r = 0). The alternative hypothesis is that there is a linear relationship (r $\neq$ 0). The histogram of the Pearson correlation shows close to normal distribution, but slightly skewing to the right. It also shows the majority of the correlation coefficient are larger-than zero (blue dash line), indicating a more positive linear relationship between mRNA and protein expression levels. The p-value histogram also shows most of the values are less-than $\alpha$ = 0.05 (red dash line), so we reject the null hypothesis as the data suggests a significant linear relationship.  
The Spearman correlation ($\rho$) makes no assumption of data normality and the relationship doesn't have to be strictly linear. The Null hypothesis is there is not an increasing/decreasing relationship between the variables ($\rho$ = 0), and the alternative hypothesis is that there is an increasing/decreasing relationship ($\rho$ $\neq$ 0). The histogram of the Spearman correlation show an almost normal distribution, with the majority of the correlation coefficient being larger-than zero (blue dash line). This indicates a more positive or increasing relationship between the mRNA and protein expression levels. The p-value histogram shows majority of the values are less-than $\alpha$ = 0.05 (red dash line), so we reject the null hypothesis as it suggests a significant linear relationship.  
Both correlation tests show statistical significance of a positive relationship between mRNA and protein expression levels. However, both the mRNA and protein expression data sets are not normally distributed, which does not meet the assumptions of the Pearson correlation; therefore, the most suitable method for this analysis would be the *Spearman correlation*.  
  
  
```{r, eval=TRUE, echo=FALSE, class.source="fold-hide", results='hide'}
dev.off() ##reset graphics##
```  
  
  
  
# Question 3  
## Relationship between HER2 & ER status in breast cancer patients based on IHC positive and negative status.  
  
  
To determine if there is a significant association between HER2 and ER status, the chi-square test will be used on the two categorical variables (HER2 & ER). The Null hypothesis is that there is no association between the two variables, meaning they are independent of each other. The alternative hypothesis is that there is an association between HER2 and ER status.  
  
First we obtain the relevant patient data that contains positive and negative status for HER2 and ER in the clinical data.  
  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# get patient ID & classification from clinical data
my_df <- clin_data %>% 
    filter(IHCHER2 %in% c("Positive", "Negative")) %>%
    filter(ER_Status_By_IHC %in% c("Negative", "Positive")) %>%
    select("Patient_ID", "IHCHER2", "ER_Status_By_IHC")
head(my_df)
```
  
Next, we create the contingency table from the data frame, and perform the chi-square test.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
# Create a contingency table for HER2 and ER status
status_table <- table(my_df$IHCHER2, my_df$ER_Status_By_IHC)
print(status_table)

# Perform chi-square test
chisq_test <- chisq.test(status_table)
print(chisq_test)

```
  
To determine whether the difference is statistically significant, we use Fisher's exact test. This will let us know that the difference is due to chance or that the two variables are independent.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
# Check if difference is statistically significant using Fisher's exact test
fisher_test <- fisher.test(status_table)
print(fisher_test)

```
  
  
**Discussion**  
Both the p-values from the chi-square test and the Fisher exact test are not small, therefore we cannot reject the Null hypothesis. It was found that there is no statistically significant association between HER2 and ER status, and it is likely that these two variables are independent.  
  
The code below help to visualize the results of the chi-square test. The first one is a mosaic plot, which displays the frequency of each combination of HER2 and ER status. This plot allows a quick visual check for independence, in which case the segments will align proportionally, and if there is an association the segments will not align proportionally. The plot shows that the segments are very close to be aligned proportionally, indicating independence of the two variables.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Create a mosaic plot to show the relationship between HER and ER
mosaicplot(status_table, main = "Mosaic Plot of HER2 and ER Status",
           xlab = "HER2 Status", ylab = "ER Status",
           color = TRUE)
```
  
  
# Question 4  
  
## Regression model with "fraction genome altered" as input and "tumor mutation burden" as the output, showing 95% CI for the regression coefficients.  
First we obtain the relevant data containing the values for "fraction genome altered" and "tumor mutation burden" from the clinical data. Then we perform the simple linear regression model, and calculate the 95% confidence interval, in the code below.  
  

```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123)
# obtain relevant data from clinical data
my_df2 <- subset(clin_data, select = c("Patient_ID", "Fraction_Genome_Altered", "TMB_nonsynonymous"))
head(my_df2)

# fit the linear regression model
my_lm <- lm(TMB_nonsynonymous ~ Fraction_Genome_Altered, data = my_df2)
summary(my_lm)

# get the 95% CI
confint(my_lm, level = 0.95)
```
  
**Discussion**  
Within the summary of the linear model, the intercept represents the expected value of the *tumor mutation burden* when the *fraction genome altered* is zero. The coefficient for *fraction genome altered* represents the average change in *tumor mutation burden* for each one-unit increase in genome alteration -- it indicates, that *tumor mutation burden* increases as genome alteration increases, and the opposite for a negative coefficient. The p-value associated with *fraction genome altered* coefficient is less-than $\alpha$ = 0.05. The Null hypothesis is that the coefficient is zero indicating there is no relationship. Therefore, a low p-value result suggests that there is a statistically significant relationship between the *fraction genome altered* (predictor variable) and *tumor mutation burden* (response variable).  
The result of 95% confidence interval indicates a more precise estimate of the coefficient when it is narrow, while a wide interval suggests greater uncertainty. For this data it seems that the resulting interval is pretty good as it is quite narrow.  
The R-squared value represents the proportion of variance in the response variable (tumor mutation burden) explained by the predictor variable (fraction genome altered). The value ranges from zero to one, where an R-squared value close to 1 indicates that a large portion of the variability in *tumor mutation burden* is explained by *fraction genome altered*. On the other hand, an R-squared value close to zero suggests that genome alteration does not explain much of the variability in tumor mutation. This value helps to understand the strength of the linear model. Higher values indicate a better fit, however it does not guarantee a model's predictive power.  
To summarize, the coefficient for *fraction genome altered* is 1.84, with a 95% confidence interval of 0.36-3.32, and a p-value of 0.0154. This suggests a positive relationship where *tumor mutation burden* increases by 1.84 units for each unit of genome alteration. With a low p-value, the relationship is statistically significant. However, the R-squared value is very small and does not indicate a good fit; therefore, it is difficult to determine whether the variability in the tumor mutation is explained by genome alteration. The graph below shows the scatter plot to visualize the data.  
  

```{r, eval=TRUE, echo=FALSE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Generate scatter plot with the regression line
ggplot(my_df2, aes(x = Fraction_Genome_Altered, y = TMB_nonsynonymous)) +
  geom_point(color = "black", size = 3) +                  # Plot points
  geom_smooth(method = "lm", color = "magenta", se = TRUE) +  # Add regression line with confidence interval
  labs(
    title = "Relationship between Genome Alteration and Tumor Mutation",
    x = "Fraction Genome Altered",
    y = "Tumor Mutation Burden"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

```
  
  
  
  
  
# Question 5  
  
## Identify genes most positively correlated with ERBB2 in normal and cancer samples.  
  
The code below will first obtain and divide the mRNA data into cancer and normal sample data. Then the Pearson correlation coefficient will be calculated, along with the p-value associated with the correlation. Further on, multiple hypothesis testing is utilized to control for false positives.  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# new variable for mRNA data
my_df3 <- mRNA_data

# separate cancer and normal data, keeping the gene column
cancer_data <- my_df3[, 1:104]
normal_data <- my_df3[, c(1, 105:207)]


### Calculate Pearson correlation & p-value

# function to calculate Pearson correlation & p-value
calc_corr <- function(df, reference_gene = "ERBB2") {
  # subset reference gene row for comparison
  ref_expression <- df %>% filter(gene == reference_gene) %>% select(-gene) %>% as.numeric()
  
  # make lists to store results
  cor_results <- list()
  
  # calculate Pearson correlation for each gene with the reference gene
  cor_results <- lapply(1:nrow(df), function(k) {
    gene_expression <- as.numeric(df[k, -1]) # Exclude gene name column
    pearson <- cor.test(gene_expression, ref_expression)
    
    # store gene, correlation, and p-value
    data.frame(
      gene = df$gene[k],
      Pearson_corr = pearson$estimate,
      p_value = pearson$p.value
    )
  })
  # make data frame for results
  my_corr_res <- do.call(rbind, cor_results)
  return(my_corr_res)
}

# Execute function for normal and cancer data
cancer_corr <- calc_corr(cancer_data)
normal_corr <- calc_corr(normal_data)

head(cancer_corr)
head(normal_corr)

```
  
  
Due to the large number of rows being analyzed, we increase the risk of obtaining false positives and raising the chance that some genes will appear statistically significant by chance. The *Benjamini-Hochberg* method will be used to control the false discovery rate (FDR) in this case. The method works by ranking all p-values in ascending order, then adjusting each one based on its rank and the total number of tests. This adjusts the p-value so that the probability of FDR is still within a threshold, often at 0.05.
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Multiple hypothesis testing correction
cancer_corr$p_adjust <- p.adjust(cancer_corr$p_value, method = "BH")
normal_corr$p_adjust <- p.adjust(normal_corr$p_value, method = "BH")

head(cancer_corr)
head(normal_corr)
```
  
  
In a Pearson correlation analysis, a **positive** correlation between a gene and the reference gene (ERBB2) show that the mRNA expression levels of the two genes tend to increase or decrease with each other. When one gene's expression level is high, the positively correlated gene is also likely to be high (and vice versa). Below are the top 10 positively correlated genes in reference to ERBB2, within the cancer and normal samples.  
  
  
```{r, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Identify top 10 correlated genes for tumor and normal
top10_cancer <- cancer_corr %>% filter(gene != "ERBB2") %>%
  arrange(desc(Pearson_corr)) %>% head(10)

top10_normal <- normal_corr %>% filter(gene != "ERBB2") %>%
  arrange(desc(Pearson_corr)) %>% head(10)

print(top10_cancer)
print(top10_normal)
```
  
  
For visualization of the results, below are the bar plots of the top 10 genes identified in the cancer and normal samples.
```{r, echo=FALSE}
# Plot top 10 correlated genes for tumor and normal samples
ggplot(top10_cancer, aes(x = gene, y = Pearson_corr, fill = gene)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Pearson_corr, 2)), vjust = -0.5, size = 4) +
  labs(title = "Top 10 Positively Correlated Genes with ERBB2",
       subtitle = "in Cancer Samples",
       x = "Gene", y = "Pearson Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggplot(top10_normal, aes(x = gene, y = Pearson_corr, fill = gene)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Pearson_corr, 2)), vjust = -0.5, size = 4) +
  labs(title = "Top 10 Positively Correlated Genes with ERBB2",
       subtitle = "in Normal Samples",
       x = "Gene", y = "Pearson Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

```