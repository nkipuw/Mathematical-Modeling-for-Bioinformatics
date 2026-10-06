---
title: "RBIF111_Homework2"
author: "Neshita Kipuw"
date: "Last update: `r format(Sys.time(), '%d %B, %Y')`"
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
```

The data set obtained is a clinical record containing the medical records of 299 patients who had heart failure, collected during their follow-up period, where each patient profile has 13 clinical features. The code below loads the data set and partitions the data into two outcome groups, "Survived" and "Death".

```{r viewDataset, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
##Data set upload & view
heart_data <- read.csv("heart_failure_clinical_records_dataset.csv", header=TRUE)
head(heart_data)

# Create data subset
subset_data <- function(data){
    #Get numeric columns that are not binary
    numeric_cols <- sapply(data, function(col) {is.numeric(col) && !all(col %in% c(0, 1))})
    #Subset the numeric columns and DEATH_EVENT column
    subset_cols <- c(which(numeric_cols), which(names(data) == "DEATH_EVENT"))
    heart_datav2 <- data[, subset_cols, drop = FALSE]
    return(heart_datav2)
}
# Execute the function
heart_datav2 <- subset_data(heart_data)
head(heart_datav2)

##Partition data set based on DEATH_EVENT
partition_data <- function(data){
    #Partition the data set into two based on DEATH_EVENT
    partition_data <- split(data, data$DEATH_EVENT)
    
    # Convert the list to a data frame
    df_0 <- partition_data[[1]]
    df_1 <- partition_data[[2]]
    
    # Add a column to identify the group
    df_0$group <- "Survived"
    df_1$group <- "Death"
    
    # Combine the two data frames
    result_df <- rbind(df_0, df_1)
    
    return(result_df)
}
# Execute the function
heart_datav3 <- partition_data(heart_datav2)
head(heart_datav3)
```

# Question 1

## Normality test for each data feature, including Shapiro-Wilks normality test and qq-plot.
The null hypothesis for Shapiro-Wilk test is that the data comes from a normal distribution. If the p-value is small, i.e., < 0.05, the null hypothesis is rejected and the data is concluded to be *not normally distributed*. The qq-plot compares the quantiles of the numerical data with the quantiles of a normal distribution. If the points of the numerical data superimposes on the line (will be *red* in this case), the data is estimated to be normally distributed.

```{r normalityTest, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Shapiro-Wilk Normality Test for numerical features of the dataset
set.seed(123)
sw_test <- function(data) {
  # Exclude the first column and 'DEATH_EVENT' column
  columns_to_test <- setdiff(names(data)[-1], c("DEATH_EVENT", "group"))
  
  for (col in columns_to_test) {
    if (is.numeric(data[[col]])) {
      cat("\nTesting normality for:", col, "\n")
      
      # Separate data for each group
      data_survived <- data[data$group == "Survived", col]
      data_death <- data[data$group == "Death", col]
      
      # Shapiro-Wilk test for Survived group
      sw_survived <- shapiro.test(data_survived)
      cat("Shapiro-Wilk test for Survived group:\n")
      print(sw_survived)
      
      # Shapiro-Wilk test for Death group
      sw_death <- shapiro.test(data_death)
      cat("Shapiro-Wilk test for Death group:\n")
      print(sw_death)
      
      # QQ plots
      par(mfrow = c(1,2))
      qqnorm(data_survived, main = paste("QQ Plot for", col, "(Survived)"))
      qqline(data_survived, col="red")
      qqnorm(data_death, main=paste("QQ Plot for", col, "(Death)"))
      qqline(data_death, col="red")
      
      # Reset plot layout
      par(mfrow = c(1,1))
    }
  }
}
# Execute the function
sw_test(heart_datav3)
```

It seems the data contain p-values that are mostly < 0.05, therefore the null hypothesis is rejected, and most of these features are considered *not normally distributed*. The qq-plot reveal that most of the data do not superimpose well with the red line, thereby not indicating a normal distribution which correlates with the results of the Shapiro-Wilk Test. Only a couple of qq-plots for platelets level and serum sodium level, show the closest superimposition of the data to the red line. However, they are still considered *not normal distribution* due to the low p-values.
Some of the columns in the data set are binary (contains value of only 0 or 1) and not suited for normality testing, therefore these columns were skipped.

# Question 2

## Box plot of each column in the data set

The values in some of the columns in the data set are quite large, and will not display well in a plot altogether due to the significant difference in scale. Therefore, the platelets with very large values will be converted to log10 scale, and the creatinine_pk and time columns will be converted to log2 scale. 
```{r logScale, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Convert some columns to log10 and log2 scale to better visualize the data
heart_datav4 <- subset(heart_datav3, select = -DEATH_EVENT)
heart_datav4$platelets <- log(heart_data$platelets)
heart_datav4$creatinine_phosphokinase <- log2(heart_data$creatinine_phosphokinase)
heart_datav4$time <- log2(heart_data$time)
#heart_datav4$serum_sodium <- log2(heart_data$serum_sodium)

# View the converted data
head(heart_datav4)
```

#### The code below will create a box plot for each numerical column (skipping the binary columns) in the data set according to the death event.
```{r boxplots, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Convert the data set to a data frame
heart_datav5 <- as.data.frame(heart_datav4)

# Create a long data frame for ggplot2
long_datav5 <- pivot_longer(heart_datav5, cols = -group, names_to = "variable", values_to = "value")

# Create boxplots for all columns
ggplot(long_datav5, aes(x = variable, y = value, fill = factor(group))) +
    geom_boxplot() +
    scale_fill_manual(values = c("Survived" = "yellow", "Death" = "magenta"), name = "group", labels = c("Death", "Survived")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 13)) +
    labs(title = "Box plots for All Columns", x = "Features", y = "Values")
```

# Question 3

## Hypothesis test comparing paired columns
To perform hypothesis test comparing paired columns of the data set, the values must be continuous. Since some of the features have binary values, these columns will not be used. The data subset generated in question 2 above will be used for the hypothesis testing.
Based on the normality test, the features in this data set is not considered normal distribution; therefore, the common t-test is not suitable and instead the Wilcoxon test (formula given below) will be peformed, which is the equivalent to t-test for non-parametric tests. 
$$w = \binom{1, y[j] < x[i]}{0, otherwise}\\
W = \displaystyle \sum_{i,j} w$$

#### The code below will perform the hypothesis testing

```{r wilcoxon, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Show the data
head(heart_datav4)

# Create function to perform Wilcoxon test for a single column
run_wilcox_test <- function(data, column) {
  survived_data <- data[data$group == "Survived", column]
  death_data <- data[data$group == "Death", column]
  
  test_result <- wilcox.test(survived_data, death_data, paired = FALSE, conf.int = TRUE)
  
  return(list(
    column = column,
    statistic = test_result$statistic,
    p_value = test_result$p.value,
    lower_conf_int = test_result$conf.int[1],
    upper_conf_int = test_result$conf.int[2]
  ))
}

# Function to run Wilcoxon tests on all numeric columns
wilcox_tests <- function(data) {
  # Get names of numeric columns (excluding 'group' and 'DEATH_EVENT')
  numeric_cols <- names(data)[sapply(data, is.numeric) & !names(data) %in% c("group")]
  
  # Perform Wilcoxon test for each numeric column
  results <- lapply(numeric_cols, function(col) {
    run_wilcox_test(data, col)
  })
  
  # Convert results to a data frame
  wilcox_df <- do.call(rbind, lapply(results, data.frame))
  return(wilcox_df)
}

# Execute the functions
wilcox_datav6 <- wilcox_tests(heart_datav4)

# Display results
print(wilcox_datav6)
```

Alternatively, the results can be sorted by p-value.
```{r sortpvalue, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
sorted_datav6 <- wilcox_datav6[order(wilcox_datav6$p_value), ]
print(sorted_datav6)
```

# Question 4

## Hypothesis test on lowest & greatest variance data

#### (Hidden) below are several functions to calculate the variance, find the highest and lowest variance, and perfom re-sampling and hypothesis tests. *Note: code hidden due to long length and may be expanded.*  

```{r functions, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-hide", fig.height=5, fig.width=8, results='show', cache=FALSE}
set.seed(123) ##for reproducibility

# Function to calculate variance and overlap
calc_var_overlap <- function(data, column) {
  survived <- data$group == "Survived"
  var_survived <- var(data[survived, column])
  var_death <- var(data[!survived, column])
  mean_survived <- mean(data[survived, column])
  mean_death <- mean(data[!survived, column])
  sd_survived <- sd(data[survived, column])
  sd_death <- sd(data[!survived, column])
  
  overlap <- min(mean_survived + sd_survived, mean_death + sd_death) - 
             max(mean_survived - sd_survived, mean_death - sd_death)
  
  return(list(variance = var_survived + var_death, overlap = overlap))
}

# Function to select column with highest variance and overlap
high_var <- function(data) {
  numeric_cols <- names(data)[sapply(data, is.numeric) & !names(data) %in% c("group", "DEATH_EVENT")]
  
  results <- sapply(numeric_cols, function(col) {
    res <- calc_var_overlap(data, col)
    res$variance * res$overlap
  })
  
  return(names(which.max(results)))
}

# Function to select column with lowest variance and overlap
low_var <- function(data) {
  numeric_cols <- names(data)[sapply(data, is.numeric) & !names(data) %in% c("group", "DEATH_EVENT")]
  
  results <- sapply(numeric_cols, function(col) {
    res <- calc_var_overlap(data, col)
    res$variance * res$overlap
  })
  
  return(names(which.min(results)))
}

# Function to perform re-sampling and hypothesis test
resample_and_test <- function(data, feature_col, outcome_col, sample_sizes, iterations) {
  results <- list()
  
  for (size in sample_sizes) {
    Pvals <- numeric(iterations)
    
    for (i in 1:iterations) {
      # Sample from both groups
      sample_index <- sample(nrow(data), size*2, replace = TRUE)
      sample_data <- data[sample_index, ]
      
      # Perform Wilcoxon test and store p-value
      test_result <- wilcox.test(formula(paste(feature_col, "~", outcome_col)), data = sample_data, exact = FALSE)
      Pvals[i] <- test_result$p.value
    }
    results[[as.character(size)]] <- Pvals
  }
  return(results)
}
```

#### The code below executes the functions above and determines the columns with the highest and lowest variances.

```{r variance, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Select the highest variance column
high_var_column <- high_var(heart_datav3)
print(paste("Highest variance column:", high_var_column))

# Select the lowest variance column
low_var_column <- low_var(heart_datav3)
print(paste("Lowest variance column:", low_var_column))
```

#### The code below executes the resampling and hypothesis testing between both groups and retains the p-values.

```{r resampling, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Set parameters
sample_sizes <- c(10, 15, 20)
iterations <- 1000

# Perform re-sampling and testing
results_high <- resample_and_test(heart_datav3, high_var_column, "DEATH_EVENT", sample_sizes, iterations)
results_low <- resample_and_test(heart_datav3, low_var_column, "DEATH_EVENT", sample_sizes, iterations)
#View(results)

```

#### The results for high variance are visualized with histogram using the code below.
```{r high_var, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Plot distribution of p-values
plot_data <- data.frame(
  p_value = unlist(results_high),
  sample_size = rep(names(results_high), each = iterations)
)

high_var_plot <- ggplot(plot_data, aes(x = p_value, fill = sample_size)) +
  geom_histogram(bins = 30, position = "dodge", alpha = 0.7) +
  facet_wrap(~ sample_size, ncol = 1) +
  labs(title = paste("P-value Histogram for", high_var_column),
       x = "P-value",
       y = "Count") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

# Get summary statistics
high_var_stats <- lapply(results_high, function(x) {
  c(mean = mean(x), median = median(x), sd = sd(x))
})

print(high_var_plot)
print(high_var_stats)
```

In question 3 platelets was calculated to have a p-value > 0.05, and here it was determined to have the highest variance. Due to the high p-value, the null hypothesis is not rejected and the two population outcome groups have the same distribution. After re-sampling 1000x and performing the test again, the p-values are still > 0.05, and the summary statistics show the median is consistent with increasing sample size. The results indicate that this is not due to chance and the two populations are confirmed to have the same distribution. 



#### The results for low variance are visualized with histogram using the code below.
```{r low_var, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, class.source = "fold-show", fig.height=5, fig.width=8, results='show', cache=FALSE}
# Plot distribution of p-values
plot_data <- data.frame(
  p_value = unlist(results_low),
  sample_size = rep(names(results_low), each = iterations)
)

low_var_plot <- ggplot(plot_data, aes(x = p_value, fill = sample_size)) +
  geom_histogram(bins = 30, position = "dodge", alpha = 0.7) +
  facet_wrap(~ sample_size, ncol = 1) +
  labs(title = paste("P-value Histogram for", low_var_column),
       x = "P-value",
       y = "Count") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

# Get summary statistics
low_var_stats <- lapply(results_low, function(x) {
  c(mean = mean(x), median = median(x), sd = sd(x))
})

print(low_var_plot)
print(low_var_stats)
```

In question 3 serum-creatinine was calculated to have a p-value much-smaller-than $\alpha$ = 0.05, and here it was determined to have the lowest variance. With a very small p-value. Therefore the null hypothesis is rejected, the two population outcome groups are deemed to not have the same distribution, and they are significantly different. After re-sampling 1000x and performing the test again, the p-values are observed to be < 0.05 and skewing to the left, indicating that the results is not due to chance and the two populations are confirmed to not have the same distribution. The summary statistics show that the re-sampling size matters in this case, as the median p-value show a decreasing trend with increasing sample size. If only the smallest sample size was re-sampled, the results would have changed the conclusion. 