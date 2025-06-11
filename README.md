# GDSC_ML

## Introduction

This project is focused on building a machine learning model to predict IC50 values, using the GDSC2 (Genomics of Drug Sensitivity in Cancer) dataset.
This dataset contains drug response data across 969 different cell lines, tested with up to 286 compounds.
To arrive at this model, I will be performing exploration and manipulation of the dataset, statistical analysis and trying different predictive models.

## Dataset exploration

I began exploring the structure of the dataset and each individual variable. Generating a summary of the continuous data glimpses on the fact that the main variables used to evaluate drug response (LN_IC50 and AUC) are non-normally distributed, and the values spread over seven standard deviations from the mean. Also, the original data (IC50) is already normalized.

```{r message=FALSE, warning=FALSE}
library(tidyverse)

gdsc <- read_csv("C:\\Users\\Administrator\\Desktop\\GDSC\\GDSC_DATASET.csv", show_col_types = FALSE)
```

```{r}
head(gdsc)
glimpse(gdsc)
sapply(gdsc, function(x) length(unique(x)))
summary(gdsc[, c("LN_IC50", "AUC", "Z_SCORE")])
```

In the following histograms, it is clearly seen that both the LN_IC50 and AUC are not normally distributed and skewed to the left - this distribution is marked in the AUC histogram.

```{r}
ggplot(gdsc, aes(x=LN_IC50)) +
  geom_histogram(binwidth = 0.5, color='black', fill="#69b3a2", alpha=0.7) +
  theme_minimal() +
  labs(title="Histogram of LN_IC50", x="LN_IC50", y="Frequency")

ggplot(gdsc, aes(x=AUC)) +
  geom_histogram(binwidth = 0.02, color='black', fill="#404080", alpha=0.7) +
  theme_minimal() +
  labs(title="Histogram of AUC", x="AUC", y="Frequency")
```

Applying a Shapiro test to LN_IC50 corroborates a non-normal distribution, and applying skewness and kurtosis show a left skew.

```{r message=FALSE, warning=FALSE}
library(moments)
skewness(gdsc$LN_IC50)
kurtosis(gdsc$LN_IC50)
```

```{r}
shapiro.test(sample(gdsc$LN_IC50, 5000))
```

In the Q-Q plot bellow, we can see how data points deviate sharply from the diagonal line, especially at the tails, further suggesting non-linear distributions.

```{r}
qqnorm(gdsc$LN_IC50)
qqline(gdsc$LN_IC50)
```

The following plots explore the frequency of pathways targeted by all compounds, and the frequency of cancer by tissue types.

```{r}
cat("Exploration of the frequency of pathways targeted by all compounds\n")
ggplot(gdsc, aes(x=TARGET_PATHWAY)) +
  geom_bar(fill = "darkgreen", color = "#1A4C6F", alpha = 0.6) +
  labs(title = "Target Pathway Frequency", x = "Categories", y = "Count") +
  coord_flip() +
  theme_minimal()
```

```{r}
cat("Exploration of the frequency of cancers by tissue\n")
gdsc %>% group_by(`GDSC Tissue descriptor 2`) %>% summarise(count = n()) %>%
  ggplot(aes(x=reorder(`GDSC Tissue descriptor 2`,+count), y = count)) +
  geom_bar(stat = "identity", fill = "#8f4e46", color = "black", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Tissue type Frequency", x = "Tissue types", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
```

Summarizing LN_IC50 values by tissue, hematologic cancers show lower LN_IC50 means, while solid tumors - mainly gastrointestinal and lung cancers - and sarcomas show higher means.

```{r}
cat("Exploring a summary of LN_IC50 by tissue\n")
summary_by_tissue <- gdsc_final %>%
  group_by(tissue_descriptor_2) %>%
  summarise(mean = mean(LN_IC50), 
            median = median(LN_IC50), 
            SD = sd(LN_IC50), 
            n = n(), 
            max = max(LN_IC50), 
            min = min(LN_IC50)) %>%
  arrange(median)
summary_by_tissue
```

### Exploration of missing values.

```{r message=FALSE, warning=FALSE}
library(naniar)
```

The missing values of genomic markers and compound targets tend to groups together. Also, many missing or unclassified values on the variables "Cancer Type" and "TCGA_DESC", are better described by the variable "GDSC Tissue descriptor 2".

```{r}
miss_var_summary(gdsc)
gg_miss_upset(gdsc)
```

```{r}
cat("The following bar plot shows the mayority of unclassified values in TCGA_DESC, described by the variable GDSC Tissue descriptor 2\n")

par(mfrow = c(1, 2))

gdsc %>% group_by(TCGA_DESC) %>% summarise(count = n()) %>%
  ggplot(aes(
  x= reorder(TCGA_DESC,+count),
  y = count,
  fill = ifelse(reorder(TCGA_DESC,+count) == "UNCLASSIFIED", "#cc0000", "#f1c232"))
  ) +
  geom_bar(stat = "identity", color = "black", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Frequency of TCGA cancer type", x = "TCGA types", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  guides(fill = FALSE)

gdsc %>%
  filter(TCGA_DESC == "UNCLASSIFIED") %>%
    group_by(`GDSC Tissue descriptor 2`) %>%
      summarise(count = n()) %>%
        ggplot(aes(x=reorder(`GDSC Tissue descriptor 2`,+count), y = count)) +
        geom_bar(stat = "identity", fill = "#f1c232", color = "black", alpha = 0.6) +
        theme_minimal() +
        labs(title = "Frequency of TCGA cancer type", x = "TCGA types", y = "Count") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
```

Most cell lines are tested by approximately all 286 drugs, and some a tested twice for 8 specific compounds, as shown in the tables bellow.

```{r}
cat("The following compounds are present twice in many cell lines\n")
gdsc %>% 
  count(DRUG_NAME, CELL_LINE_NAME, sort = TRUE) %>% 
  filter(n > 1, n < 3) %>% 
  distinct(DRUG_NAME)
```

```{r}
cat("Number of tests run by cell line")
gdsc %>% count(CELL_LINE_NAME, CELL_LINE_NAME, sort = TRUE) %>%
  ggplot(aes(x = n)) +
  geom_histogram(binwidth = 20, color="black", alpha=0.6, fill="darkblue") +
  labs(title = "Histogram of tests per cell line", x = "Number of tests", 
       y = "Number of cell lines") +
  theme_minimal()
```

With this insights, I proceeded to arrange and transform the original dataset, clearing out missing values, to facilitate further analysis. 

```{r}
gdsc_final <- gdsc %>%
  select(-COSMIC_ID, -`Cancer Type (matching TCGA label)`, -`Screen Medium`, -`Growth Properties`, -DRUG_ID, -`GDSC Tissue descriptor 1`, -TCGA_DESC) %>%
  rename(msi = `Microsatellite instability Status (MSI)`,
         cna = CNA,
         methylation = Methylation,
         gene_expression = `Gene Expression`,
         tissue_descriptor_2 = `GDSC Tissue descriptor 2`) %>%
  mutate(msi = recode(msi, "MSI-H"=1, "MSS/MSI-L"=0),
         gene_expression = recode(gene_expression, "Y"=1, "N"=0),
         cna = recode(cna, "Y"=1, "N"=0),
         Methylation = recode(methylation, "Y"=1, "N"=0)) %>%
  drop_na(c(msi, gene_expression, cna, methylation))
```

## Exploration of drug sensitivity by cell line and compound

The following exploration allows to dive deeper and dissect the data by cell line, compound and tissue.

This function plots all IC50 Z-score values for a given cell line. Bellow, an example for a cell line of breast cancer (TCGA classification: BRCA)

```{r}
Graph_IC50_ZSCORE_by_Cell_Line <- function(dataframe, cell_line){
  dataframe %>% filter(CELL_LINE_NAME == cell_line) %>%
  ggplot(aes(x = reorder(DRUG_NAME, Z_SCORE), y = Z_SCORE)) +
  geom_point() +
  labs(title = sprintf("Drug sensitivity by cell line %s", cell_line), x = cell_line, y = "Cell line IC50 Z-Score") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
}
```

```{r}
Graph_IC50_ZSCORE_by_Cell_Line(gdsc_final, "MFM-223")
```

This next function allows for plotting all LN_IC50 values for a given compound, useful for exploration of individual compounds. Bellow, an example using the drug Dactinomycin is shown, colored by tissue type.

```{r}
Graph_IC50_by_Drug <- function(dataframe, drug_name){
  dataframe %>% filter(DRUG_NAME == drug_name) %>%
  ggplot(aes(x = reorder(CELL_LINE_NAME, LN_IC50), y = LN_IC50)) +
  geom_point(aes(color = tissue_descriptor_2), size = 1.5) +
  theme_classic() +
  labs(title = sprintf("IC50 by compound %s", drug_name), x = " ", y = "LN_IC50") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
}
```

```{r}
Graph_IC50_by_Drug(gdsc_final, "Dactinomycin")
```

The following graph shows the distribution of LN_IC50 across tissue types. Bellow, a function allowing exploration of the LN_IC50 distribution for a provided compound across tissue types, with an example using Crizotinib.

```{r}
gdsc_final %>%
  ggplot(aes(x = tissue_descriptor_2, y = LN_IC50, fill = tissue_descriptor_2)) +
  geom_boxplot(outliers = F) +
  theme_minimal() +
  labs(title = "LN_IC50 distribution by tissue type", x = "Tissue type", y = "LN IC50") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")
```

```{r}
Graph_LNIC50_by_TissueType <- function(dataframe, drug_name){
  dataframe %>% filter(DRUG_NAME == drug_name) %>%
  ggplot(aes(x = tissue_descriptor_2, y = LN_IC50, fill = tissue_descriptor_2)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = sprintf("Sensitivity to %s by tissue type", drug_name), x = "Tissue types", y = "LN_IC50") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")
}
```

```{r}
Graph_LNIC50_by_TissueType(gdsc_final, "Crizotinib")
```

Lastly, this function allows for the exploration of the distribution of LN_IC50 across different target pathways, for a given cell line. The distribution for the breast cancer cell line MFM-223 is shown as an example.

```{r}
Graph_pathway_by_CellLine <- function(dataframe, cell_line){
  dataframe %>% filter(CELL_LINE_NAME == cell_line) %>%
  ggplot(aes(x = TARGET_PATHWAY, y = LN_IC50, fill = TARGET_PATHWAY)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = sprintf("Pathway sensitivity by cell line %s", cell_line), x = "Target Pathways", y = "LN_IC50") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")
}
Graph_pathway_by_CellLine(gdsc_final, "MFM-223")
```

## Correlation analysis

Here we explore the correlation between two continuous variables in our dataset, AUC and LN_IC50

We know that a lower **LN_IC50** equals a higher drug sensitivity (cells respond to lower drug doses), and that a lower **AUC** equals a stronger overall drug effect (less cell viability). This suggests that as **drug resistance increases** (higher LN_IC50), cell viability **also increases** (higher AUC). In other words, we’d expect a *strong positive correlation* between these metrics.

We used a Spearman correlation test and a scatter plot to visualize this correlation.

```{r}
round(cor(gdsc_final$AUC, gdsc_final$LN_IC50, method = "spearman"), 3)
ggplot(gdsc_final, aes(x = LN_IC50, y = AUC)) +
  geom_point(shape = 20, size = 1) +
  labs(title = "LN_IC50 vs. AUC",
       x = "LN_IC50",
       y = "AUC") +
  theme_classic()
```

The scatter plot shows the trend clearly: At low LN_IC50, AUC values are also low. As LN_IC50 increases, AUC rises, but the curve tends to flatten at higher values, suggesting diminishing returns in resistance.

The correlation reinforces the consistency of GDSC’s drug-response measurements. It also underscores that resistance mechanisms often overlap with increased cell survival pathways.

For further analysis, I developed a function for exploring the correlation between these variables in different types of cancer. Providing a tissue from the variable \`GDSC tissue descriptor 2\`, we obtain the spearman correlation and a scatter plot. Bellow are shown two examples, denoting a stronger correlation for chronic myeloid leukaemia, compared to breast cancer.

```{r}
spearmanCorr_byTissue <- function(tissue){
  corr_df <- gdsc_final |>
    filter(tissue_descriptor_2 == tissue) |>
    select(LN_IC50, AUC)
  spearman_corr <- cor(corr_df$LN_IC50, corr_df$AUC, method = "spearman")
  cat(paste("Spearman Correlation between AUC and LN_IC50 for", tissue, ":", round(spearman_corr, 3), "\n"))
  
  ggplot(corr_df, aes(x = LN_IC50, y = AUC)) +
    geom_point(shape = 20, size = 1) +
    labs(title = "LN_IC50 vs. AUC",
         x = "LN_IC50",
         y = "AUC") +
    theme_classic()
}
spearmanCorr_byTissue("breast")
spearmanCorr_byTissue("chronic_myeloid_leukaemia")
```

## Statistical analysis

### Binary variables

The variables **CNA**, **methylation**, **gene_expression,** and **microsatellite instability (MSI)** encode yes/no or high/low genomic features. To study their influence in drug response, I performed a Wilcoxon signed-ranked test on each.

```{r}
wilcoxTest_MSI <- wilcox.test(gdsc_final$LN_IC50[gdsc_final$msi == 1],
                              gdsc_final$LN_IC50[gdsc_final$msi == 0],
                              var.equal = T)
wilcoxTest_MSI
```

```{r}
wilcoxTest_CNA <- wilcox.test(gdsc_final$LN_IC50[gdsc_final$cna == 1],
                              gdsc_final$LN_IC50[gdsc_final$cna == 0], 
                              var.equal = T)
wilcoxTest_CNA
```

```{r}
wilcoxTest_geneExpression <- wilcox.test(gdsc_final$LN_IC50[gdsc_final$gene_expression == 1],
                                         gdsc_final$LN_IC50[gdsc_final$gene_expression == 0],
                                         var.equal = T)
wilcoxTest_geneExpression
```

```{r}
wilcoxTest_methylation <- wilcox.test(gdsc_final$LN_IC50[gdsc_final$Methylation == 1],
                                      gdsc_final$LN_IC50[gdsc_final$Methylation == 0],
                                      var.equal = T)
wilcoxTest_methylation
```

Even though most of the tests are statistically significant, this data is asymmetric, lacks in depth and may oversimplify the biological complexity of their influence in cancer response. Therefore, they are not included in model training.

### Categorical variables

The categorical variables in the dataset can be divided into different types.

1.  Name variables for compounds and cell lines.
2.  Variables that describe cancer type: GDSC tissue descriptors, TCGA_DESC, Cancer Type
3.  Variables that describe the compound's presumed pathway and target of action: TARGET_PATHWAY, TARGET.
4.  A variable describing the cell line's growth properties.

I studied the significance of the second and third group of variables applying a Kruskal-Wallis test.

The variables TCGA_DESC and Cancer Type presented inconsistencies between them, with many missing or unclassified values. Therefore, I used GDSC Tissue descriptor 2, which offered the least amount of missing/unclassified values, and a more specific description of the cancer type.

```{r message=FALSE, warning=FALSE}
library(FSA, include.only = "dunnTest")
```

```{r}
gdsc_final$tissue_descriptor_2 <- as.factor(gdsc_final$tissue_descriptor_2)
tissueType_kt <- kruskal.test(LN_IC50 ~ tissue_descriptor_2, data = gdsc_final)
tissueType_kt
tissueType_dunnt <- dunnTest(LN_IC50 ~ tissue_descriptor_2, data = gdsc_final,
                             method = "bonferroni")
tissueType_dunnt
```

Then, I explored both variables describing the presumed pathway and putative target of each compound.

TARGET_PATHWAY allows for exploration of the general pathways by which many drugs work; whereas TARGET is generally more specific to a single drug, and, in many cases, unknown.

```{r}
gdsc_final$TARGET_PATHWAY <- as.factor(gdsc_final$TARGET_PATHWAY)
targetPathway_kt <- kruskal.test(LN_IC50 ~ TARGET_PATHWAY, data = gdsc_final)
targetPathway_kt
targetPathway_dunnt <- dunnTest(LN_IC50 ~ TARGET_PATHWAY, data=gdsc_final,
                                method = "bonferroni")
targetPathway_dunnt

gdsc_final$TARGET <- as.factor(gdsc_final$TARGET)
target_kt <- kruskal.test(LN_IC50 ~ TARGET, data = gdsc_final)
target_kt
target_dunnt <- dunnTest(LN_IC50 ~ TARGET, data = df)
target_dunnt
```

The following function allows for the exploration of relevant target pathways for a specific tissue. It returns a Dunn test providing the data-frame and a value from Tissue descriptor 2.

```{r}
targetPathway_by_tissue_dunnt <- function(dataframe, tissue) {
  dataframe$TARGET_PATHWAY <- as.factor(dataframe$TARGET_PATHWAY)
  tissue_df <- dataframe |>
    filter(tissue_descriptor_2 == tissue)
  dunn_test <- dunnTest(LN_IC50 ~ TARGET_PATHWAY, data = tissue_df)
  dunn_test
}
```
## Linear regression model

The previous analysis contributed to understanding of each individual variable and selecting relevant variables for modeling.
I tried different approaches for linear models, shown bellow.

```{r}
#| eval: false
#| include: false
lm(LN_IC50 ~ TARGET, data = gdsc_final)
lm(LN_IC50 ~ TARGET_PATHWAY, data = gdsc_final)
lm(LN_IC50 ~ tissue_descriptor_2, data = gdsc_final)
```

In this model, I evaluate how the cell line and drug alone describe drug response.

```{r}
mlr_model <- lm(LN_IC50 ~ CELL_LINE_NAME + DRUG_NAME,
                data = gdsc_final)
summary(mlr_model)
```

Here’s what I learned from some of my initial linear modeling attempts:

1. Testing Individual Predictors: Modeling LN_IC50 against single variables showed consistently low adjusted R² values, except for the molecular target, which alone explained 67% of the variance. This suggests that knowing a drug’s specific mechanism of action is remarkably predictive of its effectiveness across cell lines.

2. Testing cell lines and drug names howed the following results:
Residual standard error: 1.177 on 227,420 DF 
Adjusted R²: 0.818 
F-statistic: 864.6 (p < 2.2e-16) 

Summarizing these results, I draw the following conclusions:
- High R² suggests the model captures some structure
- Most explanatory power comes from cell line and drug names—essentially memorizing patterns rather than uncovering generalizable biology.

##### Evaluating the model:
In the Residuals vs. Fitted values plot bellow, we can assume heteroscedasticity by the increased dispersion of residuals by the left and middle section of the plot. The Q-Q plot of the standardized residuals show strong deviation from both tails, and the Scale-Location plot shows a non-straight line, further contributing to the assumption of non-normality.

```{r}
plot(lm_model)
```
