library(tidyverse)
library(naniar)

gdsc <- read_csv("C:\\Users\\Administrator\\Desktop\\GDSC\\GDSC_DATASET.csv", 
                 show_col_types = F)

# Dataset exploration
head(gdsc)
glimpse(gdsc)
sapply(gdsc, function(x) length(unique(x)))
summary(gdsc[, c("LN_IC50", "AUC", "Z_SCORE")])

# Distribution of continous variables
ggplot(gdsc, aes(x=LN_IC50)) +
  geom_histogram(binwidth = 0.5, color='black', fill="#69b3a2", alpha=0.7) +
  theme_minimal() +
  labs(title="Histogram of LN_IC50", x="LN_IC50", y="Frequency")

ggplot(gdsc, aes(x=AUC)) +
  geom_histogram(binwidth = 0.02, color='black', fill="#404080", alpha=0.7) +
  theme_minimal() +
  labs(title="Histogram of AUC", x="AUC", y="Frequency")

# Distribution of LN_IC50 by tissue and target pathway
ggplot(gdsc, aes(x=TARGET_PATHWAY)) +
  geom_bar(fill = "darkgreen", color = "#1A4C6F", alpha = 0.6) +
  labs(title = "Target Pathway Frequency", x = "Categories", y = "Count") +
  coord_flip() +
  theme_minimal()

gdsc_final %>% group_by(tissue_descriptor_2) %>% summarise(count = n()) %>%
  ggplot(aes(x=reorder(tissue_descriptor_2,+count), y = count)) +
  geom_bar(stat = "identity", fill = "#8f4e46", color = "black", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Tissue type Frequency", x = "Tissue types", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Explore missing values
miss_var_summary(gdsc)
gg_miss_upset(gdsc)

# Reducing dataset
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

# Exploring normality
shapiro.test(sample(gdsc$LN_IC50, 5000))

library(moments)
skewness(gdsc_final$LN_IC50)
kurtosis(gdsc_final$LN_IC50)

qqnorm(gdsc$LN_IC50)
qqline(gdsc$LN_IC50)

# Exploration of drugs and cell lines
gdsc_final %>%
  group_by(tissue_descriptor_2) %>%
  summarise(mean = mean(LN_IC50), median(LN_IC50), sd(LN_IC50), n(), max(LN_IC50), min(LN_IC50)) %>%
  arrange(mean)

# Drugs tested twice by cell lines
gdsc_final %>% count(DRUG_NAME, CELL_LINE_NAME, sort = TRUE) %>% filter(n > 1, n < 3) %>% distinct(DRUG_NAME)

# Number of tests by cell line
gdsc_final %>% count(CELL_LINE_NAME, CELL_LINE_NAME, sort = TRUE)
gdsc_final %>% count(CELL_LINE_NAME, CELL_LINE_NAME, sort = TRUE) %>%
  ggplot(aes(x = n)) +
  geom_histogram(binwidth = 20, color="black", alpha=0.6, fill="darkblue") +
  theme_minimal()

# Exploration of drug sensitivity by drug, tissue type and cell line
# Function for exploring IC50 Z-score for given cell line, with example
Graph_IC50_ZSCORE_by_Cell_Line <- function(dataframe, cell_line){
  dataframe %>% filter(CELL_LINE_NAME == cell_line) %>%
    ggplot(aes(x = reorder(DRUG_NAME, Z_SCORE), y = Z_SCORE)) +
    geom_point() +
    labs(title = sprintf("Drug sensitivity by cell line %s", cell_line), x = cell_line, y = "Cell line IC50 Z-Score") +
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
}
Graph_IC50_ZSCORE_by_Cell_Line(gdsc_final, "MFM-223")

# Function for exploring LN_IC50 for given drug across all cell lines, with example
Graph_IC50_by_Drug <- function(dataframe, drug_name){
  dataframe %>% filter(DRUG_NAME == drug_name) %>%
    ggplot(aes(x = reorder(CELL_LINE_NAME, LN_IC50), y = LN_IC50)) +
    geom_point(aes(color = tissue_descriptor_2), size = 1.5) +
    theme_classic() +
    labs(title = sprintf("IC50 by compound %s", drug_name), x = " ", y = "LN_IC50") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
}
Graph_IC50_by_Drug(gdsc_final, "Dactinomycin")

# Distribution of drug response by tissue type
gdsc_final %>%
  ggplot(aes(x = tissue_descriptor_2, y = LN_IC50, fill = tissue_descriptor_2)) +
  geom_boxplot(outliers = F) +
  theme_minimal() +
  labs(title = "LN_IC50 distribution by tissue type", x = "Tissue type", y = "LN IC50") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")

# Function for exploring a given drug's response across tissue types, with example
Graph_LNIC50_by_TissueType <- function(dataframe, drug_name){
  dataframe %>% filter(DRUG_NAME == drug_name) %>%
    ggplot(aes(x = tissue_descriptor_2, y = LN_IC50, fill = tissue_descriptor_2)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = sprintf("Sensitivity to %s by tissue type", drug_name), x = "Tissue types", y = "LN_IC50") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")
}
Graph_LNIC50_by_TissueType(gdsc_final, "Crizotinib")

# Function for exploring pathway sensitivity for given cell line, with example
Graph_pathway_by_CellLine <- function(dataframe, cell_line){
  dataframe %>% filter(CELL_LINE_NAME == cell_line) %>%
    ggplot(aes(x = TARGET_PATHWAY, y = LN_IC50, fill = TARGET_PATHWAY)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = sprintf("Pathway sensitivity by cell line %s", cell_line), x = "Target Pathways", y = "LN_IC50") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")
}
Graph_pathway_by_CellLine(gdsc_final, "MFM-223")

# Correlation analysis. Spearman correlation. Scatter plot between LN_IC50 and AUC.
round(cor(gdsc_final$AUC, gdsc_final$LN_IC50, method = "spearman"), 3)
ggplot(gdsc_final, aes(x = LN_IC50, y = AUC)) +
  geom_point(shape = 20, size = 1) +
  labs(title = "LN_IC50 vs. AUC",
       x = "LN_IC50",
       y = "AUC") +
  theme_classic()

# Function for exploring spearman correlation between LN_IC50 and AUC for specific tissue, with examples.
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

#Statistical tests.
# Binary variables.
wilcoxTest_MSI <- wilcox.test(gdsc_final$LN_IC50[gdsc_final$msi == 1],
                              gdsc_final$LN_IC50[gdsc_final$msi == 0],
                              var.equal = T)
wilcoxTest_MSI

wilcoxTest_CNA <- wilcox.test(gdsc_final$LN_IC50[gdsc_final$cna == 1],
                              gdsc_final$LN_IC50[gdsc_final$cna == 0], 
                              var.equal = T)
wilcoxTest_CNA

wilcoxTest_geneExpression <- wilcox.test(gdsc_final$LN_IC50[gdsc_final$gene_expression == 1],
                                         gdsc_final$LN_IC50[gdsc_final$gene_expression == 0],
                                         var.equal = T)
wilcoxTest_geneExpression

wilcoxTest_methylation <- wilcox.test(gdsc_final$LN_IC50[gdsc_final$Methylation == 1],
                                      gdsc_final$LN_IC50[gdsc_final$Methylation == 0],
                                      var.equal = T)
wilcoxTest_methylation

#Categorical variables.
library(FSA, include.only = "dunnTest")
#Kruskal test by target pathway
gdsc_final$TARGET_PATHWAY <- as.factor(gdsc_final$TARGET_PATHWAY)
targetPathway_kt <- kruskal.test(LN_IC50 ~ TARGET_PATHWAY, data = gdsc_final)
targetPathway_kt
targetPathway_dunnt <- dunnTest(LN_IC50 ~ TARGET_PATHWAY, data=gdsc_final,
                                method = "bonferroni")
targetPathway_dunnt

# Kruskal test by tissue type
gdsc_final$tissue_descriptor_2 <- as.factor(gdsc_final$tissue_descriptor_2)
tissueType_kt <- kruskal.test(LN_IC50 ~ tissue_descriptor_2, data = gdsc_final)
tissueType_kt
tissueType_dunnt <- dunnTest(LN_IC50 ~ tissue_descriptor_2, data = gdsc_final,
                             method = "bonferroni")
tissueType_dunnt

# Function for exploring dunn test by given tissue
targetPathway_by_tissue_dunnt <- function(dataframe, tissue) {
  dataframe$TARGET_PATHWAY <- as.factor(dataframe$TARGET_PATHWAY)
  tissue_df <- dataframe |>
    filter(tissue_descriptor_2 == tissue)
  dunn_test <- dunnTest(LN_IC50 ~ TARGET_PATHWAY, data = tissue_df)
  dunn_test
}

# Kruskal test by target
gdsc_final$TARGET <- as.factor(gdsc_final$TARGET)
target_kt <- kruskal.test(LN_IC50 ~ TARGET, data = gdsc_final)
target_kt
target_dunnt <- dunnTest(LN_IC50 ~ TARGET, data = df)
target_dunnt

# Modeling data. Multiple linear regression.
mlr_model <- lm(LN_IC50 ~ CELL_LINE_NAME + DRUG_NAME,
                data = gdsc_final)
summary(mlr_model)

# Other model tested
lm(LN_IC50 ~ TARGET + tissue_descriptor_2, data = gdsc_final)

# Evaluation of the model
plot(mlr_model)

hist(standard_residuals,
     main = "Histogram of Residuals",
     xlab = "Residuals")

# Logistic regression model
library(caret)
library(pROC)
median_IC50 <- median(gdsc_final$LN_IC50, na.rm = TRUE)
gdsc_final$resistance_label <- ifelse(gdsc_final$LN_IC50 > median_IC50, 1, 0)

set.seed(21)
trainIndex <- createDataPartition(gdsc_final$LN_IC50, p = 0.8, list = FALSE)

train_data <- gdsc_final[trainIndex, ]
test_data  <- gdsc_final[-trainIndex, ]

lr_model <- glm(
  resistance_label ~ AUC + tissue_descriptor_2 + TARGET,
  data = train_data |> select(-LN_IC50),
  family = "binomial")

logpred_probs <- predict(lr_model, newdata = test_data, type = "response")
logpred_class <- ifelse(logpred_probs > 0.5, 1, 0)
confusionMatrix(factor(logpred_class), factor(test_data$resistance_label))

roc_curve_log <- roc(test_data$resistance_label, logpred_probs)
plot(roc_curve_log)

# XGBoost regression model
library(xgboost)

train_matrix <- model.matrix(LN_IC50 ~ CELL_LINE_NAME + DRUG_NAME, 
                             data = train_data)
test_matrix  <- model.matrix(LN_IC50 ~ CELL_LINE_NAME + DRUG_NAME, 
                             data = test_data)

dtrain <- xgb.DMatrix(data = train_matrix, label = train_data$LN_IC50)
dtest  <- xgb.DMatrix(data = test_matrix, label = test_data$LN_IC50)

xgb_model <- xgboost(data = dtrain, 
                     objective = "reg:squarederror",
                     nrounds = 100,
                     eta = 0.1)

xgb_predictions <- predict(xgb_model, dtest)
reg_results <- data.frame(Real = test_data$LN_IC50, Predicted = xgb_predictions)
head(reg_results)

mse <- mean((test_data$LN_IC50 - xgb_predictions)^2)
rmse <- sqrt(mse)
cat("Mean Squared Error:", mse, "\n")
cat("Root Mean Squared Error:", rmse, "\n")

importance <- xgb.importance(model = xgb_model)
xgb.plot.importance(importance)

# XGBoost classifier
X <- model.matrix(~ AUC + tissue_descriptor_2 + TARGET, data = train_data)[, -1]
train_matrix2 <- xgb.DMatrix(data = X, label = train_data$resistance_label)
Y <- model.matrix(~ AUC + tissue_descriptor_2 + TARGET, data = test_data)[, -1]
test_matrix2 <- xgb.DMatrix(data = Y, label = test_data$resistance_label)

params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  max_depth = 6,
  eta = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8)

xgb_classifier <- xgb.train(params = params,
                             data = train_matrix2,
                             nrounds = 100,
                             watchlist = list(eval = test_matrix2, train = train_matrix2),
                             early_stopping_rounds = 10)

xgb_prob <- predict(xgb_model_class, newdata = test_matrix2)
xgb_classifier_predictions <- ifelse(xgb_prob > 0.5, 1, 0)
confusionMatrix(factor(xgb_classifier_predictions), factor(test_data$resistance_label))

rocCurve_xgbClassifier <- roc(test_data$resistance_label, xgb_prob)
plot(rocCurve_xgbClassifier)

importance_matrix <- xgb.importance(model = xgb_classifier)
xgb.plot.importance(importance_matrix)