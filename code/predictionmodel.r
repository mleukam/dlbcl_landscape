## -----------------------------------------------------

## R script to make prediction model for cluster from RNAseq data

## Expects cleaned, normalized, log-transformed, expression matrix
##    Originally defined using HTseq reads from GDC
## Expects cluster assignments (here the 4-cluster model from HCPC)
## Outputs a model and a prediction

## M. Leukam
## 7/23/2019

## -----------------------------------------------------

# clear the workspace
rm(list = ls())

# load packages
library("tidyverse")
library("caret")
library("e1071")

# set seed
set.seed(818)

# load data
combined_es <- readRDS("output/combined_es_gsva4_4cluster.rds")

# format data
clust <- pData(combined_es) %>%
  rownames_to_column(var = "sample_id") %>%
  dplyr::select(sample_id, clust) %>%
  as_tibble() %>%
  print()

expr_df <- t(exprs(combined_es)) %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_id") %>%
  as_tibble() %>%
  print()

d_in <- clust %>%
  left_join(expr_df) %>%
  dplyr::select(sample_id, clust, everything()) %>%
  mutate(clust = as.factor(clust)) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_id")

# partition data
set.seed(818)
trainIndex <- createDataPartition(d_in$clust, p = .75, 
                                  list = FALSE, 
                                  times = 1)
d_train <- d_in[ trainIndex,]
d_test  <- d_in[-trainIndex,]

# set parameters for 10-fold CV, repeated 3 times
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 3,
                           verboseIter = TRUE,
                           )

# train model
knn_fit <- train(clust ~ ., 
            data = d_train,
            method = "knn", 
            trControl = fitControl,
            preprocess = c("center", "scale"),
            tunelength = 10)

predictions <- predict(knn_fit, d_test)
xtab <- table(d_test$clust, prediction)

# write out results
saveRDS(d_train, "gpfs/data/kline-lab/tcga_macs/output/train_data.rds")
saveRDS(d_test, "gpfs/data/kline-lab/tcga_macs/output/test_data.rds")
saveRDS(knn_fit, "/gpfs/data/kline-lab/tcga_macs/output/knn_fit.rds")
saveRDS(predictions, "gpfs/data/kline-lab/tcga_macs/output/predictions.rds")

