library(tidyverse)
library(caret)

source("config/parameters.R")

data <- readRDS("data/processed/data_clean.rds")

meta_cols <- c("Row","Column","Timepoint","Sequence",
               "Compound","Concentration","Cell.Type","Cell.Count")

meta_cols <- meta_cols[meta_cols %in% colnames(data)]

metadata <- data[, meta_cols]

features <- data %>%
  select(where(is.numeric)) %>%
  select(-any_of(meta_cols))

# CLEANING
features[is.nan(as.matrix(features))] <- NA
features <- features[, colSums(is.na(features)) < na_threshold*nrow(features)]

features <- features %>%
  mutate(across(everything(),
                ~ ifelse(is.na(.), median(., na.rm=TRUE), .)))

features <- features[, -nearZeroVar(features)]

# NORMALIZATION
features_scaled <- scale(features)

metadata$Compound <- as.factor(metadata$Compound)

control_name <- unique(metadata$Compound)[
  grep("vehicle|dmso", unique(metadata$Compound), ignore.case=TRUE)
][1]

control_idx <- which(metadata$Compound == control_name)

control_mean <- colMeans(features_scaled[control_idx, ])
signature <- sweep(features_scaled, 2, control_mean, "-")

saveRDS(list(metadata=metadata, signature=signature),
        "data/processed/processed_data.rds")

cat("✅ Preprocessing complete\n")
