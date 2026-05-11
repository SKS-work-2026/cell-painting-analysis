library(tidyverse)

source("config/parameters.R")

data_list <- readRDS("data/processed/processed_data.rds")

metadata <- data_list$metadata
signature <- data_list$signature

# MITO FEATURES
mito_cols <- grep("MitoTracker|mito",
                  colnames(signature), ignore.case=TRUE, value=TRUE)

mito_data <- as.matrix(signature[, mito_cols])

control_name <- unique(metadata$Compound)[
  grep("vehicle|dmso", unique(metadata$Compound), ignore.case=TRUE)
][1]

control_idx <- which(metadata$Compound == control_name)

# MITO SCORE
control_mito <- colMeans(mito_data[control_idx,,drop=FALSE])

metadata$MitoScore <- apply(mito_data, 1, function(x)
  sqrt(sum((x-control_mito)^2)))

# PHENOTYPE
depolar_cols <- grep("Mean|Median|Sum", mito_cols, value=TRUE)
frag_cols <- grep("Haralick|SER|Spot|Ridge|Contrast", mito_cols, value=TRUE)

metadata$DepolarScore <- rowMeans(mito_data[,depolar_cols,drop=FALSE])
metadata$FragmentScore <- rowMeans(mito_data[,frag_cols,drop=FALSE])

metadata$MitoClass <- ifelse(
  metadata$DepolarScore < depolar_threshold & metadata$FragmentScore < 0,
  "Depolarization",
  ifelse(metadata$FragmentScore > fragment_threshold,
         "Fragmentation","Mixed")
)

saveRDS(list(metadata=metadata, mito_data=mito_data),
        "data/processed/mito_data.rds")

cat("✅ Mito scoring complete\n")
