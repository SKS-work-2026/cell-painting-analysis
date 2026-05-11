library(readxl)

source("config/parameters.R")

raw <- read_excel(input_file, col_names = FALSE)

header_row <- which(
  apply(raw, 1, function(x)
    any(grepl("Row", x, ignore.case = TRUE)) &
      any(grepl("Column", x, ignore.case = TRUE))
  )
)[1]

data <- read_excel(input_file, skip = header_row - 1)
colnames(data) <- make.names(colnames(data))

saveRDS(data, "data/processed/data_clean.rds")

cat("✅ Data loaded\n")
