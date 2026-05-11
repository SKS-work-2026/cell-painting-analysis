library(tidyverse)
library(umap)
library(pheatmap)

source("config/parameters.R")

data_list <- readRDS("data/processed/mito_data.rds")

metadata <- data_list$metadata
mito_data <- data_list$mito_data

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ===============================
# PCA
# ===============================
pca <- prcomp(mito_data)
pca_df <- as.data.frame(pca$x) %>% bind_cols(metadata)

pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = Compound)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = compound_colors) +
  theme_minimal(base_size = base_font_size) +
  labs(title = "PCA Plot")

ggsave("results/plots/PCA.png", pca_plot, width = 8, height = 6)

# ===============================
# UMAP
# ===============================
umap_res <- umap(mito_data)
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1","UMAP2")
umap_df <- cbind(metadata, umap_df)

umap_plot <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Compound)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = compound_colors) +
  theme_classic(base_size = base_font_size)

ggsave("results/plots/UMAP.png", umap_plot, width = 8, height = 6)

# ===============================
# HEATMAP
# ===============================
mito_var <- apply(mito_data, 2, var)
top_mito <- names(sort(mito_var, decreasing=TRUE))[1:min(30,length(mito_var))]

annotation_df <- data.frame(
  Compound = metadata$Compound,
  Dose = metadata$Concentration
)

rownames(annotation_df) <- make.unique(as.character(1:nrow(mito_data)))
rownames(mito_data) <- rownames(annotation_df)

png("results/plots/Mito_Heatmap.png",1000,900)

pheatmap(
  mito_data[,top_mito],
  scale="row",
  annotation_row=annotation_df,
  fontsize = 10
)

dev.off()

# ===============================
# DOSE RESPONSE
# ===============================
dose_df <- metadata %>%
  mutate(Concentration = as.numeric(Concentration)) %>%
  filter(Concentration > 0)

dose_plot <- ggplot(dose_df,
                    aes(Concentration, MitoScore, color = Compound)) +
  geom_point(size = 2) +
  geom_line() +
  scale_x_log10() +
  scale_color_brewer(palette = compound_colors) +
  theme_minimal(base_size = base_font_size) +
  labs(title = "Mito Dose Response")

ggsave("results/plots/Mito_Dose_Response.png", dose_plot, width = 8, height = 6)

# ===============================
# PHENOTYPE MAP
# ===============================
pheno_plot <- ggplot(metadata,
                     aes(DepolarScore, FragmentScore, color = MitoClass)) +
  geom_point(size = 4) +
  scale_color_manual(values = mito_colors) +
  theme_classic(base_size = base_font_size) +
  labs(title = "Fragmentation vs Depolarization")

ggsave("results/plots/Mito_Phenotype_Map.png", pheno_plot, width = 7, height = 6)

# ===============================
# MOA CLUSTERING
# ===============================
profiles <- as.data.frame(mito_data) %>%
  cbind(metadata) %>%
  group_by(Compound) %>%
  summarise(across(where(is.numeric), mean), .groups="drop")

png("results/plots/MOA_Dendrogram.png",800,700)

plot(hclust(dist(profiles %>% select(-Compound))),
     main="MOA Clustering")

dev.off()

# ===============================
# TOXICITY RANKING
# ===============================
tox_df <- metadata %>%
  group_by(Compound) %>%
  summarise(Max_Toxicity=max(MitoScore), .groups="drop") %>%
  arrange(desc(Max_Toxicity))

write.csv(tox_df, "results/tables/Toxicity_Ranking.csv", row.names=FALSE)

tox_plot <- ggplot(tox_df,
                   aes(reorder(Compound, Max_Toxicity), Max_Toxicity)) +
  geom_bar(stat="identity", fill = bar_color) +
  coord_flip() +
  theme_classic(base_size = base_font_size) +
  labs(title="Toxicity Ranking")

ggsave("results/plots/Toxicity_Ranking.png", tox_plot, width = 8, height = 6)

cat("✅ ALL PLOTS GENERATED & SAVED\n")
``
