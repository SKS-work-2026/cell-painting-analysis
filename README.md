# cell-painting-analysis
Automated R pipeline for Cell Painting, MOA analysis, and dose-response modeling

# 🔬 Cell Painting Mitochondrial Analysis Pipeline

## 📌 Overview
This repository contains a fully automated, reproducible R-based pipeline for analysis of **Cell Painting data** generated using Harmony software (Revvity/PerkinElmer).

The workflow integrates:
- Feature preprocessing and normalization
- Mitochondrial toxicity scoring
- Mechanism of Action (MOA) clustering
- Dose-response analysis
- High-quality visualization outputs

The pipeline is designed to be **modular, scalable, and publication-ready**, supporting high-content screening datasets.

## 📁 Project Structure
cell-painting-analysis/
│
├── data/
│   ├── raw/              # Input Excel files (not committed to Git)
│   ├── processed/        # Intermediate RDS objects
│
├── scripts/
│   ├── 01_load_data.R
│   ├── 02_preprocess.R
│   ├── 03_mito_scoring.R
│   ├── 04_analysis_plots.R
│
├── config/
│   └── parameters.R      # Central configuration file
│
├── results/
│   ├── plots/            # Generated figures
│   ├── tables/           # Output tables
│
├── run_pipeline.R        # Master script
├── README.md

# Methods
Cell Painting Data Processing and Feature Standardization
High-content Cell Painting datasets were generated using the Harmony image analysis platform (Revvity/PerkinElmer) and exported in tabular format. All downstream analyses were performed in R (v4.2 or higher) using a reproducible, script-based workflow.
To ensure robust and automated data ingestion, the header row was programmatically identified by detecting rows containing both positional identifiers (“Row”, “Column”). The dataset was subsequently imported with standardized column naming to ensure compatibility across analytical steps.
Metadata fields, including well position, timepoint, compound identity, concentration, and cell count, were separated from quantitative feature measurements. Only numeric features were retained for downstream analysis.

Quality Control and Feature Filtering
To mitigate technical noise and improve robustness of downstream analyses, feature filtering was performed using the following criteria:


Missing data filtering
Features with >30% missing values were excluded. Remaining missing entries were imputed using median substitution, preserving distributional properties while minimizing bias.


Variance filtering
Near-zero variance features were removed using the caret framework to eliminate uninformative or invariant descriptors.


Feature scaling
All features were standardized using Z-score normalization to ensure comparability:
Z=X−μσZ = \frac{X - \mu}{\sigma}Z=σX−μ​


Control-based normalization
Vehicle-treated samples (e.g., DMSO) were used as internal controls. Feature values were centered relative to the mean control profile to generate a perturbation signature:
Si=Zi−ZˉcontrolS_i = Z_i - \bar{Z}_{control}Si​=Zi​−Zˉcontrol​


This normalization strategy attenuates batch effects and emphasizes compound-induced phenotypic deviations.

Definition of Mitochondrial Feature Space
Mitochondrial-specific features were identified by pattern matching descriptors containing “MitoTracker” or “mito”. This subset captures mitochondrial morphology, texture, and intensity-based measurements, enabling targeted interrogation of mitochondrial perturbations.

Quantification of Mitochondrial Perturbation
Mitochondrial Toxicity Score
Mitochondrial perturbation was quantified using the Euclidean distance between each sample and the control centroid within the mitochondrial feature space:
MitoScorei=∑j=1n(xij−μcontrol,j)2\text{MitoScore}_i = \sqrt{\sum_{j=1}^{n} (x_{ij} - \mu_{control,j})^2}MitoScorei​=j=1∑n​(xij​−μcontrol,j​)2​
where xijx_{ij}xij​ denotes the value of feature j in sample i, and μcontrol,j\mu_{control,j}μcontrol,j​ represents the corresponding mean value in control samples.
This metric provides an aggregate measure of deviation from baseline mitochondrial state, with higher values indicating stronger perturbation.

Phenotypic Feature Decomposition
To disentangle mechanistically distinct mitochondrial phenotypes, two composite scores were derived:
Depolarization Score
Calculated as the mean of intensity-related features, including mean, median, and integrated signal intensity:
DepolarScorei=1n∑j∈intensityxij\text{DepolarScore}_i = \frac{1}{n}\sum_{j \in \text{intensity}} x_{ij}DepolarScorei​=n1​j∈intensity∑​xij​
This score reflects alterations in mitochondrial membrane potential.

Fragmentation Score
Derived from structural and texture descriptors, including Haralick features, SER-based metrics, and ridge/spot detection features:
FragmentScorei=1m∑j∈structurexij\text{FragmentScore}_i = \frac{1}{m}\sum_{j \in \text{structure}} x_{ij}FragmentScorei​=m1​j∈structure∑​xij​
This metric captures morphological disruption, including mitochondrial fragmentation.

Phenotypic Classification of Mitochondrial States
Samples were classified into discrete mitochondrial phenotypes using empirically defined thresholds:


Depolarization phenotype:
DepolarScore<−1 ∧ FragmentScore≈0\text{DepolarScore} < -1 \ \wedge \ \text{FragmentScore} \approx 0DepolarScore<−1 ∧ FragmentScore≈0


Fragmentation phenotype:
FragmentScore>1\text{FragmentScore} > 1FragmentScore>1


Mixed phenotype:
Intermediate or overlapping profiles


This rule-based framework enables mechanistic stratification of compound-induced mitochondrial effects.

Dimensionality Reduction and Phenotypic Mapping
To visualize phenotypic structure in reduced-dimensional space:


Principal Component Analysis (PCA) was employed to capture major axes of variance in mitochondrial features.


Uniform Manifold Approximation and Projection (UMAP) was used for non-linear embedding, preserving both local and global relationships between samples.


These approaches facilitate identification of clustering patterns and phenotypic relationships across compounds.

Dose–Response Profiling
Dose-dependent effects were assessed by modeling mitochondrial perturbation (MitoScore) as a function of compound concentration on a logarithmic scale. Only non-zero concentrations were included. This analysis enables detection of concentration-dependent toxicity trends and compound potency.

Feature-Level Visualization and Heatmap Clustering
To interrogate feature-level perturbations, the top 30 most variable mitochondrial features were selected based on variance ranking. A heatmap was generated following row-wise scaling, enabling comparative visualization of relative feature perturbations across samples.
Hierarchical clustering was applied to identify patterns of co-regulated features and sample similarity.

Mechanism of Action (MOA) Inference
To infer mechanistic similarity among compounds:

Feature profiles were aggregated at the compound level by computing mean feature values
Pairwise Euclidean distances were calculated
Hierarchical clustering was performed using agglomerative methods

This approach groups compounds based on similarity in mitochondrial phenotypes, enabling hypothesis generation regarding shared mechanisms of action.

Compound Toxicity Ranking
Compounds were ranked based on their maximal mitochondrial perturbation:
Max Toxicitycompound=max⁡(MitoScore)\text{Max Toxicity}_{compound} = \max(\text{MitoScore})Max Toxicitycompound​=max(MitoScore)
This metric provides a straightforward prioritization of compounds with the highest mitochondrial liability.

Software and Reproducibility
All analyses were conducted using R and the following packages:
tidyverse, readxl, caret, umap, and pheatmap.
The workflow is fully modular and reproducible, with all parameters—including thresholds, visualization settings, and file paths—defined in a centralized configuration file. Intermediate datasets are stored as serialized R objects (.rds), ensuring computational efficiency and reproducibility.
