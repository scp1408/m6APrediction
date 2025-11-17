# m6APrediction: An R Package for m6A Site Prediction Using Random Forest

## Overview
`m6APrediction` is a specialized R package developed to predict N6-methyladenosine (m6A) sites from biological sequence and feature data. Built on a pre-trained random forest model, it streamlines m6A prediction with two key functionalities:
- **Batch prediction**: Efficiently processes multiple samples using a feature data frame as input.
- **Single-sample prediction**: Enables targeted prediction for individual samples by accepting manual feature inputs.

The package returns two critical outputs for each sample: a numeric probability (ranging from 0 to 1) of the site being m6A, and a categorical status ("Positive" or "Negative") determined by a user-defined threshold (default: 0.5).


## Installation
To install `m6APrediction` from GitHub, you will need either the `devtools` or `remotes` R package—both are widely used for installing GitHub-hosted R packages. Follow the steps below:

### Step 1: Install Dependencies (if not already installed)
If you do not have `devtools` or `remotes` installed on your system, run the following code first:
```r
# Install devtools (recommended; alternative: install.packages("remotes"))
install.packages("devtools")
```

### Step 2: Install `m6APrediction` from GitHub
Replace `[Your-GitHub-Username]` with your actual GitHub username and `[Your-Repo-Name]` with the name of your GitHub repository hosting the `m6APrediction` package (e.g., `john-doe/m6APrediction`). Execute one of the following commands:
```r
# Install using devtools
devtools::install_github("[Your-GitHub-Username]/[Your-Repo-Name]")

# OR install using remotes (if you prefer this package)
# remotes::install_github("[Your-GitHub-Username]/[Your-Repo-Name]")
```

> **Important Note**: Ensure your GitHub repository includes all core package files (e.g., `R/` directory with `main.R`, `inst/extdata/` with `rf_fit.rds` and `m6A_input_example.csv`, `DESCRIPTION`, and `NAMESPACE`) to avoid installation errors.


## Minimal Usage Examples
After successful installation, load the package into your R session with:
```r
library(m6APrediction)
```

The package includes a pre-trained random forest model (`rf_fit.rds`) and example feature data (`m6A_input_example.csv`) in the `inst/extdata` directory. Use the `system.file()` function to access these files (no manual path configuration required).


### Example 1: Batch Prediction with `prediction_multiple()`
This example demonstrates how to predict m6A sites for multiple samples using the included example feature data:
```r
# 1. Load the pre-trained model and example feature data from the package
trained_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
example_feature_data <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))

# 2. Run batch prediction (uses default positive threshold = 0.5)
batch_pred_results <- prediction_multiple(ml_fit = trained_model, feature_df = example_feature_data)

# 3. View key results (first 5 rows: input features + prediction outputs)
head(
  batch_pred_results[, c(
    "gc_content", "RNA_type", "RNA_region",
    "predicted_m6A_prob", "predicted_m6A_status"
  )]
)
```


### Example 2: Single-Sample Prediction with `prediction_single()`
This example shows how to predict m6A status for a single sample by manually inputting feature values (ensure values match the allowed formats):
```r
# 1. Load the pre-trained random forest model
trained_model <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))

# 2. Run single-sample prediction with manual feature inputs
single_sample_pred <- prediction_single(
  ml_fit = trained_model,
  gc_content = 0.62,          # Numeric: GC content (0–1 range)
  RNA_type = "mRNA",          # Character: Must be "mRNA", "lincRNA", "lncRNA", or "pseudogene"
  RNA_region = "CDS",         # Character: Must be "CDS", "intron", "3'UTR", or "5'UTR"
  exon_length = 10.8,         # Numeric: Length of the associated exon
  distance_to_junction = 8.9, # Numeric: Distance to the nearest splice junction
  evolutionary_conservation = 0.75, # Numeric: Evolutionary conservation score
  DNA_5mer = "ATCGG",         # Character: 5-nucleotide sequence (A/T/C/G only)
  positive_threshold = 0.5    # Numeric: Threshold for classifying "Positive" (0–1 range)
)

# 3. View the prediction result (named vector: probability + status)
print(single_sample_pred)
```
