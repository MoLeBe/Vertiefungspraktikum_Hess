# Vertiefungspraktikum Hess

This repository contains R scripts and workflows from my internship project at the Hess Lab.

## 📦 Project Description

This project involves the analysis of gene expression by segmentation of Mircroarray or RNA-seq data using Rcpp and other tools.
The code processes experimental data, performs segmentation using Rcpp, and visualizes the results using ggplot2.

## 🔁 Example Workflow
To begin, the data must be loaded into R. This can be done in one of two ways:

   - If you already have processed data, load it directly as a SummarizedExperiment object.

   - If starting from raw data, use the prepare_replicate_data_for_fit.R pipeline to process it into the correct format.


Regardless of the starting point, you should end up with **two** `SummarizedExperiment` objects—one for each condition.  
In this example, they are named `se_norm_standard` and `se_iron_standard`, representing **standard** and **iron-depleted** conditions, respectively.

### 🧼 Preprocessing

These objects are then passed through a preprocessing step using the `rifi_preprocess_no_TI()` function:

```r
# Example preprocessing (adjust 'thrsh_check' and 'bg' to your dataset)
preprocessed_std <- rifi_preprocess_no_TI(se_norm_standard, cores = 30, thrsh_check = 850, bg = 2500)
preprocessed_iron <- rifi_preprocess_no_TI(se_iron_standard, cores = 30, thrsh_check = 850, bg = 2500)
```

### 🔬 Fitting
After preprocessing, use the nls2_fit_MB() function to perform the non-linear model fitting:

```r
fits_std <- nls2_fit_MB(preprocessed_std)
fits_iron <- nls2_fit_MB(preprocessed_iron)
```