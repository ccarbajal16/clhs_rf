# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based soil sampling optimization system that combines **Conditioned Latin Hypercube Sampling (cLHS)** with **Random Forest (RF) optimization**. The workflow generates optimized spatial sampling designs for environmental/soil data collection by:

1. Starting with cLHS to create representative initial samples
2. Using simulated annealing with Random Forest MSE to optimize sample locations
3. Comparing both approaches and exporting results

## Key Commands

### Running the Analysis

```r
# Source the main script
source("clhs_rf_optimized.R")

# Execute with default parameters (100 samples, 500 iterations)
results <- run_sampling_optimization(n_samples = 100, n_iterations = 500, seed = 123)

# Custom execution
results <- run_sampling_optimization(
  raster_file = "data/predictors.tif",
  n_samples = 150,
  n_iterations = 1000,
  output_dir = "outputs/",
  seed = 42,
  target_var = NULL,  # Specify target variable name if needed
  export_results = TRUE
)
```

### Required R Packages

The script automatically loads:
- `terra` - Spatial raster data handling
- `clhs` - Conditioned Latin Hypercube Sampling
- `randomForest` - Random Forest modeling
- `ggplot2` - Visualization
- `gridExtra` - Plot arrangement

## Architecture

### Workflow Pipeline

The system follows a functional programming approach with modular components:

1. **Data Preparation** ([clhs_rf_optimized.R:7-24](clhs_rf_optimized.R#L7-L24))
   - `prepare_covariate_dataframe()` - Converts raster to dataframe
   - `build_clhs_data()` - Cleans and validates covariates for cLHS

2. **cLHS Sampling** ([clhs_rf_optimized.R:26-44](clhs_rf_optimized.R#L26-L44))
   - `select_clhs_indices()` - Performs cLHS selection
   - `simple_clhs_sampling()` - Main wrapper for cLHS workflow

3. **RF Optimization** ([clhs_rf_optimized.R:46-147](clhs_rf_optimized.R#L46-L147))
   - `calculate_rf_mse()` - Cross-validation MSE calculation
   - `rf_optimize_sampling()` - Simulated annealing optimization using RF MSE as objective function

4. **Visualization & Export** ([clhs_rf_optimized.R:149-181](clhs_rf_optimized.R#L149-L181))
   - `plot_sampling_results()` - Creates side-by-side comparison plots
   - `create_comparison_table()` - Generates performance metrics

5. **Main Entry Point** ([clhs_rf_optimized.R:183-231](clhs_rf_optimized.R#L183-L231))
   - `run_sampling_optimization()` - Orchestrates entire workflow

### Optimization Algorithm

The RF optimization uses **simulated annealing** ([clhs_rf_optimized.R:96-147](clhs_rf_optimized.R#L96-L147)):
- Temperature: Initial 1000, cooling rate 0.95
- Each iteration randomly replaces one sample point
- Acceptance criterion: Metropolis algorithm (accepts worse solutions probabilistically)
- Objective: Minimize cross-validated Random Forest MSE

### Data Flow

```
predictors.tif (input raster)
    ↓
prepare_covariate_dataframe()
    ↓
├─→ simple_clhs_sampling() → clhs_samples
│
└─→ rf_optimize_sampling()
       ├─ Initial: select_clhs_indices()
       └─ Optimize: simulated annealing loop
           └─ Evaluate: calculate_rf_mse() (5-fold CV)
    ↓
Comparison & Visualization
    ↓
CSV exports + PNG plots
```

## Input/Output Structure

### Input
- **data/predictors.tif** - Multi-band raster with environmental covariates
  - Each layer represents a different environmental variable
  - Should have no missing data in areas of interest

### Output (in `outputs/` directory)
- `clhs_sample_locations.csv` - Initial cLHS sample coordinates + covariate values
- `rf_optimized_locations.csv` - Optimized sample coordinates + covariate values
- `sampling_comparison_table.csv` - MSE comparison between methods
- `sampling_comparison_plots.png` - Side-by-side visualization

## Reproducibility

All randomization is controlled by the `seed` parameter:
- cLHS uses `seed`
- RF optimization uses `seed + 1`
- Set seed in `run_sampling_optimization()` for reproducible results
