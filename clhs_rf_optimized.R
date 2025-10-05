# ============================================================================
# SOIL SAMPLING OPTIMIZATION: cLHS + Random Forest
# ============================================================================
# This script optimizes spatial sampling designs by combining:
#   1. Conditioned Latin Hypercube Sampling (cLHS) for initial representative samples
#   2. Random Forest optimization with simulated annealing for improved designs
# ============================================================================

library(terra)         # Spatial raster data handling
library(clhs)          # Conditioned Latin Hypercube Sampling
library(randomForest)  # Random Forest modeling for optimization
library(ggplot2)       # Visualization
library(gridExtra)     # Plot arrangement

# ============================================================================
# STEP 1: DATA PREPARATION FUNCTIONS
# ============================================================================

# Converts spatial raster data into a clean dataframe
# Removes missing values and extracts XY coordinates
prepare_covariate_dataframe <- function(covariates) {
  # Convert raster to dataframe with coordinates
  cov_df <- as.data.frame(covariates, xy = TRUE, na.rm = TRUE)
  # Remove any rows with missing data
  na.omit(cov_df)
}

# Prepares covariate data specifically for cLHS algorithm
# Removes coordinates, validates columns, and filters low-variance variables
build_clhs_data <- function(cov_df) {
  # Remove X and Y coordinates (columns 1-2) - cLHS only needs covariate values
  clhs_data <- cov_df[, -c(1, 2), drop = FALSE]

  # Ensure valid column names for R
  names(clhs_data) <- make.names(names(clhs_data))

  # Convert all columns to numeric
  clhs_data <- as.data.frame(lapply(clhs_data, function(x) as.numeric(as.character(x))))

  # Identify valid columns: must have variance > threshold
  # Low variance columns provide no information for sampling
  valid_cols <- vapply(clhs_data, function(x) {
    if (all(is.na(x))) {
      return(FALSE)
    }
    var_val <- stats::var(x, na.rm = TRUE)
    !is.na(var_val) && var_val > 1e-10
  }, logical(1))

  # Return only valid columns
  clhs_data[, valid_cols, drop = FALSE]
}

# ============================================================================
# STEP 2: cLHS SAMPLING FUNCTIONS
# ============================================================================

# Performs Conditioned Latin Hypercube Sampling to select representative points
# cLHS ensures samples cover the full range of environmental conditions
select_clhs_indices <- function(cov_df, n_samples, seed = NULL) {
  cat("Starting cLHS sampling...\n")

  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
    cat("cLHS seed set to", seed, "\n")
  }

  # Prepare data for cLHS (remove coords, validate)
  clhs_data <- build_clhs_data(cov_df)

  # Check if we have valid covariates
  if (ncol(clhs_data) == 0) {
    stop("No valid covariates available for cLHS.")
  }
  cat("Using", ncol(clhs_data), "covariates for cLHS\n")

  # Run cLHS algorithm - returns row indices of selected samples
  clhs(clhs_data, size = n_samples, progress = FALSE)
}

# Main wrapper for cLHS sampling workflow
# Returns dataframe with selected sample locations and covariate values
simple_clhs_sampling <- function(covariates, n_samples, seed = NULL) {
  # Step 2.1: Convert raster to dataframe
  cov_df <- prepare_covariate_dataframe(covariates)

  # Step 2.2: Select sample indices using cLHS
  indices <- select_clhs_indices(cov_df, n_samples, seed)

  # Step 2.3: Return selected rows (includes X, Y, and all covariates)
  cov_df[indices, , drop = FALSE]
}

# ============================================================================
# STEP 3: RANDOM FOREST OPTIMIZATION FUNCTIONS
# ============================================================================

# Calculates Random Forest Mean Squared Error using cross-validation
# This MSE is used as the objective function for optimization
# Lower MSE = better sampling design (samples are more informative)
calculate_rf_mse <- function(sample_points, full_data = NULL, n_folds = 5, target_var = NULL) {
  # Step 3.1: Remove coordinates (X, Y) - only use covariate values
  sample_data <- sample_points[, -c(1, 2), drop = FALSE]

  # Step 3.2: Align columns with full dataset if provided
  if (!is.null(full_data)) {
    shared_cols <- intersect(names(sample_data), names(full_data))
    if (length(shared_cols) > 0) {
      sample_data <- sample_data[, shared_cols, drop = FALSE]
    }
  }

  # Step 3.3: Remove columns with all NA values
  sample_data <- sample_data[, vapply(sample_data, function(col) !all(is.na(col)), logical(1)), drop = FALSE]

  # Step 3.4: Validate data dimensions
  if (nrow(sample_data) < 2 || ncol(sample_data) < 2) {
    return(NA_real_)
  }

  # Step 3.5: Determine response variable (what we're trying to predict)
  if (!is.null(target_var) && target_var %in% names(sample_data)) {
    response <- target_var
  } else if (ncol(sample_data) >= 3) {
    response <- names(sample_data)[3]  # Default: third covariate
  } else {
    response <- names(sample_data)[1]  # Fallback: first covariate
  }

  # Step 3.6: Define predictor variables (all except response)
  predictors <- setdiff(names(sample_data), response)
  if (length(predictors) == 0) {
    return(NA_real_)
  }

  # Step 3.7: Setup k-fold cross-validation
  n_samples <- nrow(sample_data)
  n_folds <- max(1, min(n_folds, n_samples))  # Ensure valid fold count
  # Randomly assign each sample to a fold
  fold_ids <- sample(rep(seq_len(n_folds), length.out = n_samples))

  # Step 3.8: Perform cross-validation
  mse_values <- numeric(n_folds)
  for (fold in seq_len(n_folds)) {
    # Split data into train/test for this fold
    test_indices <- which(fold_ids == fold)
    train_indices <- setdiff(seq_len(n_samples), test_indices)

    # Skip invalid splits
    if (length(test_indices) == 0 || length(train_indices) == 0) {
      next
    }

    train_data <- sample_data[train_indices, , drop = FALSE]
    test_data <- sample_data[test_indices, , drop = FALSE]

    # Skip if data has missing values
    if (anyNA(train_data) || anyNA(test_data)) {
      next
    }

    # Train Random Forest model
    rf_model <- randomForest(stats::as.formula(paste(response, "~", paste(predictors, collapse = "+"))),
                             data = train_data, ntree = 100)

    # Predict on test set and calculate MSE
    predictions <- predict(rf_model, test_data)
    mse_values[fold] <- mean((test_data[[response]] - predictions)^2, na.rm = TRUE)
  }

  # Step 3.9: Calculate average MSE across all folds
  mse_values <- mse_values[is.finite(mse_values)]
  if (length(mse_values) == 0) {
    return(NA_real_)
  }
  mean(mse_values)
}

# Optimizes sample locations using simulated annealing with RF MSE
# Iteratively swaps sample points to minimize prediction error
rf_optimize_sampling <- function(covariates, n_samples, n_iterations = 500,
                                 seed = NULL, temperature = 1000,
                                 cooling_rate = 0.95, target_var = NULL) {
  cat("Starting RF optimization...\n")

  # Set seed for reproducibility (different from cLHS seed)
  if (!is.null(seed)) {
    set.seed(seed + 1)
    cat("RF optimization seed set to", seed + 1, "\n")
  }

  # Step 3.10: Initialize with cLHS samples
  cov_df <- prepare_covariate_dataframe(covariates)
  initial_indices <- select_clhs_indices(cov_df, n_samples, seed)
  current_sample <- cov_df[initial_indices, , drop = FALSE]
  best_sample <- current_sample

  # Step 3.11: Calculate initial MSE (baseline to improve upon)
  current_mse <- calculate_rf_mse(current_sample, cov_df,
                                  target_var = target_var)
  best_mse <- current_mse
  initial_mse <- current_mse
  cat("Initial MSE:", round(current_mse, 4), "\n")

  # Step 3.12: SIMULATED ANNEALING OPTIMIZATION LOOP
  for (i in seq_len(n_iterations)) {
    # Step 3.12.1: Create candidate solution by swapping one point
    candidate_sample <- current_sample
    replace_idx <- sample(nrow(candidate_sample), 1)  # Point to replace

    # Find all points NOT currently in the sample
    occupied_indices <- match(paste(candidate_sample$x, candidate_sample$y),
                              paste(cov_df$x, cov_df$y))
    occupied_indices <- occupied_indices[!is.na(occupied_indices)]
    remaining_points <- setdiff(seq_len(nrow(cov_df)), occupied_indices)

    if (length(remaining_points) == 0) {
      break  # No more points available to swap
    }

    # Randomly select a new point from remaining pool
    new_point_idx <- sample(remaining_points, 1)
    candidate_sample[replace_idx, ] <- cov_df[new_point_idx, ]

    # Step 3.12.2: Evaluate candidate solution
    candidate_mse <- calculate_rf_mse(candidate_sample, cov_df,
                                      target_var = target_var)

    # Step 3.12.3: Acceptance criterion (Metropolis algorithm)
    delta <- candidate_mse - current_mse
    # Accept if better (delta < 0) OR probabilistically if worse
    if (delta < 0 || runif(1) < exp(-delta / temperature)) {
      current_sample <- candidate_sample
      current_mse <- candidate_mse

      # Update best solution if this is the best so far
      if (!is.na(current_mse) && (is.na(best_mse) ||
                                    current_mse < best_mse)) {
        best_sample <- current_sample
        best_mse <- current_mse
        cat("Iteration", i, "- New best MSE:", round(best_mse, 4), "\n")
      }
    }

    # Step 3.12.4: Cool down temperature (reduces chance of accepting worse)
    temperature <- temperature * cooling_rate

    # Progress update every 100 iterations
    if (i %% 100 == 0) {
      cat("Iteration", i, "- Current MSE:", round(current_mse, 4),
          "- Best MSE:", round(best_mse, 4), "\n")
    }
  }

  # Step 3.13: Calculate final improvement percentage
  improvement <- if (is.na(initial_mse) || is.na(best_mse)) {
    NA_real_
  } else {
    ((initial_mse - best_mse) / initial_mse) * 100
  }
  cat("Optimization completed. Improvement:", round(improvement, 2), "%\n")

  # Return results including both initial and optimized samples
  list(
    optimized_samples = best_sample,
    initial_samples = cov_df[initial_indices, , drop = FALSE],
    initial_mse = initial_mse,
    final_mse = best_mse,
    improvement = improvement
  )
}

# ============================================================================
# STEP 4: VISUALIZATION FUNCTIONS
# ============================================================================

# Creates side-by-side plots comparing cLHS and RF optimized sampling designs
# Overlays sample points on environmental covariate raster
plot_sampling_results <- function(covariates, clhs_samples, rf_samples,
                                  covariate_name = NULL) {
  # Default to 'dem' layer if not specified
  if (is.null(covariate_name)) {
    covariate_name <- "dem"
  }

  # Step 4.1: Extract covariate layer for background
  cov_df <- as.data.frame(covariates[[covariate_name]],
                          xy = TRUE, na.rm = TRUE)
  names(cov_df)[3] <- "value"

  # Step 4.2: Create cLHS plot (red points)
  p1 <- ggplot() +
    geom_raster(data = cov_df, aes(x = x, y = y, fill = value)) +
    geom_point(data = clhs_samples, aes(x = x, y = y),
               color = "red", size = 1.5, shape = 21,
               fill = "white", stroke = 1) +
    scale_fill_gradientn(name = covariate_name,
                         colors = terrain.colors(255)) +
    theme_minimal() +
    labs(title = paste("cLHS Sampling Design (n =",
                       nrow(clhs_samples), ")"),
         x = "X Coordinate", y = "Y Coordinate") +
    theme(aspect.ratio = 1, legend.position = "bottom")

  # Step 4.3: Create RF optimized plot (blue points)
  p2 <- ggplot() +
    geom_raster(data = cov_df, aes(x = x, y = y, fill = value)) +
    geom_point(data = rf_samples, aes(x = x, y = y),
               color = "blue", size = 1.5, shape = 21,
               fill = "white", stroke = 1) +
    scale_fill_gradientn(name = covariate_name,
                         colors = terrain.colors(255)) +
    theme_minimal() +
    labs(title = paste("RF Optimized Design (n =",
                       nrow(rf_samples), ")"),
         x = "X Coordinate", y = "Y Coordinate") +
    theme(aspect.ratio = 1, legend.position = "bottom")

  # Return both plots as a list
  list(clhs_plot = p1, rf_plot = p2)
}

# Creates summary table comparing cLHS vs RF optimization performance
create_comparison_table <- function(rf_results) {
  # Build comparison dataframe with MSE metrics
  comparison_data <- data.frame(
    Method = c("cLHS", "RF Optimized"),
    MSE = c(round(rf_results$initial_mse, 4),
            round(rf_results$final_mse, 4)),
    Samples = c(nrow(rf_results$initial_samples),
                nrow(rf_results$optimized_samples)),
    stringsAsFactors = FALSE
  )

  # Add improvement percentage (0% for baseline cLHS)
  comparison_data$Improvement <- c(0, round(rf_results$improvement, 2))

  comparison_data
}

# ============================================================================
# STEP 5: MAIN WORKFLOW ORCHESTRATION
# ============================================================================

# Main function that runs the complete optimization workflow
# Executes all steps and exports results
run_sampling_optimization <- function(raster_file = "data/predictors.tif",
                                      n_samples = 100,
                                      n_iterations = 500,
                                      output_dir = "outputs/",
                                      seed = 123,
                                      target_var = NULL,
                                      export_results = TRUE) {
  # Initialize random seed for reproducibility
  set.seed(seed)
  cat("Set seed to", seed, "for reproducible results\n")
  cat("=== Soil Sampling Optimization Workflow ===\n")

  # WORKFLOW STEP 1: Load environmental raster data
  cat("Loading environmental data...\n")
  predictors <- rast(raster_file)
  cat("Loaded", nlyr(predictors), "layers:",
      paste(names(predictors), collapse = ", "), "\n\n")

  # WORKFLOW STEP 2: Generate cLHS samples (baseline)
  cat("Step 1: cLHS Sampling\n")
  clhs_samples <- simple_clhs_sampling(predictors, n_samples, seed = seed)
  cat("Generated", nrow(clhs_samples), "cLHS samples\n\n")

  # WORKFLOW STEP 3: Optimize samples using Random Forest
  cat("Step 2: Random Forest Optimization\n")
  rf_results <- rf_optimize_sampling(predictors, n_samples, n_iterations,
                                     seed = seed, target_var = target_var)
  cat("Generated", nrow(rf_results$optimized_samples),
      "RF optimized samples\n\n")

  # WORKFLOW STEP 4: Create visualizations
  cat("Step 3: Creating visualizations...\n")
  plots <- plot_sampling_results(predictors, clhs_samples,
                                 rf_results$optimized_samples)

  # WORKFLOW STEP 5: Generate comparison table
  cat("Step 4: Creating comparison table...\n")
  comparison_table <- create_comparison_table(rf_results)

  # WORKFLOW STEP 6: Export results (optional)
  if (export_results) {
    cat("Step 5: Exporting results...\n")

    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    # Export CSV files
    write.csv(clhs_samples,
              file.path(output_dir, "clhs_sample_locations.csv"),
              row.names = FALSE)
    write.csv(rf_results$optimized_samples,
              file.path(output_dir, "rf_optimized_locations.csv"),
              row.names = FALSE)
    write.csv(comparison_table,
              file.path(output_dir, "sampling_comparison_table.csv"),
              row.names = FALSE)

    cat("✓ Saved clhs_sample_locations.csv\n")
    cat("✓ Saved rf_optimized_locations.csv\n")
    cat("✓ Saved sampling_comparison_table.csv\n")

    # Export combined plot
    combined_plot <- grid.arrange(plots$clhs_plot, plots$rf_plot, ncol = 2)
    ggsave(file.path(output_dir, "sampling_comparison_plots.png"),
           combined_plot, width = 12, height = 6, dpi = 300)
    cat("✓ Saved sampling_comparison_plots.png\n\n")
  }

  # WORKFLOW STEP 7: Display results summary
  cat("=== RESULTS SUMMARY ===\n")
  print(comparison_table)
  cat("\nImprovement from RF optimization:",
      round(rf_results$improvement, 2), "%\n")

  # Return all results as a list
  list(
    clhs_samples = clhs_samples,
    rf_optimized_samples = rf_results$optimized_samples,
    comparison_table = comparison_table,
    plots = plots,
    improvement = rf_results$improvement
  )
}

# ============================================================================
# SCRIPT INITIALIZATION MESSAGE
# ============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("  SOIL SAMPLING OPTIMIZATION: cLHS + Random Forest\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("All functions loaded successfully!\n\n")
cat("USAGE:\n")
cat("  results <- run_sampling_optimization(\n")
cat("    raster_file = 'data/predictors.tif',\n")
cat("    n_samples = 100,\n")
cat("    n_iterations = 500,\n")
cat("    seed = 123\n")
cat("  )\n\n")
cat("WORKFLOW:\n")
cat("  Step 1: Data Preparation\n")
cat("  Step 2: cLHS Sampling (baseline)\n")
cat("  Step 3: Random Forest Optimization (simulated annealing)\n")
cat("  Step 4: Visualization\n")
cat("  Step 5: Export Results\n")
cat("═══════════════════════════════════════════════════════════════════\n")
cat("\n")

# ============================================================================
# AUTOMATIC EXECUTION
# ============================================================================
# Uncomment the line below to run automatically when script is sourced:

results <- run_sampling_optimization(n_samples = 100,
                                     n_iterations = 500,
                                     seed = 123)
