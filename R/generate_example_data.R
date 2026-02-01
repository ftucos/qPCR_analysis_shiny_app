# Generate simulated qPCR example data
# This creates a dataset for testing statistical functions with:
# - 3 samples: ctrl, treatA, treatB
# - 5 biological replicates: R1 to R5
# - 3 technical replicates per measurement
# - 2 housekeeping genes: TBP, Actin
# - 5 target genes with different detection patterns:
#   - Target1: all detected counts
#   - Target2: undetected in all reps for one sample only (treatB)
#   - Target3: undetected in all reps for all samples
#   - Target4: undetected in 2/5 bio reps for ctrl only (R1, R2)
#   - Target5: undetected in all 5 bio reps for ctrl only

library(tidyverse)

set.seed(42)

# Define structure
samples <- c("ctrl", "treatA", "treatB")
replicates <- paste0("R", 1:5)
tech_reps <- 3  # technical replicates per measurement
hk_genes <- c("TBP", "Actin")
target_genes <- c("Target1", "Target2", "Target3", "Target4", "Target5")

# Create base grid (one row per sample + tech rep + bio rep + target )
example_data <- expand.grid(
    Sample = samples,
    Replicate = replicates,
    Target = c(hk_genes, target_genes),
    TechRep = 1:tech_reps,
    stringsAsFactors = FALSE
) |>
    as_tibble() |>
    arrange(Replicate, Sample, Target, TechRep)

# Generate Cq values with realistic variation
# HK genes: stable around 20-22
# Target genes: vary by sample to create differential expression
# Technical replicates have smaller variation than biological replicates
example_data <- example_data |>
    # First create a "true" value per bio replicate, then add tech rep noise
    group_by(Sample, Replicate, Target) |>
    mutate(
        # Base Cq value for this biological replicate
        base_Cq = case_when(
            # Housekeeping genes - stable expression
            Target %in% hk_genes ~ rnorm(1, mean = 21, sd = 0.5),
            
            # Target1: all detected, differential expression between samples
            Target == "Target1" & Sample == "ctrl" ~ rnorm(1, mean = 25, sd = 0.8),
            Target == "Target1" & Sample == "treatA" ~ rnorm(1, mean = 23, sd = 0.8), # upregulated
            Target == "Target1" & Sample == "treatB" ~ rnorm(1, mean = 27, sd = 0.8), # downregulated
            
            # Target2: undetected in treatB only
            Target == "Target2" & Sample == "ctrl" ~ rnorm(1, mean = 28, sd = 1.0),
            Target == "Target2" & Sample == "treatA" ~ rnorm(1, mean = 26, sd = 1.0),
            Target == "Target2" & Sample == "treatB" ~ NA_real_, # undetected
            
            # Target3: undetected in all samples
            Target == "Target3" ~ NA_real_,
            
            # Target4: undetected in 2/5 bio reps for ctrl only (R1, R2)
            Target == "Target4" & Sample == "ctrl" & Replicate %in% c("R1", "R2") ~ NA_real_,
            Target == "Target4" & Sample == "ctrl" ~ rnorm(1, mean = 29, sd = 0.8), # detected in R3-R5
            Target == "Target4" & Sample == "treatA" ~ rnorm(1, mean = 27, sd = 0.8),
            Target == "Target4" & Sample == "treatB" ~ rnorm(1, mean = 28, sd = 0.8),
            
            # Target5: undetected in all 5 bio reps for ctrl only
            Target == "Target5" & Sample == "ctrl" ~ NA_real_,
            Target == "Target5" & Sample == "treatA" ~ rnorm(1, mean = 30, sd = 1.0),
            Target == "Target5" & Sample == "treatB" ~ rnorm(1, mean = 29, sd = 1.0)
        )
    ) |>
    ungroup() |>
    rowwise() |>
    mutate(
        # Add small technical replicate variation (SD ~0.2 Cq)
        Cq = if_else(
            is.na(base_Cq), 
            NA_real_, 
            round(base_Cq + rnorm(1, mean = 0, sd = 0.2), 2)
        )
    ) |>
    ungroup() |>
    # Convert Cq to character for app compatibility (NA becomes "Undetermined")
    mutate(
        Cq = ifelse(is.na(Cq), "Undetermined", as.character(Cq))
    ) |>
    # Drop intermediate columns and reorder
    select(Sample, Target, Cq, Replicate) |>
    # Arrange for nice output
    arrange(Replicate, Sample, Target)

# Save to data folder
write_csv(example_data, "data/example_qPCR_data.csv")

cat("Generated example data with:\n")
cat("- Samples:", paste(samples, collapse = ", "), "\n")
cat("- Biological replicates:", paste(replicates, collapse = ", "), "\n")
cat("- Technical replicates per measurement:", tech_reps, "\n")
cat("- HK genes:", paste(hk_genes, collapse = ", "), "\n")
cat("- Target genes:", paste(target_genes, collapse = ", "), "\n")
cat("- Total rows:", nrow(example_data), "\n")
cat("Saved to data/example_qPCR_data.csv\n")
