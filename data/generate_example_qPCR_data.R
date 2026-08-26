# Generate the example dataset used by the app's "Load Example" button.
#
# Seeded noise keeps the values reproducible while avoiding perfectly parallel
# condition profiles. Each Sample x Target x Replicate group contains three
# technical measurements.

samples <- c("Control", "Treatment_A", "Treatment_B", "Treatment_C")
replicates <- paste0("R", 1:5)
targets <- c(
    "ACTB",
    "TBP",
    "Target_Stable",
    "Target_Induced_A",
    "Target_Repressed_B",
    "Target_Graded",
    "Target_TechCensor",
    "Target_BioCensor",
    "Target_ReferenceCensor",
    "Target_SampleAbsent",
    "Target_TwoSamplesAbsent",
    "Target_AllAbsent",
    "Target_HighCq"
)

example_data <- expand.grid(
    Technical = 1:3,
    Target = targets,
    Sample = samples,
    Replicate = replicates,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)

set.seed(20260827)

# Independent biological-group noise prevents target and housekeeping effects
# from cancelling exactly during normalization. Smaller well-level noise makes
# technical replicates realistic without obscuring the example's main trends.
biological_groups <- unique(
    example_data[c("Sample", "Target", "Replicate")]
)
group_key <- function(sample, target, replicate) {
    paste(sample, target, replicate, sep = "\r")
}
biological_noise <- setNames(
    rnorm(nrow(biological_groups), mean = 0, sd = 0.18),
    group_key(
        biological_groups$Sample,
        biological_groups$Target,
        biological_groups$Replicate
    )
)
technical_noise <- rnorm(nrow(example_data), mean = 0, sd = 0.04)

sample_shift <- c(Control = 0, Treatment_A = 0.30,
                  Treatment_B = -0.20, Treatment_C = 0.15)
replicate_shift <- c(R1 = -0.15, R2 = 0.10, R3 = -0.05,
                     R4 = 0.12, R5 = -0.02)
biological_deviation <- c(R1 = -0.18, R2 = 0.08, R3 = -0.06,
                          R4 = 0.14, R5 = 0.02)
technical_deviation <- c(-0.12, 0, 0.12)

target_dcq <- c(
    Target_Stable = 5.5,
    Target_Induced_A = 7.5,
    Target_Repressed_B = 6,
    Target_Graded = 8,
    Target_TechCensor = 8.5,
    Target_BioCensor = 9,
    Target_ReferenceCensor = 9.5,
    Target_SampleAbsent = 10,
    Target_TwoSamplesAbsent = 10.5,
    Target_AllAbsent = 11
)

numeric_cq <- mapply(function(sample, target, replicate, technical, well_noise) {
    run_center <- 21 + sample_shift[[sample]] + replicate_shift[[replicate]]
    tech_offset <- technical_deviation[[technical]] + well_noise
    group_offset <- biological_noise[[group_key(sample, target, replicate)]]

    if (target == "ACTB") {
        return(run_center - 0.45 + group_offset + tech_offset)
    }
    if (target == "TBP") {
        return(run_center + 0.45 + group_offset + tech_offset)
    }

    if (target == "Target_HighCq") {
        return(38.8 + sample_shift[[sample]] +
                   biological_deviation[[replicate]] + group_offset +
                   tech_offset)
    }

    dcq <- target_dcq[[target]]
    if (target == "Target_Induced_A" && sample == "Treatment_A") dcq <- 4.5
    if (target == "Target_Repressed_B" && sample == "Treatment_B") dcq <- 9
    if (target == "Target_Graded") {
        dcq <- c(Control = 8, Treatment_A = 7,
                 Treatment_B = 6, Treatment_C = 5)[[sample]]
    }

    run_center + dcq + biological_deviation[[replicate]] + group_offset +
        tech_offset
}, example_data$Sample, example_data$Target,
   example_data$Replicate, example_data$Technical, technical_noise)

example_data$Cq <- sprintf("%.2f", numeric_cq)

is_case <- function(target, sample = NULL, replicate = NULL, technical = NULL) {
    selected <- example_data$Target == target
    if (!is.null(sample)) selected <- selected & example_data$Sample %in% sample
    if (!is.null(replicate)) selected <- selected & example_data$Replicate %in% replicate
    if (!is.null(technical)) selected <- selected & example_data$Technical %in% technical
    selected
}

# One undetected technical measurement; detected measurements must take priority.
example_data$Cq[is_case("Target_TechCensor", "Treatment_A", "R2", 3)] <- "No Ct"

# One fully undetected biological replicate in a non-reference condition.
example_data$Cq[is_case("Target_BioCensor", "Treatment_B", "R3")] <- "Undetermined"

# One fully undetected biological replicate in the default reference sample.
example_data$Cq[is_case("Target_ReferenceCensor", "Control", "R4")] <- ">40"

# A target absent from one whole sample, and another absent from two samples.
example_data$Cq[is_case("Target_SampleAbsent", "Treatment_C")] <- "Undetermined"
example_data$Cq[is_case(
    "Target_TwoSamplesAbsent",
    c("Treatment_B", "Treatment_C")
)] <- "No Ct"

# A target that is undetected in every sample and biological replicate.
example_data$Cq[is_case("Target_AllAbsent")] <- "Undetermined"

# A partial HK technical failure remains usable because other wells are detected.
example_data$Cq[is_case("ACTB", "Treatment_A", "R2", 1)] <- ">40"

# A complete failure of one HK gene invalidates only this biological sample.
example_data$Cq[is_case("TBP", "Treatment_B", "R4")] <- "Undetermined"

# Exercise automatic max-cycle adjustment: the replacement should become 41.
example_data$Cq[is_case("Target_HighCq", "Treatment_C", "R5", 3)] <- "40.40"

example_data <- example_data[c("Sample", "Target", "Cq", "Replicate")]

readr::write_csv(example_data, "data/example_qPCR_data.csv")
