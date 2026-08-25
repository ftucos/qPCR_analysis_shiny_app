# Utilities for replacing censored qPCR values while retaining an explicit
# censoring flag. Numeric replacement values are used for the statistical analyses
# the flag is kept for labels and exports.

# Validate a candidate maximum-cycle value.
#
# Returns one finite, positive number, or NULL when `x` is not suitable. The
# NULL return value allows callers to build explicit fallback chains with %||%.
validate_max_cycle <- function(x) {
    value <- suppressWarnings(as.numeric(x))
    if (length(value) != 1 || !is.finite(value) || value <= 0) return(NULL)
    value
}

# Choose an integer replacement value above all detected measurements.
#
# `censored` may be a scalar or one flag per value in `x`. Censored and
# non-finite values are ignored. The result is at least `minimum` and, when
# measured Cq exceeds the maximum cycle, the first integer strictly above
# the largest detected value.
next_integer_above <- function(x, minimum = 40, censored = FALSE) {
    if (!length(censored) %in% c(1, length(x))) {
        stop("`censored` must have length 1 or the same length as `x`.")
    }

    censored <- rep_len(!is.na(censored) & censored, length(x))
    detected_values <- x[!censored & is.finite(x)]

    if (length(detected_values) == 0) return(as.numeric(minimum))
    if (max(detected_values) <= minimum) return(as.numeric(minimum))
    floor(max(detected_values)) + 1
}

# Replace explicitly censored values with a numeric limit.
#
# `x` and `censored` must have equal lengths. Missing censoring flags are
# treated as FALSE, and uncensored values are returned unchanged.
replace_censored <- function(x, censored, replacement) {
    if (length(censored) != length(x)) {
        stop("`censored` must have the same length as `x`.")
    }

    censored <- !is.na(censored) & censored
    x[censored] <- replacement
    x
}

# Keep detected technical replicates, or censored values when all are censored.
#
# Rows are evaluated within `group_cols`. When a group contains any detected
# non-missing `value_col`, only those detected rows are kept. Otherwise, its
# non-missing censored rows are retained for downstream calculations. The
# returned data frame is ungrouped.
retain_detected_or_all_censored <- function(df,
                                             group_cols = c("Sample", "Target", "Replicate"),
                                             value_col = "Cq",
                                             censored_col = "Cq_censored") {
    if (nrow(df) == 0) return(df)

    df |>
        dplyr::group_by(dplyr::across(dplyr::any_of(group_cols))) |>
        dplyr::filter({
            value <- .data[[value_col]]
            censored <- .data[[censored_col]]
            has_detected <- any(!is.na(value) & !censored)
            if (has_detected) {
                !is.na(value) & !censored
            } else {
                !is.na(value) & censored
            }
        }) |>
        dplyr::ungroup()
}

# Return the names of groups whose observations are all censored.
#
# `group_col` provides group labels and `censored_col` provides logical flags.
# Missing flags are not considered censored. Empty or incompatible data frames
# produce an empty character vector.
all_censored_groups <- function(df, group_col = "Sample", censored_col = "Cq_censored") {
    if (nrow(df) == 0 || !all(c(group_col, censored_col) %in% names(df))) {
        return(character())
    }

    groups <- unique(as.character(df[[group_col]]))
    groups[vapply(groups, function(group) {
        values <- df[[censored_col]][as.character(df[[group_col]]) == group]
        length(values) > 0 && all(!is.na(values) & values)
    }, logical(1))]
}

# Format qPCR numbers without vector-wide whitespace or zero padding.
#
# Values are rounded to `digits` significant figures. Missing
# values remain NA_character_. Each element is formatted independently.
format_qpcr_number <- function(x, digits = 4) {
    # format() pads every element to the greatest precision present in a
    # vector (e.g. c(21.7, 21.72) becomes c("21.70", "21.72")). Format each
    # value independently so data entry does not create cosmetic conversions.
    unname(vapply(x, function(value) {
        if (is.na(value)) return(NA_character_)
        format(
            signif(value, digits),
            trim = TRUE,
            nsmall = 0
        )
    }, character(1)))
}

# Identify parsed Cq conversions that are useful to report to the user.
#
# `normalized_cq` contains the final Cq text that will be displayed in the raw
# data table after parsing and censoring normalization.
#
# The result is TRUE for changed, non-empty inputs involving a decimal comma,
# explicit censoring, or text that could not be parsed. Purely cosmetic numeric
# formatting differences are not reported.
meaningful_cq_conversion <- function(original, parsed_numeric, censored, normalized_cq) {
    original <- trimws(as.character(original))
    censored <- !is.na(censored) & censored
    changed <- !is.na(original) & original != normalized_cq

    decimal_comma_changed <- !censored & is.finite(parsed_numeric) &
        grepl(",", original, fixed = TRUE)
    censored_changed <- censored
    invalid_changed <- !censored & is.na(parsed_numeric) & nzchar(original)

    changed & nzchar(original) &
        (decimal_comma_changed | censored_changed | invalid_changed)
}

# Format numeric values together with their censoring direction.
#
# Right-censored values receive a ">" prefix and left-censored values a "<"
# prefix. Uncensored values have no prefix, and missing values remain missing.
format_censored_value <- function(value, censored = FALSE,
                                  direction = c("right", "left"), digits = 4) {
    direction <- match.arg(direction)
    value_label <- format_qpcr_number(value, digits)
    censored <- !is.na(censored) & censored
    prefix <- if (direction == "right") ">" else "<"

    ifelse(
        is.na(value),
        NA_character_,
        ifelse(censored, paste0(prefix, value_label), value_label)
    )
}

# Map a result metric to the column containing its censoring flag.
#
# Exponentiated metrics reuse the flag of their untransformed counterpart.
# Housekeeping metrics cannot be censored and therefore return NULL.
censoring_column_for <- function(metric) {
    metric <- sub("^exp_", "", metric)
    if (metric == "HK_mean_Cq" || grepl("^HK_mean_.+_Cq$", metric)) {
        return(NULL)
    }
    if (metric == "ref_mean_dCq") return("ref_dCq_censored")

    paste0(sub("_mean$", "", metric), "_censored")
}

# Determine the display direction for a censored result metric.
#
# Cq-scale metrics are right-censored and exponentiated metrics are
# left-censored. A negative transformation `sign` reverses that direction.
censoring_direction_for <- function(metric, sign = 1) {
    direction <- if (grepl("^exp_", metric)) "left" else "right"
    if (sign < 0) {
        direction <- if (direction == "right") "left" else "right"
    }
    direction
}
