# Return a numeric mean or NA, never NaN, when no observed values are present.
mean_or_na <- function(x) {
    if (!any(!is.na(x))) {
        return(NA_real_)
    }

    mean(x, na.rm = TRUE)
}

count_non_missing <- function(x) {
    sum(!is.na(x))
}
