
# use it for averaging Cq and dCq
mean_handle_inf <- function(x) {
    # drop NA
    x <- x[!is.na(x)]
    if (length(x) == 0) {
        return(NA)
    }
    
    # drop Inf
    x <- x[is.finite(x)]
    if (length(x) == 0) {
        return(Inf)
    } else {
        return(mean(x, na.rm = TRUE))
    }
    
}

sd_handle_inf <- function(x) {
    # drop NA
    x <- x[!is.na(x)]
    # drop Inf
    x <- x[is.finite(x)]
    sd(x)
}

# use it for averaging exponentiated dCq and ddCq values
mean_handle_0 <- function(x) {
    # drop NA
    x <- x[!is.na(x)]
    if (length(x) == 0) {
        return(NA)
    }
    
    # drop Inf
    x <- x[x != 0]
    if (length(x) == 0) {
        return(0)
    } else {
        return(mean(x, na.rm = TRUE))
    }
}

n_valid_Cq <- function(x) {
    # drop NA
    x <- x[!is.na(x)]
    
    # drop Inf
    x <- x[is.finite(x)]
    return(length(x))
}
