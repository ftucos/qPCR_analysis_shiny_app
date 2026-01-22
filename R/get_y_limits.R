get_y_limits <- function(x, min_range = 3, margin = c(0, 0), undetected_value = 40) {
    # remove NA
    x <- x[!is.na(x)]
    
    undetected_present <- any(!is.finite(x))
    
    # return fixed range if all values are undetected
    if (all(!is.finite(x))) {
        return(c(34, 40))
    }

    # remove infinite values
    x <- x[is.finite(x)]
    
    y_min <- min(x)
    y_max <- max(x)
    
    # if undetected valeus are present, expand the scale up to 40
    if(undetected_present) {
        y_max <- max(y_max, undetected_value)
    }
    
    y_range <- y_max - y_min
    
    if(y_range < min_range) {
        y_center <- (y_max + y_min) / 2
        y_min = y_center - min_range / 2
        y_max = y_center + min_range / 2
    }
    
    return(c(y_min - margin[1], y_max + tail(margin, 1)))
}

