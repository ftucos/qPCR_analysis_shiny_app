get_Cq_y_limits <- function(x, min_range = 3, margin = c(0, 0),
                            undetected_value = 40, undetected_present = FALSE) {
    
    # remove NA
    x <- x[!is.na(x)]
    
    # return fixed range if all values are undetected
    if (length(x) == 0) {
        return(c(undetected_value - 6, undetected_value))
    }

    if (undetected_present && all(dplyr::near(x, undetected_value))) {
        return(c(undetected_value - 6, undetected_value))
    }
    
    y_min <- floor(min(x))
    y_max <- ceiling(max(x))
    
    # if undetected values are present, expand the scale up to undetected_value
    if(undetected_present) {
        y_max <- max(ceiling(y_max), undetected_value)
    }
    
    y_range <- y_max - y_min
    
    if(y_range < min_range) {
        y_center <- (y_max + y_min) / 2
        y_min = y_center - min_range / 2
        y_max = y_center + min_range / 2
    }
    
    return(c(y_min - margin[1], y_max + tail(margin, 1)))
}


# for (d)dCq and exp_(d)dCq plots ---------------------------------------------------
get_y_limits <- function(values, metric = c("dCq"), undetected_value = 40) {
    
    # drop NA
    values <- values[!is.na(values)]
    
    if (metric %in% c("dCq", "ddCq")) {
        
        y_min <- floor(min(-values))
        y_max <- ceiling(max(-values))
        
    } else { # 2^-d(d)Cq
        # allways start at 0
        y_min <- 0
        
        if(all(values == 0)) {
            # all undetected
            y_max <- 4
        } else {
            # no undetected
            y_max <- max(values)
        }
    }
    
    return(c(y_min, y_max))
}
