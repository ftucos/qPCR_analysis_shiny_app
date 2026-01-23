get_Cq_y_limits <- function(x, min_range = 3, margin = c(0, 0), undetected_value = 40) {
    
    # remove NA
    x <- x[!is.na(x)]
    
    undetected_present <- any(!is.finite(x))
    
    # return fixed range if all values are undetected
    if (all(!is.finite(x))) {
        return(c(34, 40))
    }

    # remove infinite values
    x <- x[is.finite(x)]
    
    y_min <- floor(min(x))
    y_max <- ceiling(max(x))
    
    # if undetected valeus are present, expand the scale up to 40
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


# for dCq and exp_dCq plots ---------------------------------------------------
get_dCq_y_limits <- function(df_target, df_summary_target, stat_type = c("none"), metric = c("dCt")) {
    
    values = df_target[[metric]]
    
    if (stat_type != "none") {
        error_bar_high  <- glue("{metric}_{stat_type}_high") 
        error_bar_low   <- glue("{metric}_{stat_type}_low")
    }
    
    if(stat_type != "none") {
        values <- values |> 
            append(c(
                df_summary_target |> pull(error_bar_low),
                df_summary_target |> pull(error_bar_high)
            ))
    }
    
    # drop NA
    values <- values[!is.na(values)]
    
    if (metric == "dCq") {
        
        finite_values <- values[is.finite(values)]
        
        if(all(!is.finite(values))) {
            # all undetected
            y_min <- -40
            y_max <- -34
            
        } else if (any(!is.finite(values))) {
            # some undetected
            y_min <- floor(min(-finite_values) - 2)
            y_max <- max(-finite_values)
        } else {
            # no undetected
            y_min <- floor(min(-values))
            y_max <- ceiling(max(-values))
        }
        
    } else { # 2^-dCq
        # allways start at 0
        y_min <- 0
        
        print(values)
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

