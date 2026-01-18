get_y_limits <- function(x, min_range = 3) {
    y_min <- min(x, na.rm = TRUE)
    y_max <- max(x, na.rm = TRUE)
    y_range <- y_max - y_min
    
    if(y_range < min_range) {
        y_center <- (y_max + y_min) / 2
        return(c(y_center - min_range / 2,
                 y_center + min_range / 2))
    } else {
        return(c(y_min, y_max))
    }
}