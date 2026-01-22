squish_infinite_to_val <- function (x, range = c(0, 1), to_value = 40) 
    {
        force(range)
        x[x == -Inf] <- -to_value
        x[x == Inf] <- to_value
        x
    }
