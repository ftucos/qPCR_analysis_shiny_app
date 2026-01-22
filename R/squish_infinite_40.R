squish_infinite_40 <- function (x, range = c(0, 1)) 
    {
        force(range)
        x[x == -Inf] <- -40
        x[x == Inf] <- 40
        x
    }
