library(tidyverse)
library(ggpubr)

df <- data.frame(
    Sample = rep(c("A", "B", "C", "D", "E"), each = 10),
    dCq = c(
        rnorm(10, mean = 5, sd = 1),
        rnorm(10, mean = 7, sd = 1),
        rnorm(10, mean = 6, sd = 1),
        rnorm(10, mean = 5.2, sd = 1),
        rnorm(10, mean = 9, sd = 1)
    )
)

# Adjust y-positions of significance bars to minimize height
adjust_y_position <- function(res, y_base = NULL, step = NULL) {
    # Calculate bar width and sort by width (shorter bars first, then by xmin)
    res <- res |> 
        mutate(.row_id = row_number(),
               .width = xmax - xmin) |>
        arrange(.width, xmin)
    
    # Calculate step from existing y.positions if not provided
    if (is.null(step)) {
        unique_y <- sort(unique(res$y.position))
        if (length(unique_y) > 1) {
            step <- min(diff(unique_y))
        } else {
            step <- 0.5  # fallback
        }
    }
    
    # Use minimum y.position as base if not provided
    if (is.null(y_base)) {
        y_base <- min(res$y.position)
    }
    
    n <- nrow(res)
    levels <- integer(n)  # Track which level each bar is assigned to
    level_ends <- numeric()  # Track the xmax of the rightmost bar at each level
    
    for (i in seq_len(n)) {
        xmin_i <- res$xmin[i]
        
        # Find the lowest level where this bar doesn't overlap
        assigned_level <- NA
        for (lvl in seq_along(level_ends)) {
            if (level_ends[lvl] <= xmin_i) {
                # No overlap - can use this level (allows shared endpoints)
                assigned_level <- lvl
                level_ends[lvl] <- res$xmax[i]
                break
            }
        }
        
        # If no existing level works, create a new one
        if (is.na(assigned_level)) {
            level_ends <- c(level_ends, res$xmax[i])
            assigned_level <- length(level_ends)
        }
        
        levels[i] <- assigned_level
    }
    
    # Assign y-positions based on levels
    res$y.position <- y_base + (levels - 1) * step
    
    # Restore original order and remove temp columns
    res <- res |> 
        arrange(.row_id) |>
        select(-.row_id, -.width)
    
    return(res)
}

res <- df |>
    rstatix::pairwise_t_test(
        dCq ~ Sample,
        p.adjust.method = "none"
        ) %>%
    rstatix::add_xy_position(x = "Sample")

res.1 <- adjust_y_position(res)

ggplot(df, aes(x = Sample, y=  dCq)) +
    geom_jitter(width = 0.2) +
    stat_pvalue_manual(data = res.1,bracket.shorten = 0.1, tip.length = 0)
    
