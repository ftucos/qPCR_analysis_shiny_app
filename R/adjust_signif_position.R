# Functions to add and adjust significance bar positions for ggplot2

# Add xmin, xmax, y.position columns based on sample factor levels
# Works with any statistical result format (not just rstatix)
add_signif_xy_position <- function(res, samples, y_max, step = 0.5) {
    # samples should be a factor with levels in plot order
    sample_positions <- setNames(seq_along(levels(samples)), levels(samples))
    
    # Convert group1 and group2 to numeric positions
    res <- res |>
        mutate(
            xmin = sample_positions[as.character(group1)],
            xmax = sample_positions[as.character(group2)],
            y.position = y_max + step  # Initial y position
        )
    
    # Make sure xmin < xmax (swap if needed)
    res <- res |>
        mutate(
            .xmin_temp = pmin(xmin, xmax),
            .xmax_temp = pmax(xmin, xmax),
            xmin = .xmin_temp,
            xmax = .xmax_temp
        ) |>
        select(-.xmin_temp, -.xmax_temp)
    
    return(res)
}

# Adjust y-positions of significance bars to minimize height
# Non-overlapping bars (based on xmin/xmax) are placed at the same y-level.
# Shorter bars are placed at the bottom, longer bars at the top.
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

# Prepare significance data for plotting
prepare_signif_data <- function(stats_result, samples, y_max, 
                                 hide_ns = FALSE, show_p_value = FALSE,
                                 step = NULL) {
    
    # Get the comparison results
    res <- stats_result$test_res |>
        # remove rows where group1 or group2 are not defined
        # this happens for test like 2 sample ANCOVA that include also covariate a and residuals stats
        drop_na(group1, group2)

    if (is.null(res)) return(NULL)
    
    # Ensure we have required columns
    if (!all(c("group1", "group2", "Significance") %in% names(res))) {
        return(NULL)
    }
    
    # Filter out ns if requested
    if (hide_ns) {
        res <- res |> filter(Significance != "ns")
    }
    
    if (nrow(res) == 0) return(NULL)
    
    # Create label column
    if (show_p_value) {
        # Find the p-value column (could be p-value, Adj. p-value, etc.)
        p_col <- intersect(c("Adj. p-value", "p-value"), names(res))[1]
        if (!is.na(p_col)) {
            res <- res |>
                mutate(label = ifelse(
                    .data[[p_col]] < 0.001,
                    "p < 0.001",
                    paste0("p = ", format(round(.data[[p_col]], 3), nsmall = 3))
                ))
        } else {
            res <- res |> mutate(label = Significance)
        }
    } else {
        res <- res |> mutate(label = Significance)
    }
    
    # Calculate step if not provided
    if (is.null(step)) {
        y_range <- y_max - min(0, y_max)
        step <- y_range * 0.08  # 8% of the y range
    }
    
    # Add x/y positions
    res <- res |>
        add_signif_xy_position(samples, y_max, step) |>
        adjust_y_position(step = step)
    
    return(res)
}
