# Build ggplot objects for the Results and Export tabs
# Both functions consume the output of result_plot_data() plus tab-specific args

# Results tab plot (interactive, will be converted to plotly) -----------------

build_results_plot <- function(plot_data, accent_color, secondary_color) {
    
    df_target         <- plot_data$df_target
    df_summary_target <- plot_data$df_summary
    y_value           <- plot_data$y_value
    y_label           <- plot_data$y_label
    sign              <- plot_data$sign
    y_summary_value   <- plot_data$y_summary_value
    error_bar_high    <- plot_data$error_bar_high
    error_bar_low     <- plot_data$error_bar_low
    y_limits          <- plot_data$y_limits
    y_min_label       <- plot_data$y_min_label
    error_bar_label   <- plot_data$error_bar_label
    out_metric        <- plot_data$out_metric
    stat_type         <- plot_data$stat_type
    summarize_bio_reps <- plot_data$summarize_bio_reps
    
    # Plotly hover text labels
    has_replicate <- "Replicate" %in% names(df_target)
    
    df_target <- df_target |>
        # build label to show on hover
        mutate(text = if (has_replicate) {
            glue(
                "{Sample} ({Replicate})
                Target: {Target}
                {y_label}: {round(sign*.data[[y_value]], 2)}"
            )
        } else {
            glue(
                "{Sample}
                Target: {Target}
                {y_label}: {round(sign*.data[[y_value]], 2)}"
            )
        })
    
    df_summary_target <- df_summary_target |>
        mutate(text = glue(
            "{Sample}
            Target: {Target}
            Mean {y_label}: {round(sign*.data[[y_summary_value]], 2)}
            {error_bar_label}"
        ))
    
    # Reference line
    if (out_metric == "exp_ddCq") {
        p <- ggplot() +
            geom_hline(yintercept = 1, linetype = "dashed", color = "gray30")
    } else if (out_metric == "ddCq") {
        p <- ggplot() +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray30")
    } else {
        p <- ggplot()
    }
    
    # Average bars (not for dCq — zero has no meaning)
    if (out_metric != "dCq") {
        p <- p +
            geom_bar(
                data = df_summary_target,
                aes(x = Sample, y = sign * .data[[y_summary_value]], text = text),
                stat = "identity",
                fill = accent_color,
                alpha = 0.5,
                width = 0.6
            )
    }
    
    # Beeswarm points
    p <- p +
        geom_beeswarm(
            data = df_target,
            aes(
                x = Sample, y = sign * .data[[y_value]],
                text = text,
                color = point_type_label,
                shape = point_type_label
            ),
            method = "compactswarm", preserve.data.axis = TRUE
        )
    
    # Mean points for dCq
    if (out_metric == "dCq") {
        p <- p +
            geom_point(
                data = df_summary_target |>
                    mutate(point_type_label = "Mean"),
                aes(
                    x = Sample, y = sign * .data[[y_summary_value]],
                    text = text,
                    shape = point_type_label,
                    color = point_type_label
                ),
                inherit.aes = F,
                size = 4
            )
    }
    
    p <- p +
        scale_shape_manual(
            values = c("Detected" = 16, "Undetected" = 1, "Mean" = 4),
            name = "",
        ) +
        scale_color_manual(
            values = c("Detected" = secondary_color, "Undetected" = "#C03A2B", "Mean" = accent_color),
            name = "",
        ) +
        labs(x = "Sample", y = y_label, title = NULL) +
        scale_y_continuous(
            expand = expansion(mult = 0.05, add = 0),
            labels = function(x) ifelse(x == y_limits[1], y_min_label, x),
            oob = function(x, range) squish_infinite_to_val(x, range, to_value = abs(y_limits[1]))
        ) +
        coord_cartesian(ylim = y_limits) +
        theme_minimal(base_size = 14) +
        theme(
            legend.position = "bottom",
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    
    # Add error bars if requested
    if (stat_type != "none") {
        p <- p +
            geom_errorbar(
                data = df_summary_target,
                aes(
                    x    = Sample,
                    ymin = sign * .data[[error_bar_low]],
                    ymax = sign * .data[[error_bar_high]]
                ),
                inherit.aes = F,
                width = 0.4,
                color = accent_color
            )
    }
    
    # Facet by replicate if present
    if (summarize_bio_reps == "split" & ("Replicate" %in% names(df_target))) {
        p <- p + ggh4x::facet_wrap2(~Replicate, axes = "all")
    }
    
    return(p)
}


# Export tab plot (static, print-quality) ------------------------------------

build_export_plot <- function(plot_data, colors, lw, axis_text_size, signif_text_size,
                              bar_width, plot_width, plot_height,
                              show_signif_bars, stats_result,
                              hide_ns, show_exact_pvalue) {
    
    df_target         <- plot_data$df_target
    df_summary_target <- plot_data$df_summary
    y_value           <- plot_data$y_value
    y_label           <- plot_data$y_label
    sign              <- plot_data$sign
    y_summary_value   <- plot_data$y_summary_value
    error_bar_high    <- plot_data$error_bar_high
    error_bar_low     <- plot_data$error_bar_low
    y_limits          <- plot_data$y_limits
    y_min_label       <- plot_data$y_min_label
    out_metric        <- plot_data$out_metric
    stat_type         <- plot_data$stat_type
    summarize_bio_reps <- plot_data$summarize_bio_reps
    
    # For exp_dCq and exp_ddCq, force lower limit to 0
    if (out_metric %in% c("exp_dCq", "exp_ddCq")) {
        y_limits[1] <- 0
    }
    
    # Build title and subtitle
    plot_title <- plot_data$target_name
    plot_subtitle <- NULL
    if (isTRUE(show_signif_bars) && !is.null(stats_result) && !is.null(stats_result$omnibus_label)) {
        plot_subtitle <- stats_result$omnibus_label
    }
    
    # Reference line
    if (out_metric == "exp_ddCq") {
        p <- ggplot() +
            geom_hline(yintercept = 1, linetype = "dashed", color = "gray30", linewidth = lw)
    } else if (out_metric == "ddCq") {
        p <- ggplot() +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = lw)
    } else {
        p <- ggplot()
    }
    
    # Bars (not for dCq)
    if (out_metric != "dCq") {
        p <- p +
            geom_bar(
                data = df_summary_target,
                aes(x = Sample, y = sign * .data[[y_summary_value]], fill = Sample),
                stat = "identity",
                alpha = 0.7,
                width = bar_width,
                color = "black",
                linewidth = lw
            ) +
            scale_fill_manual(values = colors, guide = "none")
    }
    
    # Beeswarm points
    p <- p +
        geom_beeswarm(
            data = df_target,
            aes(
                x = Sample, y = sign * .data[[y_value]],
                shape = point_type_label
            ),
            method = "compactswarm", preserve.data.axis = TRUE,
            color = "black"
        )
    
    # Mean points for dCq
    if (out_metric == "dCq") {
        p <- p +
            geom_point(
                data = df_summary_target |> mutate(point_type_label = "Mean"),
                aes(x = Sample, y = sign * .data[[y_summary_value]], shape = point_type_label),
                color = "black",
                size = 4
            )
    }
    
    p <- p +
        scale_shape_manual(
            values = c("Detected" = 16, "Undetected" = 1, "Mean" = 4),
            name = ""
        ) +
        labs(x = NULL, y = y_label, title = plot_title, subtitle = plot_subtitle) +
        scale_y_continuous(
            expand = if (out_metric %in% c("exp_dCq", "exp_ddCq")) expansion(mult = c(0, 0.05), add = 0) else expansion(mult = 0.05, add = 0),
            labels = function(x) ifelse(x == y_limits[1], y_min_label, x),
            oob = function(x, range) squish_infinite_to_val(x, range, to_value = abs(y_limits[1]))
        ) +
        coord_cartesian(ylim = y_limits, clip = "off") +
        theme_minimal(base_size = axis_text_size) +
        theme(
            legend.position = "bottom",
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
            axis.text.y = element_text(size = axis_text_size),
            axis.line = element_line(linewidth = lw, color = "black"),
            axis.ticks = element_line(linewidth = lw, color = "black")
        )
    
    # Error bars
    if (stat_type != "none") {
        p <- p +
            geom_errorbar(
                data = df_summary_target,
                aes(
                    x = Sample,
                    ymin = sign * .data[[error_bar_low]],
                    ymax = sign * .data[[error_bar_high]]
                ),
                width = 0.4,
                color = "black",
                linewidth = lw
            )
    }
    
    # Significance bars
    if (isTRUE(show_signif_bars) && !is.null(stats_result) && is.null(stats_result$error)) {
        signif_data <- prepare_signif_data(
            stats_result,
            samples = df_target$Sample,
            # TODO: adjust max based on errorbars (or datapoints?) and define step size
            y_max = sign * max(df_summary_target[[y_summary_value]], na.rm = TRUE),
            hide_ns = isTRUE(hide_ns),
            show_p_value = isTRUE(show_exact_pvalue)
        )
        
        if (!is.null(signif_data) && nrow(signif_data) > 0) {
            signif_data <- signif_data |>
                mutate(y.position = sign * abs(y.position))
            
            p <- p +
                ggsignif::geom_signif(
                    data = signif_data,
                    aes(
                        xmin = group1,
                        xmax = group2,
                        y_position = y.position,
                        annotations = label
                    ),
                    manual = TRUE,
                    textsize = signif_text_size / ggplot2::.pt,
                    vjust = 0.5,
                    tip_length = 0, extend_line = -0.01,
                    margin_top = 0.1,
                    size = lw
                )
        }
    }
    
    # Facet by replicate if present
    if (summarize_bio_reps == "split" & ("Replicate" %in% names(df_target))) {
        p <- p + ggh4x::facet_wrap2(~Replicate, axes = "all")
    }
    
    # Force panel size
    p <- p +
        ggh4x::force_panelsizes(
            rows = unit(plot_height, "in"),
            cols = unit(plot_width, "in")
        )
    
    return(p)
}
