example_data <- data.frame(
    Sample = rep(c("Sample1", "Sample2", "Sample3"), each = 4 * 1),
    Cq     = round(rnorm(3 * 1 * 4, mean = 23.5, sd = 1.2), 2),
    Keep   = sample(c(TRUE, FALSE), prob = c(0.9, 0.1), 12, replace = TRUE)
) |>
    mutate(Keep_label = ifelse(Keep, "Included", "Excluded"))

df_summary_target <- example_data |>
    group_by(Sample) |>
    summarise(
        mean = mean(as.numeric(Cq)),
        .groups = "drop"
    )

p <- ggplot(example_data,
            aes(x = Sample, y = Cq,
                alpha = Keep_label,
                text = paste0(
                  Sample, "\n",
                  ifelse(Keep, "", " (excluded)")
                )
            )) +
  geom_quasirandom() +
  geom_point(
    data = df_summary_target,
    aes(x = Sample, y = mean, shape = "Mean"),   # <- map a label
    inherit.aes = FALSE,
    size = 4, color = "#007bc2"
  ) +
  scale_shape_manual(values = c("Mean" = 4), name = "") +  # <- show "Mean" as X
  labs(x = "Sample", y = "Cq", title = NULL) +
  scale_y_continuous(expand = c(0.1)) +
  scale_alpha_manual(
    values = c("Included" = 1, "Excluded" = 0.3),
    name = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

b <- plotly::plotly_build(p)

for (i in seq_along(b$x$data)) {
    nm <- b$x$data[[i]]$name
    if (!is.null(nm)) {
        b$x$data[[i]]$name <- sub("^\\(([^,]+),.*\\)$", "\\1", nm)
        b$x$data[[i]]$legendgroup <- b$x$data[[i]]$name
    }
}

# 2) Hide duplicate legend entries (plotly often creates multiple traces per group)
seen <- character()
for (i in seq_along(b$x$data)) {
    nm <- b$x$data[[i]]$name
    if (!is.null(nm) && nzchar(nm)) {
        if (nm %in% seen) b$x$data[[i]]$showlegend <- FALSE
        else seen <- c(seen, nm)
    }
}

b




