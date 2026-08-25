# adapted from https://stackoverflow.com/a/69035732/7793290
fix_plotly_legend <- function(gp) {
    for (i in seq_along(gp$x$data)) {
        nm <- gp$x$data[[i]]$name
        # Extract the group identifier and assign it to the name and legendgroup arguments
        gp$x$data[[i]]$name <- str_remove(nm, "^\\(") |> str_remove(",?(NA)?\\)$")
        gp$x$data[[i]]$legendgroup <- gp$x$data[[i]]$name
    }
    
    gp
}

# Disable hover labels for empty filled regions created by ggplot annotations
# (for example, the gray borderline-Cq band). With faceting, ggplotly creates
# one such trace per panel, so identify them by their trace properties rather
# than by position.
suppress_empty_filled_hover <- function(gp) {
    for (i in seq_along(gp$x$data)) {
        trace <- gp$x$data[[i]]
        text <- trace$text
        has_no_text <- is.null(text) ||
            all(is.na(text) | !nzchar(as.character(text)))

        if (identical(trace$fill, "toself") &&
            identical(trace$mode, "lines") &&
            has_no_text) {
            gp$x$data[[i]]$hoverinfo <- "skip"
        }
    }

    gp
}
