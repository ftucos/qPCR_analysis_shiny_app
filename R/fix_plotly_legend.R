# adapted from https://stackoverflow.com/a/59723151/7793290
fix_plotly_legend <- function(pplot) {
    for (i in seq_along(pplot$x$data)) {
        nm <- pplot$x$data[[i]]$name
        if (!is.null(nm)) {
            pplot$x$data[[i]]$name <- sub("^\\(([^,]+),.*\\)$", "\\1", nm)
        }
    }
    
    pplot
}