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

