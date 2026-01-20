# adapted from https://stackoverflow.com/a/69035732/7793290
fix_plotly_legend <- function(gp) {
    for (i in seq_along(gp$x$data)) {
        # Is the layer the first entry of the group?
        is_first <- grepl("^\\(.*?,1\\)", gp$x$data[[i]]$name)
        
        nm <- gp$x$data[[i]]$name
        # Extract the group identifier and assign it to the name and legendgroup arguments
        gp$x$data[[i]]$name <- gsub("^\\((.*?),\\d+\\)", "\\1", nm)
        gp$x$data[[i]]$legendgroup <- gp$x$data[[i]]$name
        # Show the legend only for the first layer of the group 
        if (!is_first) gp$x$data[[i]]$showlegend <- FALSE
    }
    
    gp
}

