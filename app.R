library(shiny)
library(bslib)
library(rhandsontable)

ui <- page_fillable(
    # custom CSS
    tags$head(
        # Note the wrapping of the string in HTML()
        tags$style(HTML("
        /* prevent horizontal scrollbars in rhandsontable */
        .handsontable {
            overflow-x: hidden;
            overflow-y: scroll
        }
        
        /* Handsontable context menu & dropdowns above bslib fullscreen overlays */
        .htContextMenu:not(.htGhostTable) {
          z-index: 30000 !important;
        }
        
        /* fix hidden row number in the overfwloing part of the table */
        .handsontable.ht_clone_left {
            height: auto !important;
        }
                
        .handsontable .wtHolder {
             height: 100% !important;
        }
        
        "))
    ),
    
    navset_tab(
        id = "main_tabs",
        #title = "qPCR Analysis Tool",
        nav_item(
            h5("qPCR Analysis Tool")
        ),
        nav_spacer(),
        # Panel 1: Input Data
        nav_panel(
            title = "Input Data",
            layout_sidebar(
                sidebar = sidebar(
                    title = "Data Options",
                    open = TRUE,
                    checkboxInput(
                        "include_replicates",
                        "Include biological replicates column",
                        value = FALSE
                    ),
                    hr(),
                    helpText("Paste your Raw Cq data directly into the table below."),

                ),
                
                # Main content area
                card(max_height = "600px",
                    full_screen = TRUE,
                    fillable = T,
                    card_header("qPCR Data Entry"),
                    rHandsontableOutput("raw_data", height = "90%"),
                    actionButton(
                        "load_example",
                        "Load Example Data",
                        width = "230px",
                        class = "btn-secondary"
                    )
                )
            )
        ),
        
        # Panel 2: Analysis
        nav_panel(
            title = "Analysis",
            layout_sidebar(
            )
        ),
        
        # Panel 3: Results
        nav_panel(
            title = "Results",
            layout_sidebar(
            )
        )
    )
)

server <- function(input, output, session) {
    
    # Reactive values to store data
    values <- reactiveValues(
        ct_data = data.frame(
            Sample = character(5),
            Target = character(5),
            Cq = numeric(5),
            stringsAsFactors = FALSE
        )
    )
    
    # Update data structure when replicate checkbox changes
    observeEvent(input$include_replicates, {
        current_data <- hot_to_r(input$raw_data)
        
        if (input$include_replicates) {
            if (!"Replicate" %in% names(current_data)) {
                values$ct_data <- cbind(current_data, Replicate = character(nrow(current_data)))
            }
        } else {
            if ("Replicate" %in% names(current_data)) {
                values$ct_data <- current_data[, !names(current_data) %in% "Replicate"]
            }
        }
    })
    
    # Load example data
    observeEvent(input$load_example, {
        example_data <- data.frame(
            Sample = rep(c("Sample1", "Sample2", "Sample3"), each = 4 * 3),              # 4 genes * 3 tech reps
            Target      = rep(rep(c("HK1", "HK2", "TG1", "TG2"), each = 3), 3),
            Cq        = round(rnorm(3 * 3 * 4, mean = 23.5, sd = 1.2), 2)
        )
        
        if (input$include_replicates) {
            example_data$Replicate <- c("R1", "R2", "R3", "R1", "R2", "R3")
        }
        
        values$ct_data <- example_data
    })
    
    # Render the rhandsontable
    output$raw_data <- renderRHandsontable({
        df <- values$ct_data
        
        rhandsontable(df,
                      rowHeaders = TRUE,
                      readOnly = FALSE,
                      contextMenu = TRUE,
                      stretchH = "all") %>%
            hot_col("Sample", type = "text") %>%
            hot_col("Target", type = "text") %>%
            hot_col("Cq", type = "numeric", format = "0.00") %>%
            hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
    })
    
    # Update values when table is edited
    observeEvent(input$raw_data, {
        values$ct_data <- hot_to_r(input$raw_data)
    })
    

    
    
}

shinyApp(ui = ui, server = server)
