# =============================================================================
# qPCR Analysis Shiny App
# =============================================================================

library(shiny)
library(bslib)
library(tidyverse)
library(rhandsontable)

# =============================================================================
# Example Data
# =============================================================================

example_data <- data.frame(
    Sample = rep(c("Sample1", "Sample2", "Sample3"), each = 4 * 3),
    Target = rep(rep(c("HK1", "HK2", "TG1", "TG2"), each = 3), 3),
    Cq     = round(rnorm(3 * 3 * 4, mean = 23.5, sd = 1.2), 2)
)

# =============================================================================
# UI Definition
# =============================================================================

ui <- page_fillable(
    # Custom CSS
    tags$head(
        tags$style(HTML("
            /* Prevent horizontal scrollbars in rhandsontable */
            .handsontable {
                overflow-x: hidden;
                overflow-y: scroll;
            }
            
            /* Handsontable context menu above bslib fullscreen overlays */
            .htContextMenu:not(.htGhostTable) {
                z-index: 30000 !important;
            }
            
            /* Fix hidden row numbers in overflowing table */
            .handsontable.ht_clone_left {
                height: auto !important;
            }
                    
            .handsontable .wtHolder {
                height: 100% !important;
            }
            
            /* Style for excluded samples */
            .sample-excluded {
                opacity: 0.5;
                text-decoration: line-through;
            }
        "))
    ),
    
    # Main navigation
    navset_tab(
        id = "main_tabs",
        
        nav_item(h5("qPCR Analysis Tool")),
        nav_spacer(),
        
        # ---------------------------------------------------------------------
        # Panel 1: Input Data
        # ---------------------------------------------------------------------
        nav_panel(
            title = "Input Data",
            layout_sidebar(
                sidebar = sidebar(
                    title = "Sample Controls",
                    open = TRUE,
                    width = 380,
                    
                    checkboxInput(
                        "include_replicates",
                        "Include biological replicates column",
                        value = FALSE
                    ),
                    
                    hr(),
                    helpText("Drag rows to reorder. Edit 'New Label' to rename."),
                    helpText("Uncheck 'Include' to exclude samples from analysis."),
                    
                    rHandsontableOutput("samples_tab", height = "400px")
                ),
                
                # Main content area
                card(
                    max_height = "500px",
                    full_screen = TRUE,
                    fillable = TRUE,
                    card_header("qPCR Data Entry"),
                    rHandsontableOutput("raw_data", height = "90%"),
                    actionButton(
                        "load_example",
                        "Load Example Data",
                        width = "230px",
                        class = "btn-outline-secondary btn-sm"
                    )
                )
            )
        ),
        
        # ---------------------------------------------------------------------
        # Panel 2: Analysis
        # ---------------------------------------------------------------------
        nav_panel(
            title = "Analysis",
            layout_sidebar()
        ),
        
        # ---------------------------------------------------------------------
        # Panel 3: Results
        # ---------------------------------------------------------------------
        nav_panel(
            title = "Results",
            layout_sidebar()
        )
    )
)

# =============================================================================
# Server Logic
# =============================================================================

server <- function(input, output, session) {
    
    # -------------------------------------------------------------------------
    # Reactive Values
    # -------------------------------------------------------------------------
    
    values <- reactiveValues(
        # Raw Cq data as entered by user
        ct_raw_df = data.frame(
            Sample = character(5),
            Target = character(5),
            Cq     = numeric(5),
            stringsAsFactors = FALSE
        ),
        
        # Initialize empty sample control table (rename, reorder, exclude)
        # this has to be used only as cached last state 
        # samples rename, reordering and exclusion should be retrieved from input$samples_tab
        cached_samples_tab = data.frame(
            Sample    = character(),
            New_Label = character(),
            Include   = logical(),
            stringsAsFactors = FALSE
        )
    )
  
    # -------------------------------------------------------------------------
    # Derived Reactive: Processed data (with renames, ordering, exclusions)
    # -------------------------------------------------------------------------
    
    # ct_processed <- reactive({
    #     req(values$ct_raw_df)
    #     req(nrow(hot_to_r(input$samples_tab)) > 0)
    #     
    #     samples_df <- hot_to_r(input$samples_tab)
    #     
    #     values$ct_raw_df |>
    #         inner_join(samples_df) |>
    #         filter(Include) |>
    #         # update and reorder sample name
    #         mutate(Sample = factor(New_Label,
    #                                levels = samples_df$New_Label)) |>
    #         arrange(New_Label) |>
    #         select(-New_Label, -Include)
    # })
    
    # -------------------------------------------------------------------------
    # Observer: update cached_samples_tab when samples change in raw_data while preserving previous edits
    # -------------------------------------------------------------------------
    
    observeEvent(input$raw_data, {
        
        # cache last state of samples_tab
        values$cached_samples_tab <- hot_to_r(input$samples_tab)
        
        current_samples  <- hot_to_r(input$raw_data)$Sample |> unique()
        # drop empty sample names
        current_samples <- current_samples[!is.na(current_samples) & nzchar(trimws(current_samples))] 
        
        previous_samples <- values$cached_samples_tab$Sample
        new_samples <- setdiff(current_samples, previous_samples)
        
        if (length(current_samples) == 0) {
            # No valid samples, keep empty
            values$cached_samples_tab <- data.frame(
                Sample    = character(),
                New_Label = character(),
                Include   = logical()
            )
        } else if (length(previous_samples) == 0) {
            # First time cached_samples_tab update
            values$cached_samples_tab <- data.frame(
                Sample    = current_samples,
                New_Label = current_samples,
                Include   = rep(TRUE, times=length(current_samples))
            )
         } else {
            values$cached_samples_tab <- data.frame(Sample = current_samples) |>
                # reapply previous edits for existing samples
                left_join(values$cached_samples_tab) |>
                # fill in defaults for new samples
                mutate(
                    New_Label = ifelse(is.na(New_Label), Sample, New_Label),
                    Include   = ifelse(is.na(Include), TRUE, Include)
                )
         }
    })
    
    # -------------------------------------------------------------------------
    # Observer: Toggle biological replicates column
    # -------------------------------------------------------------------------
    
    observeEvent(input$include_replicates, {
        if (is.null(input$raw_data)) return()
        current_data <- hot_to_r(input$raw_data)
        
        if (input$include_replicates) {
            if (!"Replicate" %in% names(current_data)) {
                values$ct_raw_df <- cbind(
                    current_data,
                    Replicate = rep("R1", nrow(current_data))
                )
            }
        } else {
            if ("Replicate" %in% names(current_data)) {
                values$ct_raw_df <- current_data[
                    , !names(current_data) %in% "Replicate", drop = FALSE
                ]
            }
        }
    })
    
    # -------------------------------------------------------------------------
    # Observer: Load example data
    # -------------------------------------------------------------------------
    
    observeEvent(input$load_example, {
        data_to_load <- example_data
        
        if (input$include_replicates) {
            data_to_load$Replicate <- rep("R1", nrow(data_to_load))
        }
        
        values$ct_raw_df <- data_to_load
    })
    
    # -------------------------------------------------------------------------
    # Observer: Update raw data from user edits
    # -------------------------------------------------------------------------
    
    observeEvent(input$raw_data, {
        if (!is.null(input$raw_data)) {
            values$ct_raw_df <- hot_to_r(input$raw_data)
        }
    })
    
    # -------------------------------------------------------------------------
    # Output: Raw data table
    # -------------------------------------------------------------------------
    
    output$raw_data <- renderRHandsontable({
        req(values$ct_raw_df)
        
        rhandsontable(
            values$ct_raw_df,
            rowHeaders = TRUE,
            readOnly = FALSE,
            contextMenu = TRUE,
            stretchH = "all"
        ) %>%
            hot_col("Sample", type = "text") %>%
            hot_col("Target", type = "text") %>%
            hot_col("Cq", type = "numeric", format = "0.00") %>%
            hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)

    })
    
    # -------------------------------------------------------------------------
    # Output: Sample control table
    # -------------------------------------------------------------------------
    output$samples_tab <- renderRHandsontable({
        req(nrow(values$cached_samples_tab) > 0)
        
        rhandsontable(
            values$cached_samples_tab,
            rowHeaders = TRUE,
            readOnly = FALSE,
            stretchH = "all",
            contextMenu = FALSE,
            manualRowMove = TRUE
        ) |>
            hot_col("Sample", readOnly = TRUE) %>%
            hot_col("New_Label", type = "text") %>%
            hot_col("Include", type = "checkbox") %>%
            hot_cols(columnSorting = FALSE)
            
    })

}

# =============================================================================
# Run App
# =============================================================================

shinyApp(ui = ui, server = server)
