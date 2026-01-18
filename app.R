# =============================================================================
# qPCR Analysis Shiny App
# =============================================================================

library(shiny)
library(bslib)
library(shinyWidgets)
library(tidyverse)
library(rhandsontable)
library(plotly)
library(ggbeeswarm)
library(scales)

# =============================================================================
# Example Data
# =============================================================================

# empty raw data table to initialize rhandsontable
empty_raw_data <- data.frame(
    Sample = character(5),
    Target = character(5),
    Cq     = character(5),
    stringsAsFactors = FALSE
)

example_data <- data.frame(
    Sample = rep(c("Sample1", "Sample2", "Sample3"), each = 4 * 3),
    Target = rep(rep(c("HK1", "HK2", "TG1", "TG2"), each = 3), 3),
    Cq     = round(rnorm(3 * 3 * 4, mean = 23.5, sd = 1.2), 2) |> as.character()
)


# helper functions ==========================================
drop_empty <- function(x) {x[!is.na(x) & nzchar(trimws(x))]}

source("R/parse_Cq.R")
source("R/get_y_limits.R")
source("R/fix_plotly_legend.R")

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
                fillable = T,
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
                    helpText("Edit 'New Label' to rename samples. Drag the row number to reorder them. Uncheck 'Include' to exclude samples from analysis."),
                    rHandsontableOutput("samples_tab")
                ),
                
                # Main content area
                card(
                    max_height = "500px",
                    full_screen = TRUE,
                    fillable = TRUE,
                    card_header("qPCR Data Entry"),
                    helpText("Paste your qPCR data below directly from Excel"),
                    rHandsontableOutput("raw_data", height = "90%"),
                    # place two buttons next to each other
                    div(
                        style = "display: flex; gap: 10px",
                        class = "d-flex justify-content-center",
                        actionButton(
                            "load_example",
                            "Load Example Data",
                            width = "150px",
                            class = "btn-outline-secondary btn-sm"
                        ),
                        actionButton(
                            "clear_data",
                            "Clear Data",
                            width = "150px",
                            class = "btn-outline-danger btn-sm"
                        )
                    )
                )
            )
        ),
        
        # ---------------------------------------------------------------------
        # Panel 2: Cq values inspection and exclusion
        # ---------------------------------------------------------------------
        nav_panel(
            title = "Cq",
            fillable = T,
            layout_sidebar(
                sidebar = sidebar(
                    title = "Cq Inspection",
                    open = TRUE,
                    width = 380,
                    
                    selectInput(
                        "target_select",
                        "Select Target",
                        choices = NULL # popolate dinamically
                    ),
                    radioGroupButtons(
                        inputId = "aggr_function",
                        label = "Aggregate Technical Replicates by:",
                        choices = c("Mean" = "mean", "Median" = "median"),
                        justified = TRUE,   # fill the row
                        individual = FALSE,  # more "segmented" feel
                        width = "100%",
                        size = "sm"
                    ),
                    hr(),
                    br(),
                    helpText("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua."),
                    br()
                ),
                # Main content area
                card(
                    full_screen = TRUE,
                    fillable = TRUE,
                    card_header(textOutput("ct_plot_title", inline = TRUE)),
                    plotlyOutput("ct_plot", height = "100%")
                )
                #tableOutput("table_log")
            )
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
    # Current theme accent color
    # -------------------------------------------------------------------------
    accent_color <- reactive({
        bslib::bs_get_variables(bslib::bs_current_theme(session), "primary")[[1]]
    })
    
    secondary_color <- reactive({
        bslib::bs_get_variables(bslib::bs_current_theme(session), "secondary")[[1]]
    })
    
    # -------------------------------------------------------------------------
    # Reactive Values
    # -------------------------------------------------------------------------
    
    values <- reactiveValues(
        # cached raw data, required to add/remove replicate column
        cached_raw_data = empty_raw_data,
        
        # Initialize empty sample control table (rename, reorder, exclude)
        # this has to be used only as cached last state 
        # samples rename, reordering and exclusion should be retrieved from input$samples_tab
        cached_samples_tab = data.frame(
            Sample    = character(),
            New_Label = character(),
            Include   = logical(),
            stringsAsFactors = FALSE
        ),
        
        # list of excluded points
        excluded_point_keys = c()
    )

    
    # -------------------------------------------------------------------------
    # Observer: Toggle biological replicates column
    # -------------------------------------------------------------------------
    
    observeEvent(input$include_replicates, {
        if (is.null(input$raw_data)) return()
        
        current_data <- hot_to_r(input$raw_data)
        
        if (input$include_replicates) {
            # add replicate column if missing
            if (!"Replicate" %in% names(current_data)) {
                values$cached_raw_data <- current_data |>
                    mutate(Replicate = "R1")
            }
        } else {
            # drop replicate column if present
            values$cached_raw_data <- select(current_data, -any_of("Replicate"))
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
        
        values$cached_raw_data <- data_to_load
    })
    
    # -------------------------------------------------------------------------
    # Observer: clean data
    # -------------------------------------------------------------------------
    
    observeEvent(input$clear_data, {
        values$cached_raw_data <- empty_raw_data
        # reset the lsit of excluded points
        values$excluded_point_keys <- c()
    })

    # -------------------------------------------------------------------------
    # Output: Raw data table
    # -------------------------------------------------------------------------
    
    output$raw_data <- renderRHandsontable({
        req(values$cached_raw_data)

        rhandsontable(
            values$cached_raw_data,
            rowHeaders = TRUE,
            readOnly = FALSE,
            contextMenu = TRUE,
            stretchH = "all",
            renderAllRows = TRUE
        ) |>
            hot_col("Sample", type = "text") |>
            hot_col("Target", type = "text") |>
            hot_col("Cq", type = "text") |> # type = "text" to allow for Undetermined or other labels
            hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
        
    })
    
    # -------------------------------------------------------------------------
    # Observer: on raw data edit:
    # 1. cache last `raw_data` and validate Cq values
    # 2. update `cached_samples_tab` when samples change in raw_data while preserving previous edits 
    # -------------------------------------------------------------------------
    
    observeEvent(input$raw_data, {

        # 1. validate Cq values --------------------
        current_data <- hot_to_r(input$raw_data)
        
        original_cq <- current_data$Cq |> as.character() |> trimws() |> replace_na("")
        parsed_cq <- parse_Cq(original_cq) |> replace_na("")
        
        # Find values that changed
        changed_mask <- original_cq != parsed_cq
        conversions <- map2_chr(original_cq[changed_mask], parsed_cq[changed_mask],
                                ~ paste0(.x, " → ", .y)) |>
                                unique()
        
        # Update parsed Cq values
        current_data$Cq <- parsed_cq
        
        # Cache updated raw data
        values$cached_raw_data <- current_data
        
        # Show warning modal if any conversions happened
        if (length(conversions) > 0) {
            showModal(modalDialog(
                title = "Cq Values Converted",
                tags$p("The following conversions were applied:"),
                tags$ul(
                    lapply(conversions, function(x) {tags$li(x)})
                ),
                easyClose = TRUE,
                footer = modalButton("OK")
            ))
        }

        # 2. cache last state of samples_tab ------------------------------
        values$cached_samples_tab <- hot_to_r(input$samples_tab)
        
        current_samples  <- hot_to_r(input$raw_data)$Sample |> unique() |> drop_empty()
        
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
            # reapply previous edits for existing samples
            values$cached_samples_tab <- data.frame(Sample = current_samples) |>
                left_join(values$cached_samples_tab) |>
                # fill in defaults for new samples
                mutate(
                    New_Label = coalesce(New_Label, Sample),
                    Include   = coalesce(Include, TRUE)
                )
        }
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
            hot_col("Sample", readOnly = TRUE) |>
            hot_col("New_Label", type = "text") |>
            hot_col("Include", type = "checkbox") |>
            hot_cols(columnSorting = FALSE)
            
    })

    
    
    # -------------------------------------------------------------------------
    # Derived Reactive: Processed data (with parsed Cq, renames, sample ordering and exclusions)
    # -------------------------------------------------------------------------
    
    cq_data <- reactive({
        req(hot_to_r(input$raw_data))
        req(nrow(hot_to_r(input$samples_tab)) > 0)

        samples_metadata <- hot_to_r(input$samples_tab)
        
        hot_to_r(input$raw_data) |>
            mutate(Key = row_number()) |> # add unique Key ID matching raw data rows
            relocate(Key, .before = 1) |>
            # join with sample metadata for renaming, reordering and exclusion
            inner_join(samples_metadata) |>
            filter(Include) |>
            # update and reorder sample name
            mutate(Sample = factor(New_Label,
                                   levels = samples_metadata$New_Label)) |>
            select(-New_Label, -Include) |>
            arrange(Sample) |>
            mutate(Cq = parse_Cq(Cq) |> as.numeric(),
                   Cq_no_und = ifelse(Cq == Inf, NA, Cq)) |>
        # mark excluded points
            mutate(
                Keep = !Key %in% values$excluded_point_keys
            )
    })
    
    # -------------------------------------------------------------------------
    # Observer: Update target selector choices
    # -------------------------------------------------------------------------

    observe({
        req(cq_data())
        targets <- cq_data()$Target |>
            unique() |>
            drop_empty()

        updateSelectInput(session, "target_select", choices = targets, selected = targets[1])
    })
    
    # -------------------------------------------------------------------------
    # Output: Plot and plot title
    # -------------------------------------------------------------------------
    
    output$ct_plot_title <- renderText({
        req(input$target_select)
        paste("Cq Values for", input$target_select)
    })
    
    output$ct_plot <- renderPlotly({
        df <- cq_data()
        req(df, nrow(df) > 0)
        req(input$target_select)
        
        df_target <- df |>
            filter(Target == input$target_select) |>
            mutate(Keep_label = ifelse(Keep, "Included", "Excluded"))
        
        df_summary_target <- df_target |>
            filter(Keep) |>
            group_by(across(
                c("Sample", "Target", any_of("Replicate"))
            )) |>
            summarize(
                mean   = mean(Cq_no_und, na.rm = TRUE),
                median = median(Cq_no_und, na.rm = TRUE),
            )
        
        # force a minumum of y-axis range of 3 units
        y_limits <- get_y_limits(df_target$Cq_no_und, min_range = 3)
        
        
        p <- ggplot(df_target,
                    aes(x = Sample, y = Cq_no_und,
                        alpha = Keep_label,
                        # label on hoover
                        text = paste0(
                            Sample, 
                            ifelse("Replicate" %in% names(df_target),
                                    paste0(" (", Replicate, ")"),
                                    ""
                            ), "\n",
                            "Target: ", Target, "\n",
                            "Cq: ", round(Cq_no_und, 2),
                            ifelse(Keep, "", " (excluded)")
                            ),
                        key = Key
                        )) +
            geom_quasirandom(
                color = secondary_color(),
            ) +
            scale_alpha_manual(
                values = c("Included" = 1, "Excluded" = 0.3),
                name = "",
            ) +
            # add mean/median points
            geom_point(data = df_summary_target, 
                       aes(x = Sample, y = .data[[input$aggr_function]],
                           shape = str_to_title(input$aggr_function)),
                       inherit.aes = F,
                       size = 4, color = accent_color()
                       ) +
            scale_shape_manual(
                values = 4,
                name = ""
            ) +
            labs(
                x = "Sample",
                y = "Cq",
                title = NULL
            ) +
            scale_y_continuous(expand = c(0.1)) +
            theme_minimal(base_size = 14) +
            theme(
                legend.position = "bottom",
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1)
            )
        
        # facet by replicate if present
        if ("Replicate" %in% names(df)) {
            p <- p + facet_wrap(~ Replicate)
        }
        
        ggplotly(p, tooltip = "text", source = "ct_plot") |>
            event_register("plotly_click") |>
            fix_plotly_legend()
        })
    
    # -------------------------------------------------------------------------
    # Observer: Handle click-to-exclude
    # -------------------------------------------------------------------------
    
    observeEvent(event_data("plotly_click", source = "ct_plot"), {
        click <- event_data("plotly_click", source = "ct_plot")
        req(click)
        
        clicked_key <- click$key
        
        # toggle inclusion/exclusion
        if(clicked_key %in% values$excluded_point_keys) {
            # currently excluded, include it (remove it from the exclusion list)
            values$excluded_point_keys <- setdiff(values$excluded_point_keys, clicked_key)
        } else {
            # currently included, exclude it (add it to the exclusion list)
            values$excluded_point_keys <- append(clicked_key, values$excluded_point_keys)
        }
    })
    
}



# =============================================================================
# Run App
# =============================================================================

shinyApp(ui = ui, server = server)
