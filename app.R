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
library(bsicons)

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
    Target = rep(rep(c("TBP", "Target1", "Actin", "Target2"), each = 3), 3),
    Cq     = round(rnorm(3 * 3 * 4, mean = 23.5, sd = 1.2), 2) |> as.character()
)


# helper functions ==========================================
drop_empty <- function(x) {x[!is.na(x) & nzchar(trimws(x))]}

source("R/parse_Cq.R")
source("R/get_y_limits.R")
source("R/fix_plotly_legend.R")
source("R/is_HK.R")

# =============================================================================
# UI Definition
# =============================================================================

ui <- page_fillable(
    # Custom CSS
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
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
            page_sidebar(
                fillable = TRUE,
                sidebar = sidebar(
                    title = "Sample Controls",
                    open = TRUE,
                    width = "380px",
                    prettySwitch(
                        "include_replicates",
                        label = "Include biological replicates column",
                        fill  = TRUE, status = "primary",
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
                        class = "d-flex justify-content-between",
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
                            class = "btn-danger btn-sm"
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
            page_sidebar(
                fillable = TRUE,
                sidebar = sidebar(
                    title = "Cq Inspection",
                    open = TRUE,
                    width = "380px",
                    selectInput(
                        "ct_target_select",
                        "Select Target",
                        choices = NULL # popolate dinamically
                    ),
                    br(),
                ),
                # Main content area
                card(
                    full_screen = TRUE,
                    fillable = TRUE,
                    card_header(textOutput("ct_plot_title", inline = TRUE)),
                    plotlyOutput("ct_plot", height = "100%")
                )
            )
        ),
        
        # ---------------------------------------------------------------------
        # Panel 3: dCq analysis
        # ---------------------------------------------------------------------
        nav_panel(
            title = "ΔCq",
            page_sidebar(
                fillable = TRUE,
                sidebar = sidebar(
                    title = "ΔCq Settings",
                    open = TRUE,
                    width = "380px",
                    h5("HK genes"),
                    pickerInput(
                        inputId = "hk_genes",
                        label = "Select HK gene(s):", 
                        choices = NULL, # populate dynamicallyLL,
                        multiple = TRUE,
                        options = pickerOptions(container = "body", 
                                                actionsBox = TRUE),
                        width = "100%"
                    ),
                    hr(),
                    h5("Visualization"),
                    selectInput(
                        "dCt_target_select",
                        "Select Target",
                        choices = NULL # popolate dinamically
                    ),
                    div(
                        #style = "display: flex; gap: 5px",
                        class = "d-flex justify-content-start",
                        prettySwitch(
                            "dCq_stats",
                            label = "Show descriptive statistics",
                            fill  = TRUE, status = "primary",
                            value = FALSE
                        ),
                        tooltip(
                            bs_icon("info-circle"),
                            "Technical replicates: SD (and SEM) here are shown for QC only (assay/measurement precision).
                            They should not be used for statistical testing or to infer biological variability.",
                            placement = "right",
                            
                        )
                    ),
                    
                    conditionalPanel(
                        condition = "input.dCq_stats === true",
                        div(
                            #style = "display: flex; gap: 5px",
                            class = "d-flex justify-content-start",
                            radioGroupButtons(
                                inputId = "dCq_stat_type",
                                label = "Summarize by:",
                                choices = c("SD", "SEM"),
                                justified = TRUE,
                                width = "120px",
                                size = "sm"
                            ),
                            tooltip(
                                bs_icon("info-circle"),
                                "For technical replicates, the preferred display is the individual points only (no SD/SEM) to avoid confusion with biological variability.
                                If an error bar is shown, SD is more common because it reflects the spread/precision of the assay; SEM reflects precision of the mean (SD/√n) and shrinks with n.",
                                placement = "right"
                            )
                        ),
                        
                        div(
                            style = "display: flex; gap: 5px",
                            class = "d-flex justify-content-start",
                            prettySwitch(
                                "dCq_propagate_hk_var",
                                label = "Propagate HK variance",
                                fill  = TRUE, status = "primary",
                                value = TRUE
                            ),
                            tooltip(
                                bs_icon("info-circle"),
                                "Default: variability is computed using both Target and HK technical replicate variation (HK variation is pooled across different HK and propagated into ΔCq).
                                This is more faithful to the measurement process.
                                “Target-only” stats are offered for reproducibility with common practice, not because they’re more correct, just simpler to calculate.",
                                placement = "right"
                            )
                        ),
                        
                    ),
                    br()
                ),
                # Main content area
                card(
                    full_screen = TRUE,
                    fillable = TRUE,
                    card_header(textOutput("dCt_plot_title", inline = TRUE)),
                    plotlyOutput("dCt_plot", height = "100%")
                )
            )
        ),
        
        # ---------------------------------------------------------------------
        # Panel 4: Results
        # ---------------------------------------------------------------------
        nav_panel(
            title = "Results",
            page_sidebar(
                fillable = TRUE,
            )
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

        updateSelectInput(session, "ct_target_select", choices = targets, selected = targets[1])
        
        selected_hk <- targets[is_HK(targets)]
        
        updatePickerInput(session, "hk_genes", choices = targets, selected = selected_hk)
    })
    
    # -------------------------------------------------------------------------
    # Output: Cq Plot and plot title
    # -------------------------------------------------------------------------
    
    output$ct_plot_title <- renderText({
        req(input$ct_target_select)
        paste("Cq Values for", input$ct_target_select)
    })
    
    output$ct_plot <- renderPlotly({
        df <- cq_data()
        req(df, nrow(df) > 0)
        req(input$ct_target_select)
        
        df_target <- df |>
            filter(Target == input$ct_target_select) |>
            mutate(Keep_label = ifelse(Keep, "Included", "Excluded"))
        
        df_summary_target <- df_target |>
            filter(Keep) |>
            group_by(across(
                c("Sample", "Target", any_of("Replicate"))
            )) |>
            summarize(
                mean   = mean(Cq_no_und, na.rm = TRUE)
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
                       aes(x = Sample, y = mean,
                           shape = "Mean"),
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
    
    # -------------------------------------------------------------------------
    # Derived Reactive: dCq
    # -------------------------------------------------------------------------
    
    dCq_data <- reactive({
        req(cq_data())
        req(nrow(cq_data()) > 0)
        req(length(input$hk_genes) > 0)
        
        HK_summary <- cq_data() |>
            filter(Keep,
                   Target %in% input$hk_genes) |>
            group_by(across(
                c("Sample", "Target", any_of("Replicate"))
                )) |>
            # summarize each HK separately
            summarize(
                HK_mean = mean(Cq_no_und, na.rm = TRUE),
                HK_sd   = sd(Cq, na.rm = TRUE),
                HK_n    = length(na.omit(Cq)), # sample size for each gene
                .groups = "drop"
            ) |>
            # summarize all HKs per sample
            group_by(across(
                c("Sample",  any_of("Replicate"))
            )) |>
            summarize(
                HK_mean      = mean(HK_mean, na.rm = TRUE),
                n_HK_genes   = sum(HK_n > 0),
                # pooled SD for independet samples, allowing different mean (same as in ANOVA)
                HK_sd_pool   = sqrt(
                    sum(HK_sd^2*(HK_n-1), na.rm = T)/
                        (sum(HK_n,  na.rm = T)-n_HK_genes)
                ), 
                HK_se_pooled = HK_sd_pool*sqrt(1/(sum(HK_n, na.rm = T))),
                .groups = "drop"
            )
        
        cq_data() |>
            filter(Keep,
                   !Target %in% input$hk_genes) |>
            left_join(HK_summary) |>
            mutate(
                dCq    = Cq_no_und - HK_mean,
            )
    })
    
    # -------------------------------------------------------------------------
    # Derived Reactive: dCq summary per individual replicate
    # -------------------------------------------------------------------------
    
    dCq_rep_summary <- reactive({
        req(dCq_data())
        req(nrow(dCq_data()) > 0)
        
        dCq_data() |>
            group_by(across(
                c("Sample", "Target", any_of("Replicate"))
            )) |>
            summarize(
                Cq_n   = length(na.omit(dCq)),
                Cq_sd  = sd(Cq_no_und, na.rm = TRUE),
                Cq_se  = Cq_sd/sqrt(Cq_n),
                
                dCq_mean = mean(dCq, na.rm = TRUE),
                # propagate SD and SE including HK variance/uncertainty.
                dCq_sd = sqrt(Cq_sd^2 + HK_sd_pool^2),
                dCq_se = sqrt(Cq_sd^2 + HK_sd_pool^2),
                
                # Compute Cq stats (SD only), without propagating HK variance/uncertainty.
                # This is statistically wrong but common in practice because it’s straightforward to compute.
                dCq_sd_no_p = sd(dCq, na.rm = TRUE),
                dCq_se_no_p = dCq_sd_no_p/sqrt(Cq_n),
                
                
                .groups = "drop"
            )
    })
    
    # -------------------------------------------------------------------------
    # Output: dCq Plot and plot title
    # -------------------------------------------------------------------------
    observeEvent(dCq_data(), {
        req(dCq_data(), nrow(dCq_data()) > 0)
        
        non_hk_genes <- dCq_data()$Target |> unique() |> drop_empty()
        
        updateSelectInput(session, "dCt_target_select", choices = non_hk_genes, selected = non_hk_genes[1])
    })
    
    output$dCt_plot_title <- renderText({
        req(input$dCt_target_select)
        paste("ΔCq Values for", input$dCt_target_select)
    })
    
    output$dCt_plot <- renderPlotly({
        df <- dCq_data()
        df_summary <- dCq_summary()
        req(df, nrow(df) > 0)
        req(df_summary, nrow(df_summary) > 0)
        req(input$dCt_target_select)
        
        df_target <- df |>
            filter(Target == input$dCt_target_select)
        
        df_summary_target <- df_summary |>
            filter(Target == input$dCt_target_select)
        
        # force a minumum of y-axis range of 3 units
        y_limits <- get_y_limits(df_target$dCq, min_range = 3)
        
        p <- ggplot(df_target,
                    aes(x = Sample, y = -dCq,
                        # label on hoover
                        text = paste0(
                            Sample, 
                            ifelse("Replicate" %in% names(df_target),
                                   paste0(" (", Replicate, ")"),
                                   ""
                            ), "\n",
                            "Target: ", Target, "\n",
                            "ΔCq: ", round(dCq, 2)
                        )
                    )) +
            geom_quasirandom(
                color = secondary_color(),
            ) +
            # add mean points
            geom_point(data = df_summary_target, 
                       aes(x = Sample, y = -dCq_mean,
                           shape = "Mean"),
                       inherit.aes = F,
                       size = 4, color = accent_color()
            ) +
            scale_shape_manual(
                values = 4,
                name = ""
            ) +
            labs(
                x = "Sample",
                y = "-ΔCq",
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
        
        ggplotly(p, tooltip = "text") |>
            fix_plotly_legend()
    })
}


# =============================================================================
# Run App
# =============================================================================

shinyApp(ui = ui, server = server)
