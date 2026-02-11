library(shiny)
library(bslib)
library(shinyWidgets)
library(bsicons)
library(rhandsontable)

library(tidyverse)
library(plotly)
library(ggbeeswarm)
library(scales)
library(glue)
library(DT)

# Example Data =================================================================

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


# helper functions =============================================================
drop_empty <- function(x) {
    x[!is.na(x) & nzchar(str_trim(x))]
}

source("R/parse_Cq.R")
source("R/get_y_limits.R")
source("R/fix_plotly_legend.R")
source("R/is_HK.R")
source("R/handle_undetected_stats.R")
source("R/squish_infinite_to_val.R")
source("R/statistical_tests.R")

# UI Definition ================================================================

ui <- page_fillable(
    # Custom CSS
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "www/custom.css")
    ),

    # Main navigation
    navset_tab(
        id = "main_tabs",
        nav_item(h5("qPCR Analysis Tool")),
        nav_spacer(),
        # Panel 1: Input Data --------------------------------------------------
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
                        fill = TRUE, status = "primary",
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
        # Panel 2: Cq values inspection and exclusion --------------------------
        nav_panel(
            title = "Cq",
            page_sidebar(
                fillable = TRUE,
                sidebar = sidebar(
                    title = "Cq Inspection",
                    open = TRUE,
                    width = "380px",
                    selectInput(
                        "select_ct_target",
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
        # Panel 3: dCq analysis ------------------------------------------------
        nav_panel(
            title = "Results",
            page_sidebar(
                fillable = TRUE,
                sidebar = sidebar(
                    title = "Settings",
                    open  = TRUE,
                    width = "380px",
                    pickerInput(
                        inputId = "hk_genes",
                        label = "Select HK gene(s)",
                        choices = NULL, # populate dynamicallyLL,
                        multiple = TRUE,
                        options = pickerOptions(
                            container  = "body",
                            actionsBox = TRUE
                        ),
                        width = "100%"
                    ),
                    selectInput(
                        "select_out_target",
                        "Select Target to visualize",
                        choices = NULL # popolate dinamically
                    ),
                    hr(),
                    radioGroupButtons(
                        inputId   = "summarize_bio_reps",
                        label     = "Biological Replicates",
                        choices   = c("Plot Individually" = "split", "Aggregate" = "aggregate"),
                        justified = TRUE,
                        width     = "100%",
                        size      = "sm"
                    ),
                    radioGroupButtons(
                        inputId  = "out_metric",
                        label    = "Plot",
                        choices  = c("-ΔCq" = "dCq", "2^-ΔCq" = "exp_dCq", "-ΔΔCq" = "ddCq", "2^-ΔΔCq" = "exp_ddCq"),
                        justified = TRUE,
                        disabled  = c("ddCq", "exp_ddCq"), # enable dinamically
                        width     = "100%",
                        size      = "sm"
                    ),
                    div(
                        # style = "display: flex; gap: 5px",
                        class = "d-flex justify-content-start align-items-center",
                        radioGroupButtons(
                            inputId   = "stat_type",
                            label     = "Error Bars:",
                            choices   = c("SEM" = "se", "SD" = "sd", "None" = "none"),
                            selected  = "none", # updated dinamically
                            justified = TRUE,
                            width     = "180px",
                            size      = "sm"
                        ),
                        conditionalPanel(
                            condition = "input.summarize_bio_reps == 'split' && input.stat_type != 'none'",
                            tooltip(
                                bs_icon("exclamation-triangle"),
                                "For technical replicates error bars are shown as QC only.
                                They should not be used for statistical testing or to infer biological variability but rather to decide if keep or repeat the assay.
                                For technical replicates, the preferred display is 'None' (individual points only) to avoid confusion with biological variability.
                                If you want to add an error bar, SD is a more common choice because it reflects the variability of the assay; SEM reflects precision of the mean estimation (SD/√n) and shrinks with n.",
                            )
                        )
                    ),
                    conditionalPanel(
                        condition = "input.summarize_bio_reps == 'split' && input.stat_type != 'none'",
                        div(
                            style = "display: flex; gap: 5px",
                            class = "d-flex justify-content-start",
                            prettySwitch(
                                "propagate_var",
                                label   = "Propagate variance",
                                fill    = TRUE,
                                status  = "primary",
                                value   = TRUE
                            ),
                            tooltip(
                                bs_icon("info-circle"),
                                "Default: variance is computed from Target + HK technical replicates (HK pooled and propagated into ΔCq). For ΔΔCq, control variance is propagated to all samples (not ignored).
                                This is more faithful to the measurement process.
                                Stats without propagation of variance are offered for reproducibility with common practice. Note that this is not the correct approach though",
                            )
                        ),
                    ),
                    # Statistical Analysis Section (only if n_bio_reps > 1 and n_samples >= 2) ----------
                    conditionalPanel(
                        condition = "output.n_bio_reps >= 2 && output.n_samples >= 2",
                        hr(),
                        tags$h6(tags$strong("Statistical Analysis")),

                        # Target metric for statistical testing
                        radioGroupButtons(
                            inputId   = "stats_metric",
                            label     = "Test on:",
                            choices   = c(
                                "-ΔCq"    = "dCq",
                                "-ΔΔCq"   = "ddCq",
                                "2^-ΔΔCq" = "exp_ddCq"
                            ),
                            selected  = "dCq",
                            justified = TRUE,
                            width     = "100%",
                            size      = "sm"
                        ),

                        # Warning for exp data
                        conditionalPanel(
                            condition = "input.stats_metric == 'exp_ddCq'",
                            div(
                                class = "alert alert-warning py-1 px-2 mb-2",
                                style = "font-size: 0.85em;",
                                bs_icon("exclamation-triangle"),
                                "It is recommended to perform statistical analysis on the log space (-ΔCq or -ΔΔCq) values rather than on exponentiated ones for more reliable results."
                            )
                        ),

                        # Omnibus test selection (choices update dynamically via pickerInput with optgroups)
                        div(
                            class = "d-flex justify-content-start align-items-center mb-2",
                            pickerInput(
                                inputId   = "stats_test",
                                label     = "Test:",
                                choices   = list(
                                    "Parametric" = c(
                                        "ANCOVA" = "ancova",
                                        "Mixed Effect Model" = "mixed_effect",
                                        "Pairwise paired t-test" = "pairwise_paired_ttest"
                                    ),
                                    "Non-parametric" = c( # "Friedman" = "friedman",
                                        "Pairwise Wilcoxon" = "pairwise_wilcoxon"
                                    )
                                ), # update dinamically
                                selected = "ancova",
                                width = "93%",
                                options = pickerOptions(container = "body")
                            ),

                            # Tip about test recomendations (shown for dCq with > 2 samples)
                            conditionalPanel(
                                condition = "input.stats_metric == 'dCq' && output.n_samples > 2",
                                tooltip(
                                    bs_icon("lightbulb"),
                                    "ANCOVA is recommended when you have a clear reference/control sample (e.g., untreated sample).
                                    Mixed-effect Model is better when the reference sample is arbitrary across replicates (e.g., comparing expression between different patients or cell lines).
                                    Other tests are available for completeness but are generally not recommended.",
                                    placement = "right"
                                )
                            ),
                            # Tip about test recomendations (shown for dCq with exactly 2 samples)
                            conditionalPanel(
                                condition = "input.stats_metric == 'dCq' && output.n_samples == 2",
                                tooltip(
                                    bs_icon("lightbulb"),
                                    "ANCOVA is recommended when you have a clear reference/control sample (e.g., untreated vs treated). Paired t-test is better when the reference sample is arbitrary across replicates (e.g. comparing the expression between 2 different tumors or cell lines).",
                                    placement = "right"
                                )
                            )
                        ),

                        # Handle unequal variance toggle (for mixed effect model)
                        conditionalPanel(
                            condition = "output.show_unequal_variance_toggle",
                            div(
                                class = "d-flex justify-content-start align-items-center",
                                prettySwitch(
                                    inputId = "stats_unequal_variance",
                                    label   = "Handle unequal variance",
                                    fill    = TRUE,
                                    status  = "primary",
                                    value   = FALSE
                                ),
                                # Tip shown only for pairwise t-test
                                conditionalPanel(
                                    condition = "input.stats_test == 'pairwise_ttest'",
                                    tooltip(
                                        bs_icon("info-circle"),
                                        "When enabled (Welch's t-test), 
                                        the comparison against the reference sample converges to a one-sample t-test 
                                        since the reference has zero variance in ΔΔCq.",
                                        placement = "right"
                                    )
                                )
                            )
                        ),

                        # Post Hoc comparison type (only shown when > 2 samples)
                        conditionalPanel(
                            condition = "output.show_multiple_comparison_type",
                            div(
                                class = "stats-radio-compact",
                                radioGroupButtons(
                                    inputId   = "stats_comparison",
                                    label     = "Multiple comparisons:",
                                    choices   = c(
                                        "Pairwise" = "pairwise",
                                        "All vs Control" = "trt.vs.ctrl"
                                    ),
                                    selected  = "pairwise",
                                    justified = TRUE,
                                    width = "100%",
                                    size = "sm"
                                )
                            )
                        ),

                        # Multiple comparison adjustment method
                        # let chose between Benjamini-Hochberg (FDR), Holm (FWER),  none
                        conditionalPanel(
                            condition = "output.show_multiple_comparison_adjust",
                            div(
                                class = "stats-radio-compact",
                                radioGroupButtons(
                                    inputId   = "stats_multiple_comparison_adjust",
                                    label     = "Multiple comparison adjustment method:",
                                    choices   = c(
                                        "Benjamini-Hochberg (FDR)" = "BH",
                                        "Holm (FWER)" = "holm",
                                        "None" = "none"
                                    ),
                                    selected  = "BH",
                                    justified = TRUE,
                                    width     = "100%",
                                    size      = "sm",
                                    direction = "vertical"
                                )
                            )
                        ),

                        # Post-hoc test display (dynamically updated)
                        conditionalPanel(
                            condition = "output.show_post_hoc_test",
                            div(
                                class = "card bg-light border-0",
                                div(
                                    class = "card-body py-2 px-3",
                                    div(
                                        class = "d-flex align-items-center gap-2",
                                        bs_icon("arrow-return-right", class = "text-primary"),
                                        tags$span(class = "text-muted small", "Post-hoc:"),
                                        tags$span(class = "fw-semibold", textOutput("stats_posthoc", inline = TRUE))
                                    )
                                )
                            )
                        )
                    )
                ),
                # Main content area
                card(
                    full_screen = TRUE,
                    fillable = TRUE,
                    card_header(textOutput("res_plot_title", inline = TRUE)),
                    plotlyOutput("res_plot", height = "100%")
                ),
                # Statistical Results Card (only shown when stats panel is active)
                conditionalPanel(
                    condition = "output.n_bio_reps >= 2 && output.n_samples >= 2",
                    card(
                        card_header(
                            textOutput("stats_card_title", inline = TRUE)
                        ),
                        # Warning for dropped Inf values
                        conditionalPanel(
                            condition = "output.stats_dropped_count > 0",
                            div(
                                class = "alert alert-warning py-2 px-3 mb-3",
                                style = "font-size: 0.85em;",
                                bs_icon("exclamation-triangle"),
                                textOutput("stats_dropped_warning", inline = TRUE)
                            )
                        ),
                        # Omnibus section (for ANCOVA, ANOVA, Mixed Effect, Kruskal-Wallis)
                        conditionalPanel(
                            condition = "output.has_omnibus_test",
                            tags$h6(
                                #class = "mb-2",
                                "Omnibus test:"
                            ),
                            # omnibus title:
                            accordion(
                                id = "omnibus_accordion",
                                class = "mb-2",
                                accordion_panel(
                                    title = div(
                                        class = "d-flex justify-content-between align-items-center",
                                        span(class = "mx-3", style = "font-size: 14px;", textOutput("stats_omnibus_label", inline = TRUE)),
                                        uiOutput("stats_omnibus_badge")
                                    ),
                                    value = "omnibus_panel",
                                    icon = NULL,
                                    DT::dataTableOutput("stats_omnibus_table")
                                )
                            )
                        ),
                        # Post-hoc or Pairwise results section
                        div(
                            tags$h6(
                                class = "mb-2",
                                textOutput("stats_comparison_title", inline = TRUE)
                            ),
                            DT::dataTableOutput("stats_results_table")
                        ),
                        # Method description
                        div(
                            class = "bg-light py-2 px-3",
                            tags$strong("Methods: "),
                            textOutput("stats_method", inline = TRUE)
                        )
                    )
                )
            )
        )
    )
)

# Server Logic =================================================================

server <- function(input, output, session) {
    # Current theme accent color -----------------------------------------------
    accent_color <- reactive({
        bslib::bs_get_variables(bslib::bs_current_theme(session), "primary")[[1]]
    })

    secondary_color <- reactive({
        bslib::bs_get_variables(bslib::bs_current_theme(session), "secondary")[[1]]
    })
    # Cache (reactiveValues) ---------------------------------------------------

    cache <- reactiveValues(
        # raw data, required to add/remove replicate column
        raw_data = empty_raw_data,

        # sample control table (rename, reorder, exclude)
        # samples rename, reordering and exclusion should be retrieved from input$samples_tab
        samples_tab = data.frame(
            Sample    = character(),
            New_Label = character(),
            Include   = logical(),
            stringsAsFactors = FALSE
        ),

        # list of excluded points
        excluded_point_keys = c(),
        select_ct_targeted  = c(),
        targets_available   = c()
    )
    # Observer: Toggle biological replicates column ----------------------------

    observeEvent(input$include_replicates, {
        if (is.null(input$raw_data)) {
            return()
        }

        current_data <- hot_to_r(input$raw_data)

        if (input$include_replicates) {
            # add replicate column if missing
            if (!"Replicate" %in% names(current_data)) {
                cache$raw_data <- current_data |>
                    mutate(Replicate = "R1")
            }
        } else {
            # drop replicate column if present
            cache$raw_data <- select(current_data, -any_of("Replicate"))
        }
    })
    # Observer: Load example data ----------------------------------------------

    observeEvent(input$load_example, {
        # Wait for rhandsontable to be initialized to prevent infinite loop
        req(input$raw_data)
        
        # Force-enable biological replicates when loading example data
        updatePrettySwitch(session, "include_replicates", value = TRUE)

        # Load example data from CSV
        cache$raw_data <- read_csv("data/example_qPCR_data.csv", show_col_types = FALSE)
    })
    # Observer: clean data -----------------------------------------------------

    observeEvent(input$clear_data, {
        if (input$include_replicates) {
            cache$raw_data <- empty_raw_data |>
                mutate("Replicate" = character(5))
        } else {
            cache$raw_data <- empty_raw_data
        }

        # reset the lsit of excluded points
        cache$excluded_point_keys <- c()
    })
    # Output: Raw data table ---------------------------------------------------

    output$raw_data <- renderRHandsontable({
        req(cache$raw_data)

        rhandsontable(
            cache$raw_data,
            rowHeaders = TRUE,
            readOnly    = FALSE,
            contextMenu = TRUE,
            stretchH    = "all",
            renderAllRows = TRUE
        ) |>
            hot_col("Sample", type = "text") |>
            hot_col("Target", type = "text") |>
            hot_col("Cq", type = "text") |> # type = "text" to allow for Undetermined or other labels
            hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
    })

    # Observer: on raw data edit: ----------------------------------------------
    # 1. cache last `raw_data` and validate Cq values
    # 2. update `cache$samples_tab` when samples change in raw_data while preserving previous edits

    observeEvent(input$raw_data, {
        # 1. validate Cq values ------------------------------------------------
        current_data <- hot_to_r(input$raw_data)

        original_cq <- current_data$Cq |>
            as.character() |>
            str_trim() |>
            replace_na("") |>
            # ignore trailing 0 in decimals
            str_remove("(?<=[0-9]\\.[0-9]{0,10})0+$")

        parsed_cq <- parse_Cq(original_cq) |> replace_na("")

        # Find values that changed
        changed_mask <- original_cq != parsed_cq

        conversions <- map2_chr(
            original_cq[changed_mask], parsed_cq[changed_mask],
            ~ paste0(.x, " → ", .y)
        ) |>
            unique()

        # Update parsed Cq values
        current_data$Cq <- parsed_cq

        # Cache updated raw data
        cache$raw_data <- current_data

        # Show warning modal if any conversions happened
        if (length(conversions) > 0) {
            showModal(modalDialog(
                title = "Cq Values Converted",
                tags$p("The following conversions were applied:"),
                tags$ul(
                    lapply(conversions, function(x) {
                        tags$li(x)
                    })
                ),
                easyClose = TRUE,
                footer = modalButton("OK")
            ))
        }

        # 2. cache last state of samples_tab -----------------------------------
        cache$samples_tab <- hot_to_r(input$samples_tab)

        current_samples <- hot_to_r(input$raw_data)$Sample |>
            unique() |>
            drop_empty()

        previous_samples <- cache$samples_tab$Sample
        new_samples <- setdiff(current_samples, previous_samples)

        if (length(current_samples) == 0) {
            # No valid samples, keep empty
            cache$samples_tab <- data.frame(
                Sample    = character(),
                New_Label = character(),
                Include   = logical()
            )
        } else if (length(previous_samples) == 0) {
            # First time samples_tab update
            cache$samples_tab <- data.frame(
                Sample    = current_samples,
                New_Label = current_samples,
                Include   = rep(TRUE, times = length(current_samples))
            )
        } else {
            # reapply previous edits for existing samples
            cache$samples_tab <- data.frame(Sample = current_samples) |>
                left_join(cache$samples_tab) |>
                # fill in defaults for new samples
                mutate(
                    New_Label = coalesce(New_Label, Sample),
                    Include   = coalesce(Include, TRUE)
                )
        }
    })
    # Output: Sample control table ---------------------------------------------
    output$samples_tab <- renderRHandsontable({
        req(nrow(cache$samples_tab) > 0)

        rhandsontable(
            cache$samples_tab,
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
    # Derived Reactive: Processed data (with parsed Cq, renames, sample ordering and exclusions) ----------

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
                levels = unique(samples_metadata$New_Label)
            )) |>
            select(-New_Label, -Include) |>
            arrange(Sample) |>
            mutate(Cq = parse_Cq(Cq) |> as.numeric()) |>
            # mark excluded points
            mutate(
                Keep = !Key %in% cache$excluded_point_keys,
                Undetected = !is.finite(Cq)
            )
    })
    # Observer: Update target selector choices ---------------------------------

    observe({
        req(cq_data())

        targets <- cq_data()$Target |>
            unique() |>
            drop_empty() |>
            sort()

        # if no change in targets, skip update
        req(!identical(targets, cache$targets_available))

        # restore previous selection if possible
        if (cache$select_ct_targeted %in% targets) {
            selected <- cache$select_ct_targeted
        } else {
            selected <- targets[1]
        }

        # cache last list of targets
        cache$targets_available <- targets
        updateSelectInput(session, "select_ct_target", choices = targets, selected = selected)

        selected_hk <- targets[is_HK(targets)]
        updatePickerInput(session, "hk_genes", choices = targets, selected = selected_hk)
    })
    # OvserveEvent: cache last target selected ---------------------------------
    observeEvent(input$select_ct_target, {
        cache$select_ct_targeted <- input$select_ct_target
    })
    # Output: Cq Plot and plot title -------------------------------------------

    output$ct_plot_title <- renderText({
        req(input$select_ct_target)
        paste("Cq Values for", input$select_ct_target)
    })

    output$ct_plot <- renderPlotly({
        df <- cq_data()
        req(df, nrow(df) > 0)
        req(input$select_ct_target)

        df_target <- df |>
            filter(Target == input$select_ct_target) |>
            mutate(
                Keep_label = ifelse(Keep, "Included", "Excluded"),
                point_type_label = ifelse(Undetected, "Undetected", "Detected"),
            )

        n_samples <- df_target$Sample |>
            unique() |>
            length()

        df_summary_target <- df_target |>
            filter(Keep) |>
            group_by(across(
                c("Sample", "Target", any_of("Replicate"))
            )) |>
            summarize(
                mean = mean_handle_inf(Cq),
                point_type_label = "Mean",
                Keep_label = NA # initialize empty keep_label to avoid duplication of the legend in ggplotly
            )

        # force a minumum of y-axis range of 3 units
        y_limits <- get_Cq_y_limits(df_target$Cq, min_range = 3)

        # Pre-compute hover text
        has_replicate <- "Replicate" %in% names(df_target)
        
        df_target <- df_target |>
            mutate(text = if (has_replicate) {
                glue(
                    "{Sample} ({Replicate})
                    Target: {Target}
                    Cq: {round(Cq, 2)}
                    {ifelse(Keep, '', '(excluded)')}"
                )
            } else {
                glue(
                    "{Sample}
                    Target: {Target}
                    Cq: {round(Cq, 2)}
                    {ifelse(Keep, '', '(excluded)')}"
                )
            })
        
        df_summary_target <- df_summary_target |>
            mutate(text = if (has_replicate) {
                glue(
                    "{Sample} ({Replicate})
                    Target: {Target}
                    Mean Cq: {round(mean, 2)}"
                )
            } else {
                glue(
                    "{Sample}
                    Target: {Target}
                    Mean Cq: {round(mean, 2)}"
                )
            })

        p <- ggplot(
            df_target,
            aes(
                x = Sample, y = Cq,
                alpha = Keep_label,
                color = point_type_label,
                shape = point_type_label,
                text = text,
                key = Key
            )
        ) +
            annotate("rect", xmin = 0.5, xmax = n_samples + 0.5, ymin = 35, ymax = 40, alpha = 0.6, fill = "#EBEBEB") +
            geom_beeswarm(method = "compactswarm", preserve.data.axis = TRUE) +
            geom_point(
                data = df_summary_target,
                aes(
                    x = Sample, y = mean,
                    text = text,
                    shape = point_type_label,
                    color = point_type_label,
                    alpha = Keep_label
                ),
                inherit.aes = F,
                size = 4
            ) +
            scale_shape_manual(
                values = c("Detected" = 16, "Undetected" = 1, "Mean" = 4),
                name = "",
            ) +
            scale_color_manual(
                values = c("Detected" = secondary_color(), "Undetected" = "#C03A2B", "Mean" = accent_color()),
                name = "",
            ) +
            scale_alpha_manual(
                values = c("Included" = 1, "Excluded" = 0.3),
                na.value = 1,
                name = "",
            ) +
            labs(
                x = "Sample",
                y = "Cq",
                title = NULL
            ) +
            coord_cartesian(ylim = y_limits) +
            scale_y_continuous(
                expand = expansion(mult = 0.05, add = 0),
                labels = function(x) ifelse(x == 40, "≥40", x), # label 40+ for undetected
                # add extra distance to squish to separate from other points
                oob = squish_infinite_to_val
            ) +
            scale_x_discrete(expand = 0) +
            theme_minimal(base_size = 14) +
            theme(
                legend.position = "bottom",
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1)
            )

        # facet by replicate if present
        if ("Replicate" %in% names(df)) {
            p <- p + facet_wrap(~Replicate)
        }

        ggplotly(p, tooltip = "text", source = "ct_plot") |>
            event_register("plotly_click") |>
            fix_plotly_legend()
    })
    # Observer: Handle click-to-exclude ----------------------------------------

    observeEvent(event_data("plotly_click", source = "ct_plot"), {
        click <- event_data("plotly_click", source = "ct_plot")
        req(click)
        # ignore click on summary points (no key column)
        req("key" %in% names(click))

        clicked_key <- click$key

        # toggle inclusion/exclusion
        if (clicked_key %in% cache$excluded_point_keys) {
            # currently excluded, include it (remove it from the exclusion list)
            cache$excluded_point_keys <- setdiff(cache$excluded_point_keys, clicked_key)
        } else {
            # currently included, exclude it (add it to the exclusion list)
            cache$excluded_point_keys <- append(clicked_key, cache$excluded_point_keys)
        }
    })
    # Derived Reactive: dCq ----------------------------------------------------

    dCq_data <- reactive({
        req(cq_data())
        req(nrow(cq_data()) > 0)
        req(length(input$hk_genes) > 0)

        HK_summary <- cq_data() |>
            filter(
                Keep,
                Target %in% input$hk_genes
            ) |>
            group_by(across(
                c("Sample", "Target", any_of("Replicate"))
            )) |>
            # summarize each HK separately
            summarize(
                HK_mean = mean(Cq, na.rm = TRUE),
                HK_sd   = sd(Cq, na.rm = TRUE),
                HK_n    = length(na.omit(Cq)), # sample size for each gene
                .groups = "drop"
            ) |>
            # summarize all HKs per sample
            group_by(across(
                c("Sample", any_of("Replicate"))
            )) |>
            summarize(
                HK_mean = mean(HK_mean, na.rm = TRUE),
                n_HK_genes = sum(HK_n > 0),
                # pooled SD for independet samples, allowing different mean (same as in ANOVA)
                HK_sd_pool = sqrt(
                    sum(HK_sd^2 * (HK_n - 1), na.rm = T) /
                        (sum(HK_n, na.rm = T) - n_HK_genes)
                ),
                HK_se_pooled = HK_sd_pool * sqrt(1 / (sum(HK_n, na.rm = T))),
                .groups = "drop"
            )

        cq_data() |>
            filter(
                Keep,
                !Target %in% input$hk_genes
            ) |>
            left_join(HK_summary) |>
            mutate(
                dCq     = Cq - HK_mean,
                exp_dCq = 2^-dCq
            )
    })
    # Derived Reactive: dCq summary per individual replicate -------------------

    dCq_rep_summary <- reactive({
        req(dCq_data())
        req(nrow(dCq_data()) > 0)

        dCq_data() |>
            group_by(across(
                c("Sample", "Target", any_of("Replicate"))
            )) |>
            summarize(
                Cq_n    = n_valid_Cq(Cq),
                Cq_sd   = sd_handle_inf(Cq),
                Cq_se   = Cq_sd / sqrt(Cq_n),
                dCq_mean = mean_handle_inf(dCq),
                # propagate SD and SE including HK variance/uncertainty.
                dCq_sd = ifelse(input$propagate_var,
                    # propagate HK SD
                    sqrt(Cq_sd^2 + HK_sd_pool^2),
                    # Compute Cq stats without propagating HK variance/uncertainty.
                    # This is statistically inaccurate but common in practice because it’s straightforward to compute.
                    Cq_sd # the same as sd of deltaCq
                ),
                dCq_se = ifelse(input$propagate_var,
                    # propagate HK SE
                    sqrt(Cq_se^2 + HK_se_pooled^2),
                    # Compute Cq stats without propagating HK variance/uncertainty.
                    # This is statistically inaccurate but common in practice because it’s straightforward to compute.
                    dCq_sd / sqrt(Cq_n)
                ),
                dCq_sd_low  = dCq_mean - dCq_sd,
                dCq_sd_high = dCq_mean + dCq_sd,
                dCq_se_low  = dCq_mean - dCq_se,
                dCq_se_high = dCq_mean + dCq_se,
                .groups = "drop"
            ) |>
            # compute exponentiated values
            mutate(
                exp_dCq_mean    = 2^-(dCq_mean),
                exp_dCq_sd_low  = 2^-(dCq_mean + dCq_sd),
                exp_dCq_sd_high = 2^-(dCq_mean - dCq_sd),
                exp_dCq_se_low  = 2^-(dCq_mean + dCq_se),
                exp_dCq_se_high = 2^-(dCq_mean - dCq_se)
            ) |>
            # mark undetected
            mutate(Undetected = Cq_n == 0)
    })
    # Derived Reactive: number of biological replicates ------------------------

    n_bio_reps <- reactive({
        req(dCq_data())
        req(nrow(dCq_data()) > 0)

        if (input$include_replicates & "Replicate" %in% names(dCq_data())) {
            dCq_data() |>
                pull("Replicate") |>
                unique() |>
                length()
        } else {
            1
        }
    })

    # make it available to javascript
    output$n_bio_reps <- reactive({
        n_bio_reps()
    })
    outputOptions(output, "n_bio_reps", suspendWhenHidden = FALSE)
    # Derived Reactive: number of samples --------------------------------------

    n_samples <- reactive({
        req(dCq_data())
        dCq_data() |>
            pull("Sample") |>
            unique() |>
            length()
    })

    # make it available to javascript
    output$n_samples <- reactive({
        n_samples()
    })
    outputOptions(output, "n_samples", suspendWhenHidden = FALSE)
    # Observe Event: disable biological replicate summary if only one replicate ----------

    observeEvent(n_bio_reps(), {
        if (n_bio_reps() > 1) {
            updateRadioGroupButtons(
                session,
                "summarize_bio_reps",
                disabledChoices = character(0)
            )
        } else {
            updateRadioGroupButtons(
                session,
                "summarize_bio_reps",
                disabledChoices = c("aggregate")
            )
        }
    })
    # Observe Event: Update default error bar type when toggling biological replicate summary ----------

    observeEvent(input$summarize_bio_reps, {
        if (input$summarize_bio_reps == "aggregate") {
            updateRadioGroupButtons(
                session,
                "stat_type",
                selected = "se"
            )
        } else {
            updateRadioGroupButtons(
                session,
                "stat_type",
                selected = "none"
            )
        }
    })
    # Derived Reactive: dCq summary aggregating replicates if present ----------

    dCq_summary <- reactive({
        req(dCq_rep_summary())
        req(n_bio_reps() > 1)

        dCq_rep_summary() |>
            group_by(across(c("Sample", "Target"))) |>
            summarize(
                dCq_n  = n_valid_Cq(dCq_mean),
                dCq_sd = sd_handle_inf(dCq_mean),
                dCq_se = dCq_sd / sqrt(dCq_n),
                # mean of means
                dCq_mean    = mean_handle_inf(dCq_mean),
                dCq_sd_low  = dCq_mean - dCq_sd,
                dCq_sd_high = dCq_mean + dCq_sd,
                dCq_se_low  = dCq_mean - dCq_se,
                dCq_se_high = dCq_mean + dCq_se,
                .groups = "drop"
            ) |>
            # compute exponentiated values
            mutate(
                exp_dCq_mean    = 2^-(dCq_mean),
                exp_dCq_sd_low  = 2^-(dCq_mean + dCq_sd),
                exp_dCq_sd_high = 2^-(dCq_mean - dCq_sd),
                exp_dCq_se_low  = 2^-(dCq_mean + dCq_se),
                exp_dCq_se_high = 2^-(dCq_mean - dCq_se)
            ) |>
            # mark undetected
            mutate(Undetected = dCq_n == 0)
    })
    # derived reactive: average dCq for reference sample -----------------------

    reference_sample_dCq <- reactive({
        req(nrow(dCq_rep_summary()) > 0)

        samples_metadata <- hot_to_r(input$samples_tab)
        # at least 2 distinct samples
        req(samples_metadata$New_Label |> unique() |> length() >= 2)

        reference_sample <- samples_metadata |>
            first() |>
            pull("New_Label")

        dCq_rep_summary() |>
            filter(Sample == reference_sample) |>
            select(Target, any_of("Replicate"),
                ref_dCq_mean = dCq_mean,
                ref_dCq_sd = dCq_sd,
                ref_dCq_se = dCq_se
            )
    })
    # Observe: allow ddCq calculation only if reference sample for slected target is valid ----------

    observeEvent(input$select_out_target, {
        req(input$select_out_target)
        req(input$select_out_target)

        current_target_ref_dCq <- reference_sample_dCq() |>
            filter(Target == input$select_out_target) |>
            pull("ref_dCq_mean")

        if (!all(is.finite(current_target_ref_dCq))) {
            updateRadioGroupButtons(
                session,
                "out_metric",
                disabledChoices = c("ddCq", "exp_ddCq")
            )
        } else {
            updateRadioGroupButtons(
                session,
                "out_metric",
                disabledChoices = character(0)
            )
        }
    })

    # Reactive: ddCq data points -----------------------------------------------

    ddCq_data <- reactive({
        req(reference_sample_dCq())

        dCq_data() |>
            left_join(reference_sample_dCq()) |>
            mutate(
                # Handle Inf - Inf = NaN case 
                ddCq     = ifelse(is.infinite(dCq) & is.infinite(ref_dCq_mean),
                                  0,
                                  dCq - ref_dCq_mean),
                exp_ddCq = 2^-ddCq
            )
    })

    ddCq_rep_summary <- reactive({
        req(reference_sample_dCq())

        dCq_rep_summary() |>
            left_join(reference_sample_dCq()) |>
            mutate(
                # Handle Inf - Inf = NaN case   
                ddCq_mean     = ifelse(is.infinite(dCq_mean) & is.infinite(ref_dCq_mean),
                                      0,
                                      dCq_mean - ref_dCq_mean),
                exp_ddCq_mean = 2^-ddCq_mean,
                ddCq_sd = ifelse(input$propagate_var,
                    # propagate control SD
                    sqrt(dCq_sd^2 + ref_dCq_sd^2),
                    dCq_sd # the same as sd of ddCq
                ),
                ddCq_se = ifelse(input$propagate_var,
                    sqrt(dCq_se^2 + ref_dCq_se^2),
                    dCq_se
                ),
                ddCq_sd_low      = ddCq_mean - ddCq_sd,
                ddCq_sd_high     = ddCq_mean + ddCq_sd,
                ddCq_se_low      = ddCq_mean - ddCq_se,
                ddCq_se_high     = ddCq_mean + ddCq_se,
                exp_ddCq_mean    = 2^-(ddCq_mean),
                exp_ddCq_sd_low  = 2^-(ddCq_mean + ddCq_sd),
                exp_ddCq_sd_high = 2^-(ddCq_mean - ddCq_sd),
                exp_ddCq_se_low  = 2^-(ddCq_mean + ddCq_se),
                exp_ddCq_se_high = 2^-(ddCq_mean - ddCq_se)
            )
    })

    ddCq_summary <- reactive({
        req(ddCq_rep_summary())
        req(n_bio_reps() > 1)

        ddCq_rep_summary() |>
            group_by(across(c("Sample", "Target"))) |>
            summarize(
                ddCq_n    = n_valid_Cq(ddCq_mean),
                ddCq_sd   = sd_handle_inf(ddCq_mean),
                ddCq_se   = ddCq_sd / sqrt(ddCq_n),
                # mean of means
                ddCq_mean    = mean_handle_inf(ddCq_mean),
                ddCq_sd_low  = ddCq_mean - ddCq_sd,
                ddCq_sd_high = ddCq_mean + ddCq_sd,
                ddCq_se_low  = ddCq_mean - ddCq_se,
                ddCq_se_high = ddCq_mean + ddCq_se,
                .groups = "drop"
            ) |>
            # compute exponentiated values
            mutate(
                exp_ddCq_mean    = 2^-(ddCq_mean),
                exp_ddCq_sd_low  = 2^-(ddCq_mean + ddCq_sd),
                exp_ddCq_sd_high = 2^-(ddCq_mean - ddCq_sd),
                exp_ddCq_se_low  = 2^-(ddCq_mean + ddCq_se),
                exp_ddCq_se_high = 2^-(ddCq_mean - ddCq_se)
            ) |>
            # mark undetected
            mutate(Undetected = ddCq_n == 0)
    })
    # Output Flags: Conditional panel visibility for statistical tests ---------

    stats_ui_flags <- reactive({
        req(input$stats_test)

        # Determine all UI flags based on selected test
        flags <- switch(input$stats_test,
            # --- ANCOVA ---
            "ancova" = list(
                show_unequal_variance_toggle    = FALSE,
                show_multiple_comparison_type   = n_samples() > 2,
                show_multiple_comparison_adjust = FALSE,
                show_post_hoc_test              = TRUE
            ),

            # --- Mixed Effect Model ---
            "mixed_effect" = list(
                show_unequal_variance_toggle    = TRUE,
                show_multiple_comparison_type   = TRUE,
                show_multiple_comparison_adjust = FALSE,
                show_post_hoc_test              = TRUE
            ),

            # --- ANOVA ---
            "anova" = list(
                show_unequal_variance_toggle    = FALSE,
                show_multiple_comparison_type   = TRUE,
                show_multiple_comparison_adjust = FALSE,
                show_post_hoc_test              = TRUE
            ),

            # --- Kruskal-Wallis ---
            "kruskal" = list(
                show_unequal_variance_toggle    = FALSE,
                show_multiple_comparison_type   = TRUE,
                show_multiple_comparison_adjust = isTRUE(input$stats_comparison == "pairwise"),
                show_post_hoc_test              = TRUE
            ),

            # --- Pairwise t-test ---
            "pairwise_ttest" = list(
                show_unequal_variance_toggle    = TRUE,
                show_multiple_comparison_type   = TRUE,
                show_multiple_comparison_adjust = TRUE,
                show_post_hoc_test              = FALSE
            ),

            # --- Pairwise paired t-test ---
            "pairwise_paired_ttest" = list(
                show_unequal_variance_toggle    = FALSE,
                show_multiple_comparison_type   = TRUE,
                show_multiple_comparison_adjust = TRUE,
                show_post_hoc_test              = FALSE
            ),

            # --- Pairwise Wilcoxon (signed-rank) ---
            "pairwise_wilcoxon" = list(
                show_unequal_variance_toggle    = FALSE,
                show_multiple_comparison_type   = TRUE,
                show_multiple_comparison_adjust = TRUE,
                show_post_hoc_test              = FALSE
            ),

            # --- Pairwise Wilcoxon-Mann-Whitney ---
            "pairwise_mann_whitney" = list(
                show_unequal_variance_toggle    = FALSE,
                show_multiple_comparison_type   = TRUE,
                show_multiple_comparison_adjust = TRUE,
                show_post_hoc_test              = FALSE
            ),

            # --- Default fallback ---
            list(
                show_unequal_variance_toggle    = FALSE,
                show_multiple_comparison_type   = FALSE,
                show_multiple_comparison_adjust = FALSE,
                show_post_hoc_test              = FALSE
            )
        )
        flags
    })

    # Expose individual flags as outputs for JavaScript conditional panels
    output$show_unequal_variance_toggle <- reactive({
        stats_ui_flags()$show_unequal_variance_toggle
    })
    outputOptions(output, "show_unequal_variance_toggle", suspendWhenHidden = FALSE)

    output$show_multiple_comparison_type <- reactive({
        stats_ui_flags()$show_multiple_comparison_type
    })
    outputOptions(output, "show_multiple_comparison_type", suspendWhenHidden = FALSE)

    output$show_multiple_comparison_adjust <- reactive({
        stats_ui_flags()$show_multiple_comparison_adjust
    })
    outputOptions(output, "show_multiple_comparison_adjust", suspendWhenHidden = FALSE)

    output$show_post_hoc_test <- reactive({
        stats_ui_flags()$show_post_hoc_test
    })
    outputOptions(output, "show_post_hoc_test", suspendWhenHidden = FALSE)
    # Observer: Update stats_test choices based on metric and sample count ----------

    observeEvent(list(input$stats_metric, n_samples(), n_bio_reps()), {
        req(input$stats_metric)
        req(n_samples() >= 2)

        metric <- input$stats_metric
        # Non-parametric tests require at least 5 biological replicates
        include_nonparam <- n_bio_reps() >= 5

        # Determine available tests and default based on metric and sample count
        if (n_samples() > 2) {
            # > 2 samples
            if (metric == "dCq") {
                # dCq: comparison across samples accounting for HK variance
                choices <- list(
                    "Parametric" = c(
                        "ANCOVA" = "ancova",
                        "Mixed Effect Model" = "mixed_effect",
                        "Pairwise paired t-test" = "pairwise_paired_ttest"
                    )
                )
                if (include_nonparam) {
                    choices[["Non-parametric"]] <- c("Pairwise Wilcoxon signed-rank (paired)" = "pairwise_wilcoxon")
                }
                default <- "ancova"
            } else {
                # ddCq or exp_ddCq: standard group comparisons
                choices <- list(
                    "Parametric" = c(
                        "ANOVA" = "anova",
                        "Pairwise t-test" = "pairwise_ttest"
                    )
                )
                if (include_nonparam) {
                    choices[["Non-parametric"]] <- c(
                        "Kruskal-Wallis" = "kruskal",
                        "Pairwise Wilcoxon-Mann-Whitney" = "pairwise_mann_whitney"
                    )
                }
                default <- "anova"
            }
        } else {
            # = 2 samples
            if (metric == "dCq") {
                choices <- list(
                    "Parametric" = c(
                        "ANCOVA" = "ancova",
                        "Paired t-test" = "pairwise_paired_ttest"
                    )
                )
                if (include_nonparam) {
                    choices[["Non-parametric"]] <- c("Wilcoxon signed-rank (paired)" = "pairwise_wilcoxon")
                }
                default <- "ancova"
            } else {
                # ddCq or exp_ddCq
                choices <- list(
                    "Parametric" = c("t-test" = "pairwise_ttest")
                )
                if (include_nonparam) {
                    choices[["Non-parametric"]] <- c("Wilcoxon-Mann-Whitney" = "pairwise_mann_whitney")
                }
                default <- "pairwise_ttest"
            }
        }

        # Freeze the reactive value to prevent stats_result from running with stale test value
        freezeReactiveValue(input, "stats_test")
        updatePickerInput(session, "stats_test", choices = choices, selected = default)
    })
    # Output: Post-hoc test description based on omnibus test and comparison ----------

    output$stats_posthoc <- renderText({
        req(input$stats_test)

        test <- input$stats_test
        comparison <- input$stats_comparison
        handle_variance <- isTRUE(input$stats_unequal_variance)

        # Determine post-hoc test based on omnibus and comparison type
        posthoc <- switch(test,
            "ancova" = if (comparison == "pairwise") "Tukey HSD" else "Dunnett",
            "mixed_effect" = if (handle_variance) {
                if (comparison == "pairwise") "Dunnett T3" else "Dunnett adjusted for unequal var."
            } else {
                if (comparison == "pairwise") "Tukey HSD" else "Dunnett"
            },
            "anova" = if (comparison == "pairwise") "Tukey HSD" else "Dunnett",
            "kruskal" = "Dunn's test",
            ""
        )

        posthoc
    })

    # Reactive: Run Statistical Test -------------------------------------------

    stats_result <- reactive({
        req(input$stats_test)
        req(input$select_out_target)
        req(n_bio_reps() >= 2)
        req(n_samples() >= 2)

        # statistical test settings
        test       <- input$stats_test
        response   <- input$stats_metric
        equal_var  <- !isTRUE(input$stats_unequal_variance)
        comparison <- input$stats_comparison
        p_adjust   <- input$stats_multiple_comparison_adjust



        # Prepare data based on metric
        if (response == "dCq") {
            # For dCq tests, we need dCq with reference sample info
            data <- dCq_rep_summary() |>
                filter(Target == input$select_out_target) |>
                left_join(reference_sample_dCq(), by = c("Target", "Replicate")) |>
                # Rename to match expected column names in statistical functions
                rename(dCq = dCq_mean, ref_dCq = ref_dCq_mean) |>
                drop_na(dCq, ref_dCq)

        } else {
            # For ddCq or exp_ddCq, use ddCq data
            data <- ddCq_rep_summary() |>
                filter(Target == input$select_out_target) |>
                # Rename to match expected column names in statistical functions
                rename(ddCq = ddCq_mean, exp_ddCq = exp_ddCq_mean) |>
                drop_na(ddCq, exp_ddCq)
        }
      
        # for ANCOVA and paired t-test drop entire replicate run
        if (test %in% c("ancova", "pairwise_paired_ttest")) {
            # Find replicates with any Inf/undetected values
            reps_with_inf <- data |>
                filter(!is.finite(.data[[response]])) |>
                pull(Replicate) |>
                unique()
            
            n_dropped_reps <- length(reps_with_inf)
            n_dropped_points <- 0  # no individual data points dropped
            
            # Drop entire replicates containing undetected values
            if (n_dropped_reps > 0) {
                data <- data |>
                    filter(!Replicate %in% reps_with_inf)
            }
        }
        # replace Inf with 999 for non-parametric tests (to avoid error in finding the max/min when more than 1 Inf is present)
        else if (test %in% c("kruskal", "pairwise_wilcoxon", "pairwise_mann_whitney")) {
            inf_replacemenet <- 999
            data <- data |>
                mutate(!!response := if_else(!is.finite(.data[[response]]), inf_replacemenet, .data[[response]]))
            
            n_dropped_points <- 0
            n_dropped_reps <- 0
        }
        # Filter out Inf values (undetected) for other parametric tests 
        else {
            n_before <- nrow(data)
            data_filtered <- data |>
                filter(is.finite(.data[[response]]))
            n_dropped_points <- n_before - nrow(data_filtered)
            n_dropped_reps <- 0
        
            # Use filtered data for analysis
            data <- data_filtered
        }

        # calculate samples after removing NA and eventual Inf for the selected target
        n_samples <- data$Sample |> unique() |> length()
        req(n_samples > 1)
        
        # Run appropriate test based on selection
        tryCatch({
            result <- switch(test,
                "ancova" = run_ancova(data, response = response, comparison = comparison),
                "mixed_effect" = run_mixed_effect(data, response = response, comparison = comparison, equal.var = equal_var),
                "anova" = run_anova(data, response = response,comparison = comparison),
                "kruskal" = run_kruskal(data, response = response, comparison = comparison, p.adjust = p_adjust),
                "pairwise_ttest" = run_pairwise_ttest(data, response = response, comparison = comparison, equal.var = equal_var, p.adjust = p_adjust),
                "pairwise_paired_ttest" = run_pairwise_paired_ttest(data, response = response, comparison = comparison, p.adjust = p_adjust),
                "pairwise_wilcoxon" = run_pairwise_wilcoxon(data, response = response, comparison = comparison, p.adjust = p_adjust),
                "pairwise_mann_whitney" = run_pairwise_mann_whitney(data, response = response, comparison = comparison, p.adjust = p_adjust)
            )
            # Add dropped count to result for warning display
            result$n_dropped_points <- n_dropped_points
            result$n_dropped_reps <- n_dropped_reps
            result
        }, error = function(e) {
            list(error = as.character(e$message), n_dropped_points = n_dropped_points, n_dropped_reps = n_dropped_reps)
            # print the error in the consol
            print(e)
        })
    })

    # Output: Stats card title with current target ------------------------------

    output$stats_card_title <- renderText({
        req(input$select_out_target)
        paste("Statistical Results for", input$select_out_target)
    })

    # Output: Flag for omnibus test type (for conditional UI) ------------------

    output$has_omnibus_test <- reactive({
        req(input$stats_test)
        input$stats_test %in% c("ancova", "mixed_effect", "anova", "kruskal")
    })
    outputOptions(output, "has_omnibus_test", suspendWhenHidden = FALSE)

    # Output: Dropped count for Inf values (for conditional warning) -----------
    output$stats_dropped_count <- reactive({
        req(stats_result())
        max(stats_result()$n_dropped_points, stats_result()$n_dropped_reps)
    })
    outputOptions(output, "stats_dropped_count", suspendWhenHidden = FALSE)

    # Output: Warning text for dropped Inf values ------------------------------
    output$stats_dropped_warning <- renderText({
        req(stats_result())
        req(stats_result()$n_dropped_points > 0 || stats_result()$n_dropped_reps > 0)
        
        n_dropped_points <- stats_result()$n_dropped_points
        n_dropped_reps <- stats_result()$n_dropped_reps
        
        if (n_dropped_reps > 0) {
            glue("{n_dropped_reps} replicate run{ifelse(n_dropped_reps > 1, 's', '')} with undetected values excluded from the analysis.")
        } else {
            glue("{n_dropped_points} data point{ifelse(n_dropped_points > 1, 's', '')} with undetected values (Inf) excluded from the analysis.")
        }
    })

    # Output: Omnibus badge (brief p-value indicator) --------------------------

    output$stats_omnibus_badge <- renderUI({
        req(stats_result()$omnibus)
        req(stats_result()$omnibus_pvalue)
        
        p <- stats_result()$omnibus_pvalue
        badge_class <- ifelse(p <= 0.05, "badge bg-success", "badge bg-secondary")
        badge_text  <- ifelse(p < 0.001,
                                "p < 0.001",
                                paste0("p = ", signif(p, digits = 2))
                            )
        
        div(class = badge_class, badge_text)
    })

    # Output: Omnibus test label -----------------------------------------------

    output$stats_omnibus_label <- renderText({
        req(stats_result()$omnibus)
        req(stats_result()$omnibus_label)
        stats_result()$omnibus_label |>
        # strip p-value that will be added in a label
            str_remove(", p = [0-9\\-e\\.]+$")
            
    })

    # Output: Omnibus test table -----------------------------------------------

    output$stats_omnibus_table <- DT::renderDataTable({
        req(stats_result())
        req(stats_result()$omnibus)

        stats_result()$omnibus |>
            DT::datatable(
                # remove additional elements such as paging, search etc.
                options = list(
                    layout = list(
                        topStart = NULL,
                        topEnd = NULL,
                        bottomStart = NULL,
                        bottomEnd = NULL
                    ),
                    paging = FALSE,
                    searching = FALSE,
                    ordering = FALSE,
                    info = FALSE
                ),
                rownames = FALSE,
                selection = "none",
                class = "compact stripe"
            ) |>
            DT::formatSignif(
                columns = intersect(
                    names(stats_result()$omnibus),
                    c("Sum Sq", "Mean Sq", "Df", "Numerator Df", "Denominator Df", "F-value", "Chi-sq", "p-value")
                ),
                digits = 2
            )
    })

    # Output: Comparison section title -----------------------------------------

    output$stats_comparison_title <- renderText({
        req(stats_result())
        
        if (!is.null(stats_result()$post_hoc_method)) {
            paste("Post-hoc:", stats_result()$post_hoc_method)
        } else if (!is.null(stats_result()$test_label)) {
            stats_result()$test_label
        } else {
            warning("Unknown Test")
            "Warning: Unknown Test"
        }
    })

    # Output: Results table (post-hoc or pairwise) -----------------------------

    output$stats_results_table <- DT::renderDataTable({
        req(stats_result())

        # Get the appropriate results table
        results_df <- if (!is.null(stats_result()$post_hoc)) {
            stats_result()$post_hoc
        } else if (!is.null(stats_result()$test_res)) {
            stats_result()$test_res
        } else {
            return(NULL)
        }

        # Identify numeric columns for rounding
        numeric_cols <- results_df |> select_if(is.numeric) |> colnames()
        
        results_df |>
            DT::datatable(
                options = list(
                    layout = list(
                        topStart = NULL,
                        topEnd = NULL,
                        bottomStart = NULL,
                        bottomEnd = NULL
                    ),
                    paging = FALSE,
                    searching = FALSE,
                    ordering = TRUE,
                    info = FALSE
                ),
                rownames = FALSE,
                selection = "none",
                class = "compact stripe"
            ) |>
            DT::formatRound(
                columns = numeric_cols,
                digits = 4
            ) |>
            DT::formatStyle(
                columns = "Significance",
                color = DT::styleEqual(
                    c("***", "**", "*", "ns"),
                    c("#198754", "#28a745", "#5cb85c", "#6c757d")
                ),
                fontWeight = "bold"
            )
    })

    # Output: Method description -----------------------------------------------

    output$stats_method <- renderText({
        req(stats_result())
        req(stats_result()$method)
        stats_result()$method
    })

    # Output: Results Plot -----------------------------------------------------
    observeEvent(dCq_data(), {
        req(dCq_data(), nrow(dCq_data()) > 0)

        non_hk_genes <- dCq_data()$Target |>
            unique() |>
            drop_empty()
        
        # Preserve current selection if it's still a valid (non-HK) gene
        current <- input$select_out_target
        selected <- if (!is.null(current) && current %in% non_hk_genes) current else non_hk_genes[1]
        
        updateSelectInput(session, "select_out_target", choices = non_hk_genes, selected = selected)
    })


    output$res_plot_title <- renderText({
        req(input$select_out_target)
        req(input$out_metric)
        y_label <- case_match(
            input$out_metric,
            "dCq" ~ "-ΔCq",
            "exp_dCq" ~ "2^-ΔCq",
            "ddCq" ~ "ΔΔCq",
            "exp_ddCq" ~ "2^-ΔΔCq",
        )

        paste(y_label, " Values for", input$select_out_target)
    })

    output$res_plot <- renderPlotly({
        req(input$select_out_target)
        req(input$stat_type)
        req(input$out_metric)
        req(nrow(dCq_data()) > 0)
        req(nrow(dCq_rep_summary()) > 0)

        if (input$summarize_bio_reps == "split" & input$out_metric %in% c("dCq", "exp_dCq")) {
            df_target <- dCq_data() |>
                filter(Target == input$select_out_target)

            y_value <- input$out_metric

            df_summary_target <- dCq_rep_summary() |>
                filter(Target == input$select_out_target)
        } else if (input$summarize_bio_reps == "aggregate" & input$out_metric %in% c("dCq", "exp_dCq")) {
            req(n_bio_reps() > 1)
            req(nrow(dCq_summary()) > 0)

            df_target <- dCq_rep_summary() |>
                filter(Target == input$select_out_target)

            y_value <- glue("{input$out_metric}_mean")

            df_summary_target <- dCq_summary() |>
                filter(Target == input$select_out_target)
        } else if (input$summarize_bio_reps == "split" & input$out_metric %in% c("ddCq", "exp_ddCq")) {
            df_target <- ddCq_data() |>
                filter(Target == input$select_out_target)

            y_value <- input$out_metric

            df_summary_target <- ddCq_rep_summary() |>
                filter(Target == input$select_out_target)
        } else if (input$summarize_bio_reps == "aggregate" & input$out_metric %in% c("ddCq", "exp_ddCq")) {
            req(n_bio_reps() > 1)
            req(nrow(ddCq_summary()) > 0)

            df_target <- ddCq_rep_summary() |>
                filter(Target == input$select_out_target)

            y_value <- glue("{input$out_metric}_mean")

            df_summary_target <- ddCq_summary() |>
                filter(Target == input$select_out_target)
        }


        y_label <- case_match(
            input$out_metric,
            "dCq" ~ "-ΔCq",
            "exp_dCq" ~ "2^-ΔCq",
            "ddCq" ~ "-ΔΔCq",
            "exp_ddCq" ~ "2^-ΔΔCq",
        )

        # used to invert sign when plotting -dCq
        if (input$out_metric %in% c("dCq", "ddCq")) {
            sign <- -1
        } else {
            sign <- 1
        }

        y_summary_value <- glue("{input$out_metric}_mean")

        if (input$stat_type == "none") {
            error_bar_high  <- NULL
            error_bar_low   <- NULL
            error_bar_label <- ""
        } else {
            error_bar_high  <- glue("{input$out_metric}_{input$stat_type}_high")
            error_bar_low   <- glue("{input$out_metric}_{input$stat_type}_low")
            error_bar_label <- glue("{str_to_upper(input$stat_type)}: ({round(df_summary_target[[error_bar_low]], 2)}; {round(df_summary_target[[error_bar_high]], 2)})")
        }

        # compute plot y limits -----------------------------------------------
        values <- df_target[[y_value]]

        # append Error bars values to compute y limits
        if (input$stat_type != "none") {
            values <- values |>
                append(c(
                    df_summary_target |> pull(error_bar_low),
                    df_summary_target |> pull(error_bar_high)
                ))
        }

        y_limits <- get_y_limits(values, metric = input$out_metric)

        undetected_present <- any(!is.finite(values[!is.na(values)]))
        if (input$out_metric == "dCq" & undetected_present) {
            y_min_label <- glue("≤{y_limits[1]}")
        } else {
            y_min_label <- glue("{y_limits[1]}")
        }

        # plot -----------------------------------------------------------------
        # add text label for hover label
        has_replicate <- "Replicate" %in% names(df_target)
        
        df_target <- df_target |>
            mutate(point_type_label = ifelse(Undetected, "Undetected", "Detected")) |>
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


        # add reference sample line
        if (input$out_metric == "exp_ddCq") {
            p <- ggplot() +
                geom_hline(yintercept = 1, linetype = "dashed", color = "gray30")
        } else if (input$out_metric == "ddCq") {
            p <- ggplot() +
                geom_hline(yintercept = 0, linetype = "dashed", color = "gray30")
        } else {
            p <- ggplot()
        }

        # don't plot average bar for dCq value (0 has no meaning)
        if (input$out_metric != "dCq") {
            p <- p +
                geom_bar(
                    data = df_summary_target,
                    aes(
                        x = Sample, y = sign * .data[[y_summary_value]],
                        text = text
                    ),
                    stat = "identity",
                    fill = accent_color(),
                    alpha = 0.5,
                    width = 0.6
                )
        }

        p <- p +
            geom_beeswarm(
                data = df_target,
                aes(
                    x = Sample, y = sign * .data[[y_value]],
                    text = text,
                    color = point_type_label,
                    shape = point_type_label
                    # label on hover
                ),
                method = "compactswarm", preserve.data.axis = TRUE
            )

        if (input$out_metric == "dCq") {
            p <- p +
                # add mean points
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
                values = c("Detected" = secondary_color(), "Undetected" = "#C03A2B", "Mean" = accent_color()),
                name = "",
            ) +
            labs(
                x = "Sample",
                y = y_label,
                title = NULL
            ) +
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

        # add error bars if requested
        if (input$stat_type != "none") {
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
                    color = accent_color()
                )
        }

        # facet by replicate if present
        if (input$summarize_bio_reps == "split" & ("Replicate" %in% names(df_target))) {
            p <- p + facet_wrap(~Replicate)
        }

        ggplotly(p, tooltip = "text") |>
            fix_plotly_legend()
    })
}

# Run App ======================================================================

shinyApp(ui = ui, server = server)
# TODO: remove hover on rectangle
