library(shiny)
library(bslib)
library(shinyWidgets)
library(bsicons)
library(rhandsontable)
library(colourpicker)
library(tidyverse)
library(plotly)
library(ggbeeswarm)
library(ggsignif)
library(ggh4x)
library(scales)
library(glue)
library(DT)
library(openxlsx)

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
source("R/censoring.R")
source("R/get_y_limits.R")
source("R/fix_plotly_legend.R")
source("R/is_HK.R")
source("R/missing_value_stats.R")
source("R/statistical_tests.R")
source("R/adjust_signif_position.R")
source("R/build_result_plots.R")

# UI Definition ================================================================

ui <- page_fillable(

    # Custom CSS
    tags$head(
        includeCSS("www/custom.css")
    ),
    
    # Main navigation
    navset_tab(
        id = "main_tabs",
        nav_item(h5("qPCR Analysis Tool", class = "ms-2")),
        nav_spacer(),
        # Panel 1: Input Data --------------------------------------------------
        nav_panel(
            title = "Input Data",
            page_sidebar(
                fillable = TRUE,
                sidebar = sidebar(
                    title = "Data Controls",
                    open = TRUE,
                    width = "380px",
                    prettySwitch(
                        "include_replicates",
                        label = "Include biological replicates column",
                        fill = TRUE, status = "primary",
                        value = FALSE
                    ),
                    hr(),
                    tags$h6(tags$strong("Samples")),
                    helpText("Edit 'New Label' to rename samples. Drag the row number to reorder them. Uncheck 'Include' to exclude samples from analysis."),
                    rHandsontableOutput("samples_tab"),
                    hr(),
                    tags$h6(tags$strong("Targets")),
                    helpText("Edit 'New Label' to rename targets. Drag the row number to reorder them. Uncheck 'Include' to exclude targets from analysis."),
                    rHandsontableOutput("targets_tab"),
                    hr(),
                    div(
                        class = "d-flex align-items-end gap-2",
                        numericInput(
                            "max_cycle",
                            label = "Undetected replacement cycle",
                            value = 40,
                            min = 1,
                            step = 1,
                            width = "220px"
                        ),
                        div(
                            class = "pb-2",
                            tooltip(
                                bs_icon("info-circle"),
                                tags$span(
                                    "Replacing non-detects with a maximum-cycle value can produce biased estimates.",
                                    tags$br(),
                                    tags$strong("Reference: "),
                                    "McCall, Matthew N et al. ‘On non-detects in qPCR data.’ ",
                                    tags$em("Bioinformatics (Oxford, England)"),
                                    " 30.16 (2014): 2310–2316. ",
                                    tags$a(
                                        "doi:10.1093/bioinformatics/btu239",
                                        href = "https://doi.org/10.1093/bioinformatics/btu239",
                                        target = "_blank",
                                        rel = "noopener noreferrer"
                                    )
                                ),
                                placement = "right"
                            )
                        )
                    ),
                    helpText("Undetected Cq values are replaced with this cycle.")
                ),
                
                # Main content area
                card(
                    max_height = "500px",
                    full_screen = TRUE,
                    fillable = TRUE,
                    card_header("qPCR Data Entry"),
                    helpText("Paste your qPCR data below directly from Excel"),
                    # height is actually max height
                    rHandsontableOutput("raw_data", height = "320px"),
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
                    selectInput(
                        "reference_sample",
                        "Reference sample (ΔΔCq / ANCOVA)",
                        choices = NULL
                    ),
                    prettySwitch(
                        inputId = "use_target_references",
                        label = "Different reference for each target",
                        value = FALSE,
                        fill = TRUE,
                        status = "primary"
                    ),
                    helpText(
                        "When enabled, select a target and then choose its reference. Turning this off makes the currently selected reference global and clears the per-target choices."
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
                        inputId   = "y_axis_mode",
                        label     = "Y Axis",
                        choices   = c("Auto" = "auto", "Free" = "free", "Custom" = "custom"),
                        selected  = "auto",
                        justified = TRUE,
                        width     = "60%",
                        size      = "sm"
                    ),
                    conditionalPanel(
                        condition = "input.y_axis_mode == 'custom'",
                        div(
                            class = "d-flex gap-2",
                            numericInput(
                                "y_axis_min",
                                label = "Y min",
                                value = NULL, step = 0.1,
                                width = "100px"
                            ),
                            numericInput(
                                "y_axis_max",
                                label = "Y max",
                                value = NULL, step = 0.1,
                                width = "100px"
                            )
                        )
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
                            div(class = "ms-2",
                                tooltip(
                                    bs_icon("exclamation-triangle"),
                                    "For technical replicates error bars are shown as QC only.
                                    They should not be used for statistical testing or to infer biological variability but rather to decide if keep or repeat the assay.
                                    For technical replicates, the preferred display is 'None' (individual points only) to avoid confusion with biological variability.
                                    If you want to add an error bar, SD is a more common choice because it reflects the variability of the assay; SEM reflects precision of the mean estimation (SD/√n) and shrinks with n.",
                                )
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
                            div(class = "ms-2",
                                tooltip(
                                    bs_icon("info-circle"),
                                    "Default: variance is computed from Target + HK technical replicates (HK pooled and propagated into ΔCq). For ΔΔCq, control variance is propagated to all samples (not ignored).
                                    This is more faithful to the measurement process.
                                    Stats without propagation of variance are offered for reproducibility with common practice. Note that this is not the correct approach though",
                                )
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

                        # Preferred analysis in Cq space.
                        conditionalPanel(
                            condition = "input.stats_metric == 'dCq'",
                            div(
                                class = "alert alert-success py-1 px-2 mb-2",
                                style = "font-size: 0.85em;",
                                bs_icon("lightbulb"),
                                "The preferred approach is to perform statistical analysis in Cq space.
                                ANCOVA and mixed-effect models account directly for differences between biological replicates, so a prior ΔΔCq transformation is not required (ANCOVA-derived sample effect sizes correspond to ΔΔCq. See ",
                                tooltip(
                                    tags$a(
                                        "Yuan et al. (2006)",
                                        href = "https://doi.org/10.1186/1471-2105-7-85",
                                        target = "_blank",
                                        rel = "noopener noreferrer"
                                    ),
                                    tags$span(
                                        tags$blockquote(
                                            class = "blockquote mb-2",
                                            style = "font-size: 0.9em;",
                                            "\"Since Ct is the observed value from experimental procedures, it should be the subject of statistical analysis.\""
                                        ),
                                        tags$blockquote(
                                            class = "blockquote mb-2",
                                            style = "font-size: 0.9em;",
                                            "\"An ANCOVA (analysis of covariance) model was proposed, and the ΔΔCt can be derived from analysis of effects of variables.\""
                                        ),
                                        tags$strong("Reference: "),
                                        "Yuan, J.S., et al. ‘Statistical analysis of real-time PCR data.’ ",
                                        tags$em("BMC Bioinformatics"),
                                        " 7, 85 (2006). ",
                                        tags$a(
                                            "doi:10.1186/1471-2105-7-85",
                                            href = "https://doi.org/10.1186/1471-2105-7-85",
                                            target = "_blank",
                                            rel = "noopener noreferrer"
                                        )
                                    ),
                                    placement = "right",
                                    options = list(customClass = "citation-tooltip")
                                ),
                                ").
                                ANCOVA is recommended when you have a clear reference/control sample (e.g., untreated sample).
                                Mixed-effect Model is better when the reference sample is arbitrary across replicates (e.g., comparing expression between different patients or cell lines)."
                            )
                        ),

                        # Warning for exp data
                        conditionalPanel(
                            condition = "input.stats_metric == 'exp_ddCq'",
                            div(
                                class = "alert alert-warning py-1 px-2 mb-2",
                                style = "font-size: 0.85em;",
                                bs_icon("exclamation-triangle"),
                                "Perform statistical analysis in Cq space (-ΔCq or -ΔΔCq), and use 2^-ΔΔCq to visualize linear fold changes. See ",
                                tooltip(
                                    tags$a(
                                        "Taylor et al. (2019)",
                                        href = "https://doi.org/10.1016/j.tibtech.2018.12.002",
                                        target = "_blank",
                                        rel = "noopener noreferrer"
                                    ),
                                    tags$span(
                                        tags$blockquote(
                                            class = "blockquote mb-2",
                                            style = "font-size: 0.9em;",
                                            "\"qPCR measurements are made on the log scale (Cq value) with statistical analysis performed in Cq space (i.e., ΔΔCq values or using log-transformed relative normalized expression [ΔCq]), while expression levels are reported after linear transformation of the ΔΔCq results.\""
                                        ),
                                        tags$strong("Reference: "),
                                        "Taylor, S.C., et al. ‘The Ultimate qPCR Experiment: Producing Publication Quality, Reproducible Data the First Time.’ ",
                                        tags$em("Trends in Biotechnology"),
                                        " 37(7), 761–774 (2019). ",
                                        tags$a(
                                            "doi:10.1016/j.tibtech.2018.12.002",
                                            href = "https://doi.org/10.1016/j.tibtech.2018.12.002",
                                            target = "_blank",
                                            rel = "noopener noreferrer"
                                        )
                                    ),
                                    placement = "right",
                                    options = list(customClass = "citation-tooltip")
                                ),
                                ".",
                            ),

                        ),

                        # Explain why pairwise Welch tests are the fallback for
                        # reference-normalized metrics.
                        conditionalPanel(
                            condition = "input.stats_metric == 'ddCq' || input.stats_metric == 'exp_ddCq'",
                            div(
                                class = "alert alert-info py-1 px-2 mb-2",
                                style = "font-size: 0.85em;",
                                bs_icon("info-circle"),
                                "When testing on -ΔΔCq or 2^-ΔΔCq, the reference sample has zero variance because all of its values are fixed at 0 or 1, respectively.
                                Welch's ANOVA cannot accommodate a zero-variance group.
                                Repeated Welch's t-tests can be used instead; each comparison against the reference reduces to a one-sample t-test against its fixed value."
                            )
                        ),
                        
                        # Omnibus test selection (choices update dynamically via pickerInput with optgroups)
                        div(
                            class = "d-flex justify-content-start align-items-center mb-2",
                            pickerInput(
                                inputId   = "stats_test",
                                label     = "Test:",
                                choices   = NULL, # update dinamically
                                width = "100%",
                                options = pickerOptions(container = "body")
                            ),

                            # Warning for Paired t-test with 2 samples
                            conditionalPanel(
                                condition = "input.stats_metric == 'dCq' && output.n_samples == 2 && input.stats_test == 'paired_ttest'",
                                div(class = "ms-2",
                                    tooltip(
                                        bs_icon("lightbulb"),
                                        "In standard conditions, the paired t-test yields identical results to the mixed-effect model (random intercept).
                                        Differences arise only in presence of missing data points (paired t-test uses only complete observations).",
                                        placement = "right"
                                    )
                                )
                            )
                        ),
                        
                        # Handle unequal variance for tests that support it.
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
                                # Explain the fixed-reference special case for independent t-tests.
                                conditionalPanel(
                                    condition = "input.stats_test == 'repeated_ttest' || input.stats_test == 'ttest'",
                                    tooltip(
                                        bs_icon("info-circle"),
                                        "Keeping this option enabled.
                                        Welch’s t-test is more robust when sample variances differ and remains appropriate for comparisons against the reference sample (zero variance).
                                        The option to disable it is provided for completeness but is generally not recommended.",
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
                                        "All vs Reference" = "trt.vs.ctrl"
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
                    uiOutput("failed_hk_sample_warning"),
                    plotlyOutput("res_plot", height = "100%")
                ),
                # Statistical Results Card (only shown when stats panel is active)
                conditionalPanel(
                    condition = "output.n_bio_reps >= 2 && output.n_samples >= 2",
                    card(
                        card_header(
                            textOutput("stats_card_title", inline = TRUE)
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
                                        span(class = "mx-3", style = "font-size: 16px;", textOutput("stats_omnibus_label", inline = TRUE)),
                                        uiOutput("stats_omnibus_badge"),
                                    ),
                                    value = "omnibus_panel",
                                    icon = NULL,
                                    DT::dataTableOutput("stats_omnibus_table"),
                                    conditionalPanel(
                                        condition = "output.stats_has_extra && output.stats_extra_in_omnibus",
                                        hr(),
                                        tags$h6(
                                            textOutput("stats_extra_title_omnibus"),
                                        ),
                                        tags$div(
                                            style = "max-width:450px; width:100%;",
                                            DT::dataTableOutput("stats_extra_table_omnibus", )
                                        )
                                    )
                                )
                            )
                        ),
                        
                        # Post-hoc or Pairwise comparison results section
                        div(
                            tags$h6(
                                class = "mb-2",
                                textOutput("stats_comparison_title", inline = TRUE)
                            ),
                            DT::dataTableOutput("stats_results_table"),
                            conditionalPanel(
                                condition = "output.stats_has_extra && !output.stats_extra_in_omnibus",
                                hr(),
                                tags$h6(
                                    textOutput("stats_extra_title_comparison"),
                                ),
                                tags$div(
                                    style = "max-width:450px; width:100%;",
                                    DT::dataTableOutput("stats_extra_table_comparison")
                                )
                            )

                        ),
                        uiOutput("all_undetected_comparison_warning"),
                        uiOutput("bio_rep_undetected_comparison_warning"),
                        # Method description
                        div(
                            class = "bg-light py-2 px-3",
                            tags$strong("Methods: "),
                            textOutput("stats_method", inline = TRUE)
                        ),
                    card_footer(
                        class = "d-flex justify-content-end",
                        downloadButton("download_stats_xlsx", "Export Stats (XLSX)", class = "btn-sm btn-outline-success")
                    )
                )
                )
            )
        ),
        # Panel 4: Export Tab --------------------------------------------------
        nav_panel(
            title = "Export",
            page_sidebar(
                fillable = TRUE,
                sidebar = sidebar(
                    title = "Plot Export Settings",
                    open = TRUE,
                    width = "380px",
                    
                    tags$h6(tags$strong("Plot Styling")),
                    # Linewidth slider (0 to 1 pt)
                    sliderInput(
                        "export_linewidth",
                        label = "Line width (pt)",
                        min = 0, max = 1, value = 0.5, step = 0.1
                    ),
                    
                    # Column width slider
                    sliderInput(
                        "export_bar_width",
                        label = "Column width",
                        min = 0.3, max = 0.9, value = 0.6, step = 0.05
                    ),
                    
                    # Point size slider
                    sliderInput(
                        "export_point_size",
                        label = "Point size",
                        min = 1, max = 5, value = 2, step = 0.5
                    ),
                    
                    # Sample colors (dynamic - rendered by server) in accordion
                    accordion(
                        id = "sample_colors_accordion",
                        open = FALSE,
                        accordion_panel(
                            title = "Sample Colors",
                            icon = icon("palette"),
                            uiOutput("sample_color_inputs")
                        )
                    ),
                    
                    # Axis text size slider (5 to 14 pt)
                    sliderInput(
                        "export_axis_text_size",
                        label = "Axis text size (pt)",
                        min = 5, max = 14, value = 10, step = 1
                    ),
                    hr(),
                    # Plot dimensions
                    tags$h6(tags$strong("Plot Dimensions")),
                    div(
                        class = "d-flex gap-2",
                        numericInput(
                            "export_plot_width",
                            label = "Width (cm)",
                            value = 4, min = 1, max = 30, step = 0.1,
                            width = "100px"
                        ),
                        numericInput(
                            "export_plot_height",
                            label = "Height (cm)",
                            value = 4, min = 2, max = 30, step = 0.1,
                            width = "100px"
                        )
                    ),
                    hr(),
                    # Significance bar display options
                    conditionalPanel(
                        condition = "output.n_bio_reps >= 2 && output.n_samples >= 2 && input.summarize_bio_reps != 'split'",
                        tags$h6(tags$strong("Significance Bars")),
                        prettySwitch(
                            "show_signif_bars",
                            label = "Show significance bars",
                            fill = TRUE, status = "primary",
                            value = FALSE
                        ),
                        conditionalPanel(
                            condition = "input.show_signif_bars",
                            prettySwitch(
                                "hide_ns_bars",
                                label = "Hide non-significant (ns)",
                                fill = TRUE, status = "primary",
                                value = TRUE
                            ),
                            prettySwitch(
                                "show_exact_pvalue",
                                label = "Show exact p-values",
                                fill = TRUE, status = "primary",
                                value = FALSE
                            ),
                            sliderInput(
                                "export_signif_text_size",
                                label = "Significance text size (pt)",
                                min = 5, max = 14, value = 8, step = 1
                            )
                        )
                    ),

                    
                    # Download buttons moved to respective cards
                ),
                
                # Main content area
                layout_columns(
                    col_widths = c(12),
                    row_heights = c("1fr", "1fr"),
                    
                    # Static plot preview
                    card(
                        full_screen = TRUE,
                        card_header(uiOutput("plot_preview_title")),
                        plotOutput("export_plot", height = "100%"),
                        card_footer(
                            class = "d-flex justify-content-end gap-2",
                            downloadButton("download_plot_png", "Plot (PNG)", class = "btn-sm btn-outline-primary"),
                            downloadButton("download_plot_pdf", "Plot (PDF)", class = "btn-sm btn-outline-primary")
                        )
                    ),
                    
                    # Data preview tabs
                    navset_card_tab(
                        id = "data_preview_tabs",
                        full_screen = TRUE,
                        nav_panel(
                            title = "Raw Cq",
                            DT::dataTableOutput("preview_raw_cq")
                        ),
                        nav_panel(
                            title = "Technical Replicates",
                            DT::dataTableOutput("preview_technical")
                        ),
                        nav_panel(
                            title = "Bio Rep Averages",
                            DT::dataTableOutput("preview_bio_rep")
                        ),
                        nav_panel(
                            title = "Summary",
                            conditionalPanel(
                                condition = "output.n_bio_reps >= 2",
                                DT::dataTableOutput("preview_summary")
                            )
                        ),
                        
                        nav_spacer(),
                        
                        nav_item(
                            tags$div( # required to escame the formatting inherited by the nav_item
                                downloadButton("download_data_xlsx", "Data (XLSX)", class = "btn-sm btn-outline-success border-0")
                            )
                        ),
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
        # Server-side backing data used to initialize or replace the raw table.
        # Live user edits remain in input$raw_data until a server update is needed.
        raw_data = empty_raw_data,

        # sample control table (rename, reorder, exclude)
        # samples rename, reordering and exclusion should be retrieved from input$samples_tab
        samples_tab = data.frame(
            Sample    = character(),
            New_Label = character(),
            Include   = logical(),
            stringsAsFactors = FALSE
        ),

        # target control table (rename, reorder, exclude)
        # targets rename, reordering and exclusion should be retrieved from input$targets_tab
        targets_tab = data.frame(
            Target    = character(),
            New_Label = character(),
            Include   = logical(),
            stringsAsFactors = FALSE
        ),
        
        # list of excluded points
        excluded_point_keys = c(),
        selected_ct_target  = c(),
        targets_available   = c(),
        max_cycle           = 40,
        # The global reference is reused by every target unless per-target mode is enabled.
        global_reference_sample = NULL,
        # Named character vector: target names map to their explicitly selected references.
        target_reference_samples = character(),
        # Remember the last invalid reference context to avoid reopening the same modal.
        reference_prompt_key = NULL
    )

    # Prefer the live browser state, falling back to the cached raw data
    # while the widget is initializing.
    current_raw_data <- reactive({
        if (!is.null(input$raw_data)) {
            return(hot_to_r(input$raw_data))
        }

        cache$raw_data
    })

    # If pasted data exceed the configured max cycle, raise the replacement to the
    # first integer above the largest detected Cq.
    max_cycle_value <- reactive({
        requested <- validate_max_cycle(input$max_cycle) %||%
            validate_max_cycle(cache$max_cycle) %||%
            40

        raw_data <- current_raw_data()

        if (is.null(raw_data) || !"Cq" %in% names(raw_data)) {
            return(requested)
        }

        parsed_cq <- parse_Cq(raw_data$Cq)

        next_integer_above(
            parsed_cq$Cq,
            minimum = requested,
            censored = parsed_cq$Cq_censored
        )
    })

    observeEvent(list(input$max_cycle, input$raw_data), {
        replacement <- max_cycle_value()
        old_value <- suppressWarnings(as.numeric(input$max_cycle))
        cache$max_cycle <- replacement

        # Update the max_cycle input if the input is invalid or raw_data exceeds the requested max cycle.
        if (length(old_value) == 0 || is.na(old_value) || !isTRUE(all.equal(old_value, replacement))) {
            updateNumericInput(session, "max_cycle", value = replacement)
            if (length(old_value) == 1 && is.finite(old_value) && replacement > old_value) {
                showNotification(
                    glue("Undetected replacement increased to {replacement}, the first integer above the largest detected Cq."),
                    type = "message",
                    duration = 6
                )
            }
        }

        # Store every censored value using the normalized ">cycle" display form
        # so the table and downstream exports show the same threshold.
        if (!is.null(input$raw_data)) {
            current_data <- current_raw_data()
            parsed_cq <- parse_Cq(current_data$Cq)
            censored_rows <- parsed_cq$Cq_censored
            censored_cq_value <- paste0(">", format_qpcr_number(replacement))

            # Update the backing data only when its displayed values changed.
            # This intentionally re-renders the table without doing so on every edit.
            if (any(censored_rows) && any(as.character(current_data$Cq[censored_rows]) != censored_cq_value)) {
                cache$raw_data <- current_data |>
                    mutate(
                        Cq = if_else(
                            censored_rows,
                            censored_cq_value,
                            as.character(Cq)
                        )
                    )
            }
        }
    }, ignoreInit = TRUE)
    # Observer: Toggle biological replicates column ----------------------------
    
    observeEvent(input$include_replicates, {
        if (is.null(input$raw_data)) {
            return()
        }
        
        current_data <- current_raw_data()
        
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
    outputOptions(output, "raw_data", suspendWhenHidden = FALSE)
    
    # Observer: on raw data edit: ----------------------------------------------
    # 1. validate Cq values and push conversions back to the cached `raw_data`
    # 2. update `cache$samples_tab` when samples change in raw_data while preserving previous edits
    
    observeEvent(input$raw_data, {
        # 1. validate Cq values ------------------------------------------------
        current_data <- current_raw_data()
        
        original_cq <- current_data$Cq |>
            as.character() |>
            str_trim() |>
            replace_na("") |>
            # ignore trailing 0 in decimals
            str_remove("(?<=[0-9]\\.[0-9]{0,10})0+$")
        
        parsed_cq <- parse_Cq(original_cq)
        replacement <- max_cycle_value()
        normalized_cq <- ifelse(
            parsed_cq$Cq_censored,
            paste0(">", format_qpcr_number(replacement)),
            ifelse(
                is.na(parsed_cq$Cq),
                "",
                format_qpcr_number(parsed_cq$Cq, digits = 10)
            )
        )

        # Find values that changed
        changed_mask <- original_cq != normalized_cq

        # Equivalent numeric spellings are normalized silently. Only surface
        # semantic conversions (for example, "Undetermined" to ">40").
        meaningful_change <- meaningful_cq_conversion(
            original = original_cq,
            parsed_numeric = parsed_cq$Cq,
            censored = parsed_cq$Cq_censored,
            normalized_cq = normalized_cq
        )

        conversions <- map2_chr(
            original_cq[meaningful_change], normalized_cq[meaningful_change],
            ~ paste0(.x, " → ", .y)
        ) |>
            unique()
        
        # Push normalized values back only when validation changed something;
        # updating the backing data on every keystroke would re-render the table.
        if (any(changed_mask)) {
            cache$raw_data <- current_data |>
                mutate(Cq = normalized_cq)
        }
        
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
        
        current_samples <- current_data$Sample |>
            unique() |>
            drop_empty()
        
        previous_samples <- cache$samples_tab$Sample
        # new_samples <- setdiff(current_samples, previous_samples)
        
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
        
        # 3. cache last state of targets_tab -----------------------------------
        if (!is.null(input$targets_tab)) {
            cache$targets_tab <- hot_to_r(input$targets_tab)
        }
        
        current_targets <- current_data$Target |>
            unique() |>
            drop_empty()
        
        previous_targets <- cache$targets_tab$Target
        
        if (length(current_targets) == 0) {
            cache$targets_tab <- data.frame(
                Target    = character(),
                New_Label = character(),
                Include   = logical()
            )
        } else if (length(previous_targets) == 0) {
            cache$targets_tab <- data.frame(
                Target    = current_targets,
                New_Label = current_targets,
                Include   = rep(TRUE, times = length(current_targets))
            )
        } else {
            cache$targets_tab <- data.frame(Target = current_targets) |>
                left_join(cache$targets_tab) |>
                mutate(
                    New_Label = coalesce(New_Label, Target),
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
    # Output: Target control table ---------------------------------------------
    output$targets_tab <- renderRHandsontable({
        req(nrow(cache$targets_tab) > 0)
        
        rhandsontable(
            cache$targets_tab,
            rowHeaders = TRUE,
            readOnly = FALSE,
            stretchH = "all",
            contextMenu = FALSE,
            manualRowMove = TRUE
        ) |>
            hot_col("Target", readOnly = TRUE) |>
            hot_col("New_Label", type = "text") |>
            hot_col("Include", type = "checkbox") |>
            hot_cols(columnSorting = FALSE)
    })
    # Derived Reactive: Processed data (with parsed Cq, sample/target renames, ordering and exclusions) ----
    
    cq_data <- reactive({
        raw_data <- current_raw_data()
        req(raw_data)
        req(nrow(hot_to_r(input$samples_tab)) > 0)
        req(!is.null(input$targets_tab), nrow(hot_to_r(input$targets_tab)) > 0)
        
        samples_metadata <- hot_to_r(input$samples_tab)
        targets_metadata <- hot_to_r(input$targets_tab)

        raw_data |>
            # Parse values and censoring status together, then apply the numeric replacement.
            # In Raw data, censored values are stored as ">{max_cycle}" strings, but in the processed data they are replaced with the numeric replacement value.
            # The ">{max_cycle}" strings are preserved in the `Cq_display` column for table display and export.
            parse_Cq_data() |>
            mutate(
                Cq = replace_censored(
                    Cq,
                    censored = Cq_censored,
                    replacement = max_cycle_value()
                ),
                Cq_display = format_censored_value(Cq, Cq_censored)
            ) |>
            mutate(Key = row_number()) |> # add unique Key ID matching raw data rows
            relocate(Key, .before = 1) |>
            # join with sample metadata for renaming, reordering and exclusion
            inner_join(samples_metadata, by = "Sample") |>
            filter(Include) |>
            # update and reorder sample name
            mutate(Sample = factor(New_Label,
                                   levels = unique(samples_metadata$New_Label)
            )) |>
            select(-New_Label, -Include) |>
            # join with target metadata for renaming, reordering and exclusion
            inner_join(targets_metadata, by = "Target") |>
            filter(Include) |>
            mutate(Target = factor(New_Label,
                                    levels = unique(targets_metadata$New_Label)
            )) |>
            select(-New_Label, -Include) |>
            arrange(Sample) |>
            # mark excluded points
            mutate(
                Keep = !Key %in% cache$excluded_point_keys
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
        if (length(cache$selected_ct_target) == 1 && cache$selected_ct_target %in% targets) {
            selected <- cache$selected_ct_target
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
        cache$selected_ct_target <- input$select_ct_target
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
                point_type_label = ifelse(Cq_censored, "Undetected", "Detected"),
            )
        
        n_samples <- df_target$Sample |>
            unique() |>
            length()
        
        df_summary_target <- df_target |>
            filter(Keep) |>
            retain_detected_or_all_censored() |>
            group_by(across(
                c("Sample", "Target", any_of("Replicate"))
            )) |>
            summarize(
                mean = mean(Cq),
                mean_censored = all(Cq_censored),
                mean_display = format_censored_value(mean, mean_censored),
                point_type_label = "Mean",
                Keep_label = NA # initialize empty keep_label to avoid duplication of the legend in ggplotly
            )
        
        # Reuse the same replacement cycle for limits, shading, and labels.
        replacement_cycle <- max_cycle_value()

        # force a minumum of y-axis range of 3 units
        y_limits <- get_Cq_y_limits(
            df_target$Cq,
            min_range = 3,
            undetected_value = replacement_cycle,
            undetected_present = any(df_target$Cq_censored, na.rm = TRUE)
        )
        
        # Pre-compute hover text
        has_replicate <- "Replicate" %in% names(df_target)
        
        df_target <- df_target |>
            mutate(text = if (has_replicate) {
                glue(
                    "{Sample} ({Replicate})
                    Target: {Target}
                    Cq: {Cq_display}
                    {ifelse(Keep, '', '(excluded)')}"
                )
            } else {
                glue(
                    "{Sample}
                    Target: {Target}
                    Cq: {Cq_display}
                    {ifelse(Keep, '', '(excluded)')}"
                )
            })
        
        df_summary_target <- df_summary_target |>
            mutate(text = if (has_replicate) {
                glue(
                    "{Sample} ({Replicate})
                    Target: {Target}
                    Mean Cq: {mean_display}"
                )
            } else {
                glue(
                    "{Sample}
                    Target: {Target}
                    Mean Cq: {mean_display}"
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
            annotate("rect", xmin = 0.5, xmax = n_samples + 0.5, ymin = replacement_cycle - 5, ymax = replacement_cycle, alpha = 0.6, fill = "#EBEBEB") +
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
                labels = function(x) ifelse(
                    x == replacement_cycle,
                    paste0(">", replacement_cycle),
                    x
                ),
                oob = scales::oob_keep
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
            suppress_empty_filled_hover() |>
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
    # Derived Reactive: housekeeping-gene status per biological sample --------

    hk_sample_status <- reactive({
        data <- cq_data()
        req(data, nrow(data) > 0, length(input$hk_genes) > 0)

        all_samples <- data |>
            filter(Keep) |>
            distinct(across(any_of(c("Sample", "Replicate"))))

        detected_hk_by_sample <- data |>
            filter(Keep, Target %in% input$hk_genes) |>
            group_by(across(any_of(c("Sample", "Target", "Replicate")))) |>
            summarize(
                HK_detected = any(!is.na(Cq) & !Cq_censored),
                .groups = "drop"
            ) |>
            group_by(across(any_of(c("Sample", "Replicate")))) |>
            summarize(
                n_detected_HK_genes = sum(HK_detected),
                .groups = "drop"
            )

        all_samples |>
            left_join(detected_hk_by_sample) |>
            mutate(
                n_detected_HK_genes = coalesce(n_detected_HK_genes, 0),
                HK_valid = n_detected_HK_genes == length(input$hk_genes)
            )
    })

    failed_hk_samples <- reactive({
        hk_sample_status() |>
            filter(!HK_valid)
    })

    output$failed_hk_sample_warning <- renderUI({
        failed_samples <- failed_hk_samples()
        req(nrow(failed_samples) > 0)

        sample_labels <- if ("Replicate" %in% names(failed_samples)) {
            paste0(failed_samples$Sample, " (", failed_samples$Replicate, ")")
        } else {
            as.character(failed_samples$Sample)
        }

        div(
            class = "alert alert-warning py-2 px-3 m-2 d-flex align-items-start gap-2",
            style = "font-size: 0.9em;",
            bs_icon("exclamation-triangle"),
            tags$span(glue(
                "Excluded sample(s) because at least one selected housekeeping gene had no detected Cq: {toString(sample_labels)}."
            ))
        )
    })

    # Derived Reactive: dCq ----------------------------------------------------
    
    dCq_data <- reactive({
        req(cq_data())
        req(nrow(cq_data()) > 0)
        req(length(input$hk_genes) > 0)

        valid_samples <- hk_sample_status() |>
            filter(HK_valid) |>
            select(any_of(c("Sample", "Replicate")))

        HK_per_gene <- cq_data() |>
            filter(
                Keep,
                Target %in% input$hk_genes
            ) |>
            semi_join(valid_samples) |>
            retain_detected_or_all_censored() |>
            group_by(across(c("Sample", "Target", any_of("Replicate")))) |>
            # summarize each HK separately
            summarize(
                HK_mean = mean(Cq),
                HK_sd   = sd(Cq),
                HK_n    = n(),
                .groups = "drop"
            )
        
        # Aggregate all HKs per sample
        HK_summary <- HK_per_gene |>
            group_by(across(c("Sample", any_of("Replicate")))) |>
            summarize(
                HK_mean = mean(HK_mean),
                n_HK_genes = n(),
                # pooled SD for independet samples, allowing different mean (same as in ANOVA)
                HK_sd_pool = sqrt(
                    sum(HK_sd^2 * (HK_n - 1), na.rm = T) /
                        (sum(HK_n, na.rm = T) - n_HK_genes)
                ),
                HK_se_pooled = HK_sd_pool * sqrt(1 / (sum(HK_n, na.rm = T))),
                .groups = "drop"
            )
        
        # Per-HK gene average columns (only if >1 HK gene)
        if (length(input$hk_genes) > 1) {
            HK_wide <- HK_per_gene |>
                select(c("Sample", any_of("Replicate"), "Target", "HK_mean")) |>
                pivot_wider(
                    names_from = Target,
                    values_from = HK_mean,
                    names_glue = "{.value}_{Target}_Cq"
                )
        }
        
        result <- cq_data() |>
            filter(
                Keep,
                !Target %in% input$hk_genes
            ) |>
            retain_detected_or_all_censored() |>
            inner_join(HK_summary) |>
            mutate(
                dCq     = Cq - HK_mean,
                dCq_censored = Cq_censored,
                exp_dCq = 2^-dCq
            )
        
        # Join per-HK gene columns if >1 HK gene
        if (length(input$hk_genes) > 1) {
            result <- result |> left_join(HK_wide)
        }
        
        result
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
                # Number of technical values retained for the calculation.
                Cq_n    = n(),
                Cq_mean = mean(Cq),
                Cq_censored = all(Cq_censored),
                Cq_sd   = sd(Cq),
                Cq_se   = Cq_sd / sqrt(Cq_n),
                HK_mean_Cq = mean(HK_mean),
                # carry forward individual HK gene averages (constant within group)
                across(starts_with("HK_mean_") & ends_with("_Cq"), mean),
                dCq_mean = mean(dCq),
                dCq_censored = all(dCq_censored),
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
            )
    })
    # Derived Reactive: number of biological replicates ------------------------
    
    n_bio_reps <- reactive({
        req(dCq_rep_summary())
        req(nrow(dCq_rep_summary()) > 0)
        req(input$select_out_target)
        
        if (input$include_replicates & "Replicate" %in% names(dCq_rep_summary())) {
            dCq_rep_summary() |>
                filter(!is.na(dCq_mean)) |>
                filter(Target == input$select_out_target) |>
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
        req(dCq_rep_summary())
        req(input$select_out_target) 
        
        dCq_rep_summary() |>
            filter(Target == input$select_out_target) |>
            pull("Sample") |>
            unique() |>
            length()
    })

    # Warn in the statistical-analysis section when replacement values are the
    # only target information available for more than one sample. Any apparent
    # difference between these samples comes from their HK values, not from a
    # measured difference in the target.
    all_undetected_samples <- reactive({
        req(dCq_rep_summary(), input$select_out_target)
        dCq_rep_summary() |>
            filter(Target == input$select_out_target) |>
            all_censored_groups(group_col = "Sample", censored_col = "Cq_censored")
    })

    # Pairwise comparisons in which both samples contain
    # fully undetected biological replicate require caution
    # even when their other replicates were detected.
    samples_with_undetected_bio_reps <- reactive({
        req(dCq_rep_summary(), input$select_out_target)

        dCq_rep_summary() |>
            filter(Target == input$select_out_target) |>
            group_by(Sample) |>
            summarize(
                has_undetected_bio_rep = any(Cq_censored, na.rm = TRUE),
                .groups = "drop"
            ) |>
            filter(has_undetected_bio_rep) |>
            pull(Sample) |>
            as.character()
    })

    undetected_bio_rep_comparison_pairs <- reactive({
        affected_samples <- samples_with_undetected_bio_reps()
        if (length(affected_samples) < 2) return(character())

        pairs <- combn(affected_samples, 2, simplify = FALSE)
        fully_undetected <- all_undetected_samples()

        # Pairs in which both samples are fully undetected are already covered
        # by the stronger warning below.
        pairs <- Filter(
            function(pair) !all(pair %in% fully_undetected),
            pairs
        )
        vapply(pairs, paste, collapse = " vs ", FUN.VALUE = character(1))
    })

    output$all_undetected_comparison_warning <- renderUI({
        samples <- all_undetected_samples()
        req(length(samples) > 1)

        div(
            class = "alert alert-warning py-2 px-3 m-2 d-flex align-items-start gap-2",
            style = "font-size: 0.9em;",
            bs_icon("exclamation-triangle"),
            tags$span(glue(
                "{toString(samples)} are undetected in every biological replicate for {input$select_out_target}.
                Comparisons among these samples should not be considered meaningful.
                Any apparent difference in normalized expression reflect HK differences rather than differences in the target itself:
                in a sample with more abundant RNA (lower HK Cq), the target would need to be very scarce to remain undetected.
                In a sample with less abundant RNA (higher HK Cq), an undetected result is relatively easier to obtain even when the target's normalized expression could still be appreciable."
            ))
        )
    })

    output$bio_rep_undetected_comparison_warning <- renderUI({
        comparison_pairs <- undetected_bio_rep_comparison_pairs()
        req(length(comparison_pairs) > 0)

        div(
            class = "alert alert-warning py-2 px-3 m-2 d-flex align-items-start gap-2",
            style = "font-size: 0.9em;",
            bs_icon("exclamation-triangle"),
            tags$span(glue(
                "The following comparison(s) involve two samples that both contain undetected values: {toString(comparison_pairs)}.
                Comparisons among these samples should be interpreted with caution.
                Apparent differences in normalized expression of undetected values reflect HK differences rather than differences in the target itself:
                in a sample with more abundant RNA (lower HK Cq), the target would need to be very scarce to remain undetected.
                In a sample with less abundant RNA (higher HK Cq), an undetected result is relatively easier to obtain even when the target's normalized expression could still be appreciable."
            ))
        )
    })

    # Samples with a non-missing numeric value for selecting the test family
    n_numeric_samples <- reactive({
        req(input$select_out_target) 
        req(input$stats_metric)

        if(input$stats_metric == "dCq") {
            req(dCq_rep_summary())

            dCq_rep_summary() |>
                filter(Target == input$select_out_target) |>
                filter(!is.na(dCq_mean)) |>
                pull("Sample") |>
                unique() |>
                length()
        } else {
            req(ddCq_rep_summary())

            stat_metric <- paste0(input$stats_metric, "_mean")
            ddCq_rep_summary() |>
                filter(Target == input$select_out_target) |>
                filter(!is.na(.data[[stat_metric]])) |>
                pull("Sample") |>
                unique() |>
                length()
        }
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
            # Force selection to 'split' if 'aggregate' was selected
            if (input$summarize_bio_reps == "aggregate") {
                updateRadioGroupButtons(
                    session,
                    "summarize_bio_reps",
                    selected = "split",
                    disabledChoices = c("aggregate")
                )
            } else {
                updateRadioGroupButtons(
                    session,
                    "summarize_bio_reps",
                    disabledChoices = c("aggregate")
                )
            }
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
            # Turn off significance bars in split mode
            shinyWidgets::updatePrettySwitch(
                session, "show_signif_bars", value = FALSE
            )
        }
    })
    # Observer: Update Y Axis mode choices based on aggregation mode -----------
    
    observeEvent(input$summarize_bio_reps, {
        if (input$summarize_bio_reps == "aggregate") {
            # Hide 'Free' and reset to 'auto' if 'free' was selected
            updateRadioGroupButtons(
                session,
                "y_axis_mode",
                disabledChoices = "free",
                selected = if (input$y_axis_mode == "free") "auto" else input$y_axis_mode
            )
        } else {
            updateRadioGroupButtons(
                session,
                "y_axis_mode",
                disabledChoices = character(0)
            )
        }
    })
    # Observer: Pre-fill custom Y axis inputs with auto values ----------------
    
    observeEvent(input$y_axis_mode, {
        if (input$y_axis_mode == "custom") {
            d <- tryCatch(result_plot_data(), error = function(e) NULL)
            if (!is.null(d)) {
                y_lim <- d$y_limits
                updateNumericInput(session, "y_axis_min", value = round(y_lim[1], 2))
                updateNumericInput(session, "y_axis_max", value = round(y_lim[2], 2))
            }
        }
    })
    # Derived Reactive: dCq summary aggregating replicates if present ----------
    
    dCq_summary <- reactive({
        req(dCq_rep_summary())
        req(n_bio_reps() > 1)
        
        dCq_rep_summary() |>
            group_by(across(c("Sample", "Target"))) |>
            summarize(
                dCq_n  = n(),
                dCq_undetected_n = sum(dCq_censored, na.rm = TRUE),
                dCq_sd = sd(dCq_mean),
                dCq_se = dCq_sd / sqrt(dCq_n),
                # mean of means
                dCq_mean    = mean(dCq_mean),
                dCq_censored = any(dCq_censored, na.rm = TRUE),
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
            )
    })
    # Reference sample selection and average dCq -------------------------------

    reference_sample_status <- reactive({
        req(dCq_rep_summary(), hk_sample_status(), input$select_out_target)

        # Start from every expected sample/replicate, including runs removed
        # from dCq_data() because their housekeeping gene was not detected.
        expected_runs <- hk_sample_status() |>
            transmute(
                Sample = as.character(Sample),
                across(any_of("Replicate")),
                HK_valid
            )

        observed_target_runs <- dCq_rep_summary() |>
            filter(as.character(Target) == input$select_out_target) |>
            transmute(
                Sample = as.character(Sample),
                across(any_of("Replicate")),
                dCq_mean,
                dCq_censored
            )

        run_keys <- intersect(
            c("Sample", "Replicate"),
            intersect(names(expected_runs), names(observed_target_runs))
        )

        expected_runs |>
            left_join(observed_target_runs, by = run_keys) |>
            group_by(Sample) |>
            summarize(
                has_undetected = any(dCq_censored %in% TRUE),
                has_failed_hk = any(!(HK_valid %in% TRUE)),
                has_missing_target = any(HK_valid %in% TRUE & is.na(dCq_mean)),
                reference_complete = !has_undetected &&
                    !has_failed_hk && !has_missing_target,
                .groups = "drop"
            )
    })

    # Leaving per-target mode promotes the reference currently visible in the
    # selector to the global reference and discards every target override.
    observeEvent(input$use_target_references, {
        req(input$reference_sample)

        if (!isTRUE(input$use_target_references)) {
            cache$global_reference_sample <- input$reference_sample
            cache$target_reference_samples <- character()
            cache$reference_prompt_key <- NULL
        }
    }, ignoreInit = TRUE, priority = 100)

    # Populate the reference selector and restore the cached global or
    # target-specific selection whenever the target or reference mode changes.
    observeEvent(list(
        dCq_rep_summary(),
        input$select_out_target,
        input$use_target_references
    ), {
        status <- reference_sample_status()
        req(nrow(status) >= 2)

        samples <- as.character(status$Sample)
        labels <- case_when(
            status$has_undetected &
                (status$has_failed_hk | status$has_missing_target) ~
                paste0(samples, " (undetected and missing replicates)"),
            status$has_undetected ~
                paste0(samples, " (contains undetected)"),
            status$has_failed_hk | status$has_missing_target ~
                paste0(samples, " (missing replicates)"),
            .default = samples
        )
        choices <- stats::setNames(samples, labels)
        samples_metadata <- hot_to_r(input$samples_tab)
        ordered <- as.character(samples_metadata$New_Label)
        first_available <- ordered[ordered %in% samples][1]

        global_reference <- isolate(cache$global_reference_sample)
        if (is.null(global_reference) || !global_reference %in% samples) {
            global_reference <- first_available
            cache$global_reference_sample <- global_reference
        }

        if (isTRUE(input$use_target_references)) {
            target <- as.character(input$select_out_target)
            target_references <- isolate(cache$target_reference_samples)
            selected <- unname(target_references[target])
            if (length(selected) == 0 || is.na(selected) || !selected %in% samples) {
                selected <- global_reference
                cache$target_reference_samples[target] <- selected
            }
        } else {
            selected <- global_reference
        }

        updateSelectInput(session, "reference_sample", choices = choices, selected = selected)
    })

    reference_sample_complete <- reactive({
        req(input$reference_sample, input$select_out_target)
        status <- reference_sample_status() |>
            filter(as.character(Sample) == input$reference_sample)
        nrow(status) == 1 && status$reference_complete
    })

    complete_reference_samples <- reactive({
        reference_sample_status() |>
            filter(reference_complete) |>
            pull(Sample) |>
            as.character()
    })

    # Resolve one explicit reference for every included target. In per-target
    # mode, targets without an override inherit the cached global reference.
    reference_assignments <- reactive({
        targets <- dCq_rep_summary() |>
            pull(Target) |>
            as.character() |>
            unique()

        global_reference <- cache$global_reference_sample %||% input$reference_sample
        req(global_reference)

        assigned_references <- rep(global_reference, length(targets))
        if (isTRUE(input$use_target_references)) {
            target_references <- cache$target_reference_samples[targets]
            has_override <- !is.na(target_references) & nzchar(target_references)
            assigned_references[has_override] <- unname(target_references[has_override])
        }

        tibble(
            Target = targets,
            Reference_Sample = assigned_references
        )
    })

    # Render the current per-target assignments as a compact table. The active
    # target is highlighted so the change that triggered the popup is clear.
    reference_assignment_table <- function(assignments = reference_assignments()) {
        rows <- lapply(seq_len(nrow(assignments)), function(index) {
            target <- assignments$Target[index]
            tags$tr(
                class = if (identical(
                    target, as.character(input$select_out_target)
                )) "table-primary" else NULL,
                tags$td(target),
                tags$td(assignments$Reference_Sample[index])
            )
        })

        div(
            class = "table-responsive",
            tags$table(
                class = "table table-sm table-striped align-middle mb-0",
                tags$thead(tags$tr(
                    tags$th("Target"),
                    tags$th("Reference sample")
                )),
                tags$tbody(rows)
            )
        )
    }

    show_reference_assignments <- function() {
        showModal(modalDialog(
            title = "Reference samples by target",
            tags$p("Current per-target reference assignments:"),
            reference_assignment_table(),
            easyClose = TRUE,
            footer = modalButton("Close")
        ))
    }

    # Cache an explicit selector change. In per-target mode, show the complete
    # assignment table only when the stored value for the active target really
    # changed; restoring a cached value after switching targets stays silent.
    observeEvent(input$reference_sample, {
        req(input$reference_sample)

        if (isTRUE(input$use_target_references)) {
            req(input$select_out_target)
            target <- as.character(input$select_out_target)
            previous_reference <- unname(cache$target_reference_samples[target])
            assignment_changed <- length(previous_reference) == 0 ||
                is.na(previous_reference) ||
                !identical(previous_reference, input$reference_sample)

            cache$target_reference_samples[target] <- input$reference_sample

            selected_status <- reference_sample_status() |>
                filter(as.character(Sample) == input$reference_sample)
            selected_is_complete <- nrow(selected_status) == 1 &&
                selected_status$reference_complete

            # Incomplete selections use the recommendation modal below, which
            # includes this same table and avoids opening two competing dialogs.
            if (assignment_changed && selected_is_complete) {
                show_reference_assignments()
            }
        } else {
            cache$global_reference_sample <- input$reference_sample
        }
    }, ignoreInit = TRUE, priority = 100)

    # Build one reference row for every expected biological replicate. A
    # censored, missing, or HK-failed reference run has no usable dCq value, but
    # the other runs for the same target remain available for normalization.
    reference_sample_dCq <- reactive({
        req(nrow(dCq_rep_summary()) > 0, hk_sample_status())

        expected_reference_runs <- reference_assignments() |>
            left_join(
                hk_sample_status() |>
                    transmute(
                        Reference_Sample = as.character(Sample),
                        across(any_of("Replicate")),
                        HK_valid
                    ),
                by = "Reference_Sample",
                relationship = "many-to-many"
            )

        reference_rows <- dCq_rep_summary() |>
            transmute(
                Target = as.character(Target),
                Reference_Sample = as.character(Sample),
                across(any_of("Replicate")),
                dCq_mean,
                dCq_censored,
                dCq_sd,
                dCq_se
            )

        reference_run_keys <- intersect(
            c("Target", "Reference_Sample", "Replicate"),
            intersect(names(expected_reference_runs), names(reference_rows))
        )

        expected_reference_runs |>
            left_join(reference_rows, by = reference_run_keys) |>
            group_by(Target, Reference_Sample) |>
            mutate(
                reference_value_available = HK_valid %in% TRUE &
                    !is.na(dCq_mean) & !(dCq_censored %in% TRUE),
                # Exported at target level: FALSE means at least one biological
                # replicate could not use the assigned reference.
                Reference_Valid = all(reference_value_available),
                ref_dCq_mean = if_else(
                    reference_value_available, dCq_mean, NA_real_
                ),
                # Missing/HK-failed runs are NA, not censored. This flag is
                # TRUE only when the reference target itself was undetected.
                ref_dCq_censored = coalesce(dCq_censored, FALSE),
                ref_dCq_sd = if_else(
                    reference_value_available, dCq_sd, NA_real_
                ),
                ref_dCq_se = if_else(
                    reference_value_available, dCq_se, NA_real_
                )
            ) |>
            ungroup() |>
            select(
                Target,
                Reference_Sample,
                Reference_Valid,
                any_of("Replicate"),
                ref_dCq_mean,
                ref_dCq_censored,
                ref_dCq_sd,
                ref_dCq_se
            )
    })

    # Target-level validity remains internal; exports only need the selected
    # reference name beside the corresponding ddCq result.
    reference_target_status <- reactive({
        reference_sample_dCq() |>
            distinct(Target, Reference_Sample, Reference_Valid)
    })
    # Keep reference-dependent metrics selectable even when some reference
    # runs are unavailable. Those runs become NA and are omitted downstream.

    observeEvent(list(input$select_out_target, input$reference_sample, reference_sample_status()), {
        req(input$select_out_target, input$reference_sample)

        updateRadioGroupButtons(
            session,
            "out_metric",
            disabledChoices = character(0)
        )
        updateRadioGroupButtons(
            session,
            "stats_metric",
            disabledChoices = character(0)
        )
    })

    # Recommend a complete reference as soon as an incomplete one is selected,
    # while still allowing the user to keep it. Only unavailable biological
    # replicates will then be excluded from reference-dependent analyses.
    observeEvent(
        list(input$reference_sample, input$select_out_target,
             input$use_target_references,
             reference_sample_status()),
        {
            req(input$reference_sample, input$select_out_target)

            # Wait until the selector and cache agree. This avoids prompting
            # for the previous target's reference while a cached selection is
            # being restored after a target switch.
            active_cached_reference <- if (isTRUE(input$use_target_references)) {
                unname(cache$target_reference_samples[
                    as.character(input$select_out_target)
                ])
            } else {
                cache$global_reference_sample
            }
            req(
                length(active_cached_reference) == 1,
                !is.na(active_cached_reference),
                identical(
                    as.character(input$reference_sample),
                    as.character(active_cached_reference)
                )
            )
            req(!reference_sample_complete())

            prompt_key <- paste(
                input$select_out_target,
                input$reference_sample,
                if (isTRUE(input$use_target_references)) "per-target" else "global",
                sep = "::"
            )

            # do not prompt again if the same reference has already been prompted for this target
            req(!identical(cache$reference_prompt_key, prompt_key))
            cache$reference_prompt_key <- prompt_key
            alternatives <- setdiff(
                complete_reference_samples(), input$reference_sample
            )

            selected_status <- reference_sample_status() |>
                filter(as.character(Sample) == input$reference_sample)
            issue <- case_when(
                selected_status$has_undetected &&
                    (selected_status$has_failed_hk |
                        selected_status$has_missing_target) ~
                    "has undetected target values and missing biological replicates",
                selected_status$has_undetected ~
                    "has undetected target values in one or more biological replicates",
                selected_status$has_failed_hk ~
                    "has one or more biological replicates with an undetected housekeeping gene",
                .default =
                    "is missing target values in one or more biological replicates"
            )

            consequence <- glue(
                "If you keep it, affected replicates will have invalid ΔΔCq values and will be excluded from ΔΔCq analyses and ANCOVA for {input$select_out_target}."
            )

            assignment_summary <- if (isTRUE(input$use_target_references)) {
                tagList(
                    tags$hr(),
                    tags$p(tags$strong("Current per-target reference assignments:")),
                    reference_assignment_table()
                )
            }

            body <- if (length(alternatives) > 0) {
                reference_scope <- if (isTRUE(input$use_target_references)) {
                    glue("This changes the cached reference only for {input$select_out_target}.")
                } else {
                    "This changes the global reference used by every target, analysis, and export."
                }
                tagList(
                    tags$p(glue(
                        "'{input$reference_sample}' {issue} for {input$select_out_target}. We recommend switching to a complete reference sample."
                    )),
                    tags$p(consequence),
                    tags$p(class = "text-muted", reference_scope),
                    assignment_summary,
                    selectInput(
                        "reference_sample_modal",
                        "New reference sample",
                        choices = alternatives,
                        selected = alternatives[1]
                    )
                )
            } else {
                tagList(
                    tags$p(glue(
                        "'{input$reference_sample}' {issue} for {input$select_out_target}. No complete alternative reference is available."
                    )),
                    tags$p(consequence),
                    assignment_summary
                )
            }

            showModal(modalDialog(
                title = "Incomplete reference sample",
                body,
                easyClose = TRUE,
                footer = if (length(alternatives) > 0) {
                    tagList(
                        modalButton("Keep current reference"),
                        actionButton(
                            "confirm_reference_sample",
                            "Switch reference",
                            class = "btn-primary"
                        )
                    )
                } else {
                    modalButton("Keep current reference")
                }
            ))
        },
        ignoreInit = TRUE
    )

    observeEvent(input$confirm_reference_sample, {
        req(input$reference_sample_modal)
        selected_reference <- input$reference_sample_modal
        updateSelectInput(session, "reference_sample", selected = selected_reference)
        cache$reference_prompt_key <- NULL
        removeModal()
    })
    
    # Reactive: ddCq data points -----------------------------------------------
    
    ddCq_data <- reactive({
        req(reference_sample_dCq())

        data <- dCq_data() |>
            mutate(Target = as.character(Target)) |>
            left_join(reference_assignments(), by = "Target") |>
            left_join(
                reference_target_status() |>
                    select(Target, Reference_Sample, Reference_Valid),
                by = c("Target", "Reference_Sample")
            )
        reference_join_cols <- intersect(
            c("Target", "Reference_Sample", "Replicate"),
            names(data)
        )

        data |>
            left_join(
                reference_sample_dCq() |> select(-Reference_Valid),
                by = reference_join_cols
            ) |>
            mutate(
                ddCq = dCq - ref_dCq_mean,
                ddCq_censored = dCq_censored | coalesce(ref_dCq_censored, FALSE),
                exp_ddCq = 2^-ddCq
            )
    })
    
    ddCq_rep_summary <- reactive({
        req(reference_sample_dCq())

        data <- dCq_rep_summary() |>
            mutate(Target = as.character(Target)) |>
            left_join(reference_assignments(), by = "Target") |>
            left_join(
                reference_target_status() |>
                    select(Target, Reference_Sample, Reference_Valid),
                by = c("Target", "Reference_Sample")
            )
        reference_join_cols <- intersect(
            c("Target", "Reference_Sample", "Replicate"),
            names(data)
        )

        data |>
            left_join(
                reference_sample_dCq() |> select(-Reference_Valid),
                by = reference_join_cols
            ) |>
            mutate(
                ddCq_mean = dCq_mean - ref_dCq_mean,
                ddCq_censored = dCq_censored | coalesce(ref_dCq_censored, FALSE),
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
            group_by(across(c(
                "Sample", "Target", "Reference_Sample", "Reference_Valid"
            ))) |>
            summarize(
                ddCq_n    = count_non_missing(ddCq_mean),
                ddCq_undetected_n = sum(ddCq_censored, na.rm = TRUE),
                ddCq_sd   = sd(ddCq_mean, na.rm = TRUE),
                ddCq_se   = ddCq_sd / sqrt(ddCq_n),
                # mean of means
                ddCq_mean    = mean_or_na(ddCq_mean),
                ddCq_censored = any(ddCq_censored, na.rm = TRUE),
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
            )
    })
    # Output Flags: Conditional panel visibility for statistical tests ---------
    
    stats_ui_flags <- reactive({
        req(input$stats_test)
        
        # Determine all UI flags based on selected test
        flags <- switch(input$stats_test,
                        # --- ANCOVA ---
                        "ancova" = list(
                            show_unequal_variance_toggle    = FALSE,
                            show_multiple_comparison_type   = TRUE,
                            show_multiple_comparison_adjust = FALSE,
                            show_post_hoc_test              = TRUE
                        ),

                        # --- ANCOVA 2 Sample ---
                        "ancova_2_sample" = list(
                            show_unequal_variance_toggle    = FALSE,
                            show_multiple_comparison_type   = FALSE,
                            show_multiple_comparison_adjust = FALSE,
                            show_post_hoc_test              = FALSE
                        ),
                        
                        # --- Mixed Effect Model ---
                        "mixed_effect" = list(
                            show_unequal_variance_toggle    = TRUE,
                            show_multiple_comparison_type   = TRUE,
                            show_multiple_comparison_adjust = FALSE,
                            show_post_hoc_test              = TRUE
                        ),


                        # --- ANCOVA 2 Sample ---
                        "mixed_effect_2_sample" = list(
                            show_unequal_variance_toggle    = TRUE,
                            show_multiple_comparison_type   = FALSE,
                            show_multiple_comparison_adjust = FALSE,
                            show_post_hoc_test              = FALSE
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
                        "repeated_ttest" = list(
                            show_unequal_variance_toggle    = TRUE,
                            show_multiple_comparison_type   = TRUE,
                            show_multiple_comparison_adjust = TRUE,
                            show_post_hoc_test              = FALSE
                        ),

                        # --- t-test ---
                        "ttest" = list(
                            show_unequal_variance_toggle    = TRUE,
                            show_multiple_comparison_type   = FALSE,
                            show_multiple_comparison_adjust = FALSE,
                            show_post_hoc_test              = FALSE
                        ),
                        
                        # --- Pairwise paired t-test ---
                        "repeated_paired_ttest" = list(
                            show_unequal_variance_toggle    = FALSE,
                            show_multiple_comparison_type   = TRUE,
                            show_multiple_comparison_adjust = TRUE,
                            show_post_hoc_test              = FALSE
                        ),

                        # --- Paired t-test ---
                        "paired_ttest" = list(
                            show_unequal_variance_toggle    = FALSE,
                            show_multiple_comparison_type   = FALSE,
                            show_multiple_comparison_adjust = FALSE,
                            show_post_hoc_test              = FALSE
                        ),
                        
                        # --- Pairwise Wilcoxon (signed-rank) ---
                        "repeated_wilcoxon" = list(
                            show_unequal_variance_toggle    = FALSE,
                            show_multiple_comparison_type   = TRUE,
                            show_multiple_comparison_adjust = TRUE,
                            show_post_hoc_test              = FALSE
                        ),

                        # --- Wilcoxon ---
                        "wilcoxon" = list(
                            show_unequal_variance_toggle    = FALSE,
                            show_multiple_comparison_type   = FALSE,
                            show_multiple_comparison_adjust = FALSE,
                            show_post_hoc_test              = FALSE
                        ),
                        
                        # --- Pairwise Wilcoxon-Mann-Whitney ---
                        "repeated_mann_whitney" = list(
                            show_unequal_variance_toggle    = FALSE,
                            show_multiple_comparison_type   = TRUE,
                            show_multiple_comparison_adjust = TRUE,
                            show_post_hoc_test              = FALSE
                        ),

                        # --- Mann-Whitney ---
                        "mann_whitney" = list(
                            show_unequal_variance_toggle    = FALSE,
                            show_multiple_comparison_type   = FALSE,
                            show_multiple_comparison_adjust = FALSE,
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
    
    observeEvent(list(input$stats_metric, n_samples(), n_numeric_samples(), n_bio_reps()), {
        req(input$stats_metric)
        req(n_samples() >= 2)
        req(n_bio_reps() >= 2)
        
        metric <- input$stats_metric
        # Non-parametric tests require at least 5 biological replicates
        include_nonparam <- n_bio_reps() >= 5

        # Number of samples
        n_samples <- n_samples()
        # Number of viable samples for parametric test
        n_numeric_samples <- n_numeric_samples()
        
        # Determine available tests and default based on metric and sample count
        # parametric choices
        choices <- list()
        default <- NULL

        # parametric tests choices
        if (n_numeric_samples > 2) {
            # > 2 samples
            if (metric == "dCq") {
                # dCq: comparison across samples accounting for HK variance
                choices[["Parametric"]] <- c(
                        "ANCOVA" = "ancova",
                        "Mixed Effect Model" = "mixed_effect",
                        "Repeated paired t-test" = "repeated_paired_ttest"
                    )
                default <- "ancova"
            } else {
                # ddCq or exp_ddCq: standard group comparisons
                choices[["Parametric"]] <- c(
                        "Repeated t-test" = "repeated_ttest",
                        "ANOVA" = "anova"
                    )
                default <- "repeated_ttest"
            }
        } else {
            # = 2 samples
            if (metric == "dCq") {
                choices[["Parametric"]] <- c(
                        "ANCOVA" = "ancova_2_sample",
                        "Mixed Effect Model" = "mixed_effect_2_sample",
                        "Paired t-test" = "paired_ttest"
                    )
                default <- "ancova_2_sample"
            } else {
                # ddCq or exp_ddCq
                choices[["Parametric"]] <- c("t-test" = "ttest")
                default <- "ttest"
            }
        }

        # Non-parametric test choices
        if (include_nonparam) {
            if (n_samples > 2) {
                if (metric == "dCq") {
                    choices[["Non-parametric"]] <- c("Repeated Wilcoxon signed-rank (paired)" = "repeated_wilcoxon")
                } else {
                    choices[["Non-parametric"]] <- c(
                        "Kruskal-Wallis" = "kruskal",
                        "Repeated Wilcoxon-Mann-Whitney" = "repeated_mann_whitney"
                    )
                }
            } else {
                if (metric == "dCq") {
                    choices[["Non-parametric"]] <- c("Wilcoxon signed-rank (paired)" = "wilcoxon")
                } else {
                    choices[["Non-parametric"]] <- c("Wilcoxon-Mann-Whitney" = "mann_whitney")
                }
            }
        }
        
        updatePickerInput(session, "stats_test", choices = choices, selected = default)
    })

    # Default to unequal variance whenever an independent t-test is selected:
    # repeated pairwise Welch tests for >2 samples, or Welch's test for 2 samples.
    observeEvent(input$stats_test, {
        req(input$stats_test)
        updatePrettySwitch(
            session,
            "stats_unequal_variance",
            value = input$stats_test %in% c("repeated_ttest", "ttest")
        )
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
        req(input$stats_metric)
        req(input$stats_test)
        req(input$select_out_target)
        req(n_bio_reps() >= 2)
        req(n_samples() >= 2)
        
        # statistical test settings
        test       <- input$stats_test
        response   <- input$stats_metric
        
        # Validate that the current test is compatible with the current metric.
        # After stats_metric changes, there is a brief moment where stats_test
        # still holds the old value (until the observer updates it).  Silently
        # abort here so we never run an incompatible test/metric combination.
        dCq_tests     <- c("ancova", "ancova_2_sample", "mixed_effect",
                           "mixed_effect_2_sample", "paired_ttest",
                           "repeated_paired_ttest", "repeated_wilcoxon", "wilcoxon")
        non_dCq_tests <- c("anova", "ttest", "repeated_ttest",
                           "kruskal", "repeated_mann_whitney", "mann_whitney")
        valid_tests   <- if (response == "dCq") dCq_tests else non_dCq_tests
        req(test %in% valid_tests)
        equal_var  <- !isTRUE(input$stats_unequal_variance)
        comparison <- input$stats_comparison
        p_adjust   <- input$stats_multiple_comparison_adjust
        
        # Prepare data based on metric
        if (response == "dCq") {
            data <- dCq_rep_summary() |>
                filter(Target == input$select_out_target) |>
                rename(dCq = dCq_mean)

            if (test %in% c("ancova", "ancova_2_sample")) {
                data <- data |>
                    left_join(reference_sample_dCq(), by = c("Target", "Replicate")) |>
                    rename(ref_dCq = ref_dCq_mean) |>
                    drop_na(dCq, ref_dCq)
            } else {
                data <- data |> drop_na(dCq)
            }
        } else {
            data <- ddCq_rep_summary() |>
                filter(Target == input$select_out_target) |>
                rename(ddCq = ddCq_mean, exp_ddCq = exp_ddCq_mean) |>
                drop_na(ddCq, exp_ddCq)
        }

        censor_flag <- paste0(response, "_censored")
        n_censored_points <- if (censor_flag %in% names(data)) {
            sum(data[[censor_flag]], na.rm = TRUE)
        } else {
            0
        }

        # Make the selected reference the control level for Dunnett/control
        # comparisons as well as for the ΔΔCq calculation.
        data <- data |>
            mutate(Sample = forcats::fct_relevel(Sample, input$reference_sample)) |>
            droplevels()

        # Recalculate available samples after filtering
        n_samples <- data$Sample |> unique() |> length()
        req(n_samples > 1)

        # Run appropriate test based on selection
        tryCatch({
            result <- switch(test,
                             "ancova" = run_ancova(data, response = response, comparison = comparison),
                             "ancova_2_sample" = run_ancova_2_sample(data, response = response),
                             "mixed_effect" = run_mixed_effect(data, response = response, comparison = comparison, equal.var = equal_var),
                             "mixed_effect_2_sample" = run_mixed_effect_2_sample(data, response = response, equal.var = equal_var),    
                             "anova" = run_anova(data, response = response,comparison = comparison),
                             "kruskal" = run_kruskal(data, response = response, comparison = comparison, p_adjust_method = p_adjust),
                             "repeated_ttest" = run_repeated_ttest(data, response = response, comparison = comparison, equal.var = equal_var, p_adjust_method = p_adjust),
                             "ttest" = run_ttest(data, response = response, equal.var = equal_var),
                             "repeated_paired_ttest" = run_repeated_paired_ttest(data, response = response, comparison = comparison, p_adjust_method = p_adjust),
                             "paired_ttest" = run_paired_ttest(data, response = response),
                             "repeated_wilcoxon" = run_repeated_wilcoxon(data, response = response, comparison = comparison, p_adjust_method = p_adjust),
                             "wilcoxon" = run_wilcoxon(data, response = response),
                             "repeated_mann_whitney" = run_repeated_mann_whitney(data, response = response, comparison = comparison, p_adjust_method = p_adjust),
                             "mann_whitney" = run_mann_whitney(data, response = response)
            )
            result$n_censored_points <- n_censored_points
            result$n_dropped_points <- 0
            result$n_dropped_reps <- 0
            result$n_dropped_samples <- 0
            result
        }, error = function(e) {
            list(error = as.character(e$message), n_censored_points = n_censored_points,
                 n_dropped_points = 0, n_dropped_reps = 0, n_dropped_samples = 0)
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
    
    # Output: Omnibus badge (brief p-value indicator) --------------------------
    
    output$stats_omnibus_badge <- renderUI({
        req(stats_result()$omnibus_res)
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
        req(stats_result()$omnibus_res)
        req(stats_result()$omnibus_label)
        stats_result()$omnibus_label |>
            # strip p-value that will be added in a label
            str_remove(", p = [0-9\\-e\\.]+$")
        
    })
    
    # Output: Omnibus test table -----------------------------------------------
    
    output$stats_omnibus_table <- DT::renderDataTable({
        req(stats_result())
        req(stats_result()$omnibus_res)
        
        stats_result()$omnibus_res |>
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
                columns = select(stats_result()$omnibus_res, where(is.numeric)) %>% names(),
                digits = 2
            ) 
    })
    
    
    # Output: Flag for additional test panel (for conditional UI) ------------------
    output$stats_has_extra <- reactive({
        req(stats_result())
        !is.null(stats_result()$extra_res)
    })
    outputOptions(output, "stats_has_extra", suspendWhenHidden = FALSE)
    
    output$stats_extra_in_omnibus <- reactive({
        req(stats_result())
        stats_result()$extra_position == "omnibus"
    })
    outputOptions(output, "stats_extra_in_omnibus", suspendWhenHidden = FALSE)
    
    # Output: Optional additional res panel ---------------------------------------------
    
    stats_extra_title_render <- renderText({
        req(stats_result())
        req(stats_result()$extra_label)
        
        stats_result()$extra_label
    })
    
    output$stats_extra_title_omnibus    <- stats_extra_title_render
    output$stats_extra_title_comparison <- stats_extra_title_render

    stats_extra_table_render <- DT::renderDataTable({
        req(stats_result())
        req(stats_result()$extra_res)
        
        stats_result()$extra_res |>
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
            )
    })
    
    output$stats_extra_table_omnibus    <- stats_extra_table_render
    output$stats_extra_table_comparison <- stats_extra_table_render
    
    # Output: Comparison section title -----------------------------------------
    
    output$stats_comparison_title <- renderText({
        req(stats_result())
        req(stats_result()$test_label)
        
        if (!is.null(stats_result()$omnibus_res)) {
            paste("Post-hoc:", stats_result()$test_label)
        } else {
            stats_result()$test_label
        }
    })
    
    # Output: Results table (post-hoc or pairwise) -----------------------------
    
    output$stats_results_table <- DT::renderDataTable({
        req(stats_result())
        
        # Get the comparison results table
        results_df <- stats_result()$test_res
        req(results_df)
        
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
            DT::formatSignif(
                columns = setdiff(
                    select(results_df, where(is.numeric)) %>% names(),
                    c("n1", "n2")
                    ),
                digits = 2
            ) |>
            DT::formatRound(
                columns = intersect(
                    select(results_df, where(is.numeric)) %>% names(),
                    c("n1", "n2")
                ),
                digits = 0) |>
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
            "dCq" ~ "log₂ relative normalized expression (-ΔCq)",
            "exp_dCq" ~ "relative normalized expression (2^-ΔCq)",
            "ddCq" ~ "log₂ fold change (-ΔΔCq)",
            "exp_ddCq" ~ "fold change (2^-ΔΔCq)",
        )
        
        paste(input$select_out_target, " ", y_label)
    })
    
    # Shared Reactive: Plot data preparation (used by Results + Export) --------
    result_plot_data <- reactive({
        req(input$select_out_target)
        req(input$stat_type)
        req(input$out_metric)
        req(nrow(dCq_data()) > 0)
        req(nrow(dCq_rep_summary()) > 0)
        
        # Select data based on metric and bio-rep summarisation mode
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
            "dCq" ~ "Log₂ norm. expression (-ΔCq)",
            "exp_dCq" ~ "Normalized expression (2^-ΔCq)",
            "ddCq" ~ "Log₂ fold change (-ΔΔCq)",
            "exp_ddCq" ~ "Fold change (2^-ΔΔCq)"
        )
        
        # Sign inversion for -dCq / -ddCq display
        sign <- if (input$out_metric %in% c("dCq", "ddCq")) -1 else 1
        
        y_summary_value <- glue("{input$out_metric}_mean")
        
        # Error bar column names
        if (input$stat_type == "none") {
            error_bar_high  <- NULL
            error_bar_low   <- NULL
            error_bar_label <- ""
        } else {
            error_bar_high <- glue("{input$out_metric}_{input$stat_type}_high")
            error_bar_low  <- glue("{input$out_metric}_{input$stat_type}_low")
            error_bar_label <- glue("{str_to_upper(input$stat_type)}: ({round(df_summary_target[[error_bar_low]], 2)}; {round(df_summary_target[[error_bar_high]], 2)})")
        }
        
        # Compute y limits
        values <- df_target[[y_value]]
        if (input$stat_type != "none") {
            values <- values |>
                append(c(
                    df_summary_target |> pull(error_bar_low),
                    df_summary_target |> pull(error_bar_high)
                ))
        }
        
        y_limits <- get_y_limits(
            values,
            metric = input$out_metric,
            undetected_value = max_cycle_value()
        )

        point_censored_col <- censoring_column_for(y_value)
        summary_censored_col <- censoring_column_for(y_summary_value)
        req(point_censored_col %in% names(df_target))
        req(summary_censored_col %in% names(df_summary_target))

        point_censored <- df_target[[point_censored_col]]
        point_censored[is.na(point_censored)] <- FALSE
        summary_censored <- df_summary_target[[summary_censored_col]]
        summary_censored[is.na(summary_censored)] <- FALSE
        undetected_present <- any(point_censored)
        y_min_label <- as.character(y_limits[1])

        # Mark point types and add censoring symbols to plot hover text.
        df_target <- df_target |>
            mutate(
                point_type_label = ifelse(point_censored, "Undetected", "Detected"),
                plot_value_display = format_censored_value(
                    sign * .data[[y_value]], point_censored,
                    censoring_direction_for(y_value, sign), digits = 2
                )
            )
        df_summary_target <- df_summary_target |>
            mutate(
                summary_value_display = format_censored_value(
                    sign * .data[[y_summary_value]], summary_censored,
                    censoring_direction_for(y_summary_value, sign), digits = 2
                )
            )
        
        list(
            df_target          = df_target,
            df_summary         = df_summary_target,
            y_value            = y_value,
            y_summary_value    = y_summary_value,
            y_label            = y_label,
            sign               = sign,
            error_bar_high     = error_bar_high,
            error_bar_low      = error_bar_low,
            y_limits           = y_limits,
            y_min_label        = y_min_label,
            undetected_present = undetected_present,
            error_bar_label    = error_bar_label,
            out_metric         = input$out_metric,
            stat_type          = input$stat_type,
            summarize_bio_reps = input$summarize_bio_reps,
            target_name        = input$select_out_target,
            y_axis_mode        = input$y_axis_mode %||% "auto",
            custom_y_limits    = if ((input$y_axis_mode %||% "auto") == "custom")
                                     c(input$y_axis_min, input$y_axis_max) else NULL
        )
    })
    
    output$res_plot <- renderPlotly({
        d <- result_plot_data()
        p <- build_results_plot(d, accent_color(), secondary_color())
        
        ggplotly(p, tooltip = "text") |>
            fix_plotly_legend()
    })
    
    # Export Tab Server Logic ==================================================
    
    # Reactive: Get unique samples for color pickers
    export_samples <- reactive({
        req(cq_data())
        levels(cq_data()$Sample)
    })
    
    # Output: Dynamic sample color inputs
    output$sample_color_inputs <- renderUI({
        samples <- export_samples()
        req(length(samples) > 0)

        # fall back color if accent is not workign
        color <- if (!is.null(accent_color())) accent_color() else "#027BC3"
        
        # Generate a color palette for samples
        default_colors <- rep(color, length(samples))
        
        tagList(
            lapply(seq_along(samples), function(i) {
                colourInput(
                    inputId = paste0("sample_color_", i),
                    label = samples[i],
                    value = default_colors[i],
                    showColour = "both",
                    palette = "square"
                )
            })
        )
    })
    
    sample_colors <- reactive({
        samples <- export_samples()
        req(length(samples) > 0)
        
        color <- if (!is.null(accent_color())) accent_color() else "#027BC3"
        default_colors <- rep(color, length(samples))
        
        colors <- sapply(seq_along(samples), function(i) {
            input[[paste0("sample_color_", i)]] %||% default_colors[i]
        })
        setNames(colors, samples)
    })
    
    # Update default plot dimensions based on sample count (1 cm per sample, 4 cm height)
    observeEvent(n_samples(), {
        req(n_samples() >= 1)
        updateNumericInput(session, "export_plot_width", value = n_samples())
    })
    
    # Reactive: Generate export plot
    export_plot_obj <- reactive({
        d <- result_plot_data()
        
        build_export_plot(
            plot_data      = d,
            colors         = sample_colors(),
            lw             = input$export_linewidth %||% 0.5,
            point_size     = input$export_point_size %||% 2,
            axis_text_size = input$export_axis_text_size %||% 10,
            signif_text_size = input$export_signif_text_size %||% 8,
            bar_width      = input$export_bar_width %||% 0.6,
            plot_width     = input$export_plot_width %||% 4,
            plot_height    = input$export_plot_height %||% 4,
            show_signif_bars = isTRUE(input$show_signif_bars),
            stats_result   = stats_result(),
            hide_ns        = isTRUE(input$hide_ns_bars),
            show_exact_pvalue = isTRUE(input$show_exact_pvalue)
        )
    })
    
    # Output: Export plot
    output$export_plot <- renderPlot({
        export_plot_obj()
    })
    # Export Data Reactives (shared by preview tables and XLSX download) ========
    
    # Format each exported metric in place. Internal numeric values and
    # censoring flags are not exported as additional columns.
    finalize_export_metrics <- function(df, metrics, digits = 4) {
        metrics <- metrics[metrics %in% names(df)]
        for (metric in metrics) {
            censor_col <- censoring_column_for(metric)

            # Housekeeping-gene means should never be censored: samples without a
            # detected HK have already been excluded from ΔCq processing.
            if (is.null(censor_col)) {
                next
            }

            censored <- if (censor_col %in% names(df)) {
                df[[censor_col]]
            } else {
                rep(FALSE, nrow(df))
            }
            censored[is.na(censored)] <- FALSE

            df[[metric]] <- format_censored_value(
                df[[metric]], censored,
                censoring_direction_for(metric), digits = digits
            )
        }
        df
    }

    # Label the target measurement itself. Reference censoring may make ddCq
    # unavailable, but it must not label a detected target row as undetected.
    add_undetected_label <- function(df) {
        df |>
            mutate(
                Undetected = if_else(
                    Cq_censored %in% TRUE,
                    "Undetected",
                    ""
                )
            )
    }
    
    export_raw_cq <- reactive({
        raw_data <- current_raw_data()
        req(raw_data)
        req(nrow(hot_to_r(input$samples_tab)) > 0)
        req(nrow(hot_to_r(input$targets_tab)) > 0)
        
        samples_metadata <- hot_to_r(input$samples_tab)
        targets_metadata <- hot_to_r(input$targets_tab)

        raw_data |>
            parse_Cq_data() |>
            mutate(
                Cq = replace_censored(
                    Cq,
                    censored = Cq_censored,
                    replacement = max_cycle_value()
                )
            ) |>
            mutate(Key = row_number()) |>
            # join sample metadata
            left_join(
                samples_metadata |> rename(Sample_Label = New_Label, Sample_Include = Include),
                by = "Sample"
            ) |>
            # join target metadata
            left_join(
                targets_metadata |> rename(Target_Label = New_Label, Target_Include = Include),
                by = "Target"
            ) |>
            # apply renames (use original name as fallback for unmatched)
            mutate(
                Sample = coalesce(Sample_Label, Sample),
                Target = coalesce(Target_Label, Target),
                Excluded = coalesce(!Sample_Include, FALSE) |
                           coalesce(!Target_Include, FALSE) |
                           Key %in% cache$excluded_point_keys
            ) |>
            select(-Key, -Sample_Label, -Sample_Include, -Target_Label, -Target_Include) |>
            finalize_export_metrics("Cq") |>
            select(Sample, Target, any_of("Replicate"),
                   Cq, Excluded)
    })
    
    export_technical <- reactive({
        req(ddCq_data())
        df <- ddCq_data() |>
            rename(
                HK_mean_Cq = HK_mean,
                ref_mean_dCq = ref_dCq_mean
            )

        hk_metrics <- names(df)[grepl("^HK_mean_.+_Cq$", names(df)) & names(df) != "HK_mean_Cq"]

        pre_ddCq_metrics <- c(
            "Cq", "HK_mean_Cq", hk_metrics, "dCq", "exp_dCq",
            "ref_mean_dCq"
        )
        ddCq_metrics <- c("ddCq", "exp_ddCq")
        metrics <- c(pre_ddCq_metrics, ddCq_metrics)

        df |>
            add_undetected_label() |>
            finalize_export_metrics(metrics) |>
            select(any_of("Replicate"), Sample, Target, Undetected,
                   any_of(pre_ddCq_metrics),
                   Reference_Sample,
                   any_of(ddCq_metrics))
    })
    
    export_bio_rep <- reactive({
        req(ddCq_rep_summary())
        
        df <- ddCq_rep_summary() |>
            mutate(Cq_n = as.integer(Cq_n)) |>
            rename(
                ref_mean_dCq = ref_dCq_mean
            )
        
        # Base columns always shown
        # Include individual HK gene average columns (HK_mean_<gene>_Cq) if present
        hk_indiv_cols <- names(df)[grepl("^HK_mean_.+_Cq$", names(df)) & names(df) != "HK_mean_Cq"]
        
        pre_ddCq_metrics <- c(
            "Cq_mean", "HK_mean_Cq", hk_indiv_cols,
            "dCq_mean", "exp_dCq_mean", "ref_mean_dCq"
        )
        ddCq_metrics <- c("ddCq_mean", "exp_ddCq_mean")
        metrics <- c(pre_ddCq_metrics, ddCq_metrics)
        base_cols <- c("Replicate", "Sample", "Target",
                       "Undetected", "Cq_n",
                       pre_ddCq_metrics,
                       "Reference_Sample",
                       ddCq_metrics)
        
        # Only include dispersion when single replicate (or no Replicate column)
        # and user has error bars enabled
        if (n_bio_reps() <= 1 && input$stat_type != "none") {
            stat <- input$stat_type  # "sd" or "se"
            disp_cols <- c(
                paste0("Cq_", stat),
                paste0("dCq_", stat),
                paste0("exp_dCq_", stat, c("_low", "_high")),
                paste0("ddCq_", stat),
                paste0("exp_ddCq_", stat, c("_low", "_high"))
            )
            base_cols <- c(base_cols, disp_cols)
        }
        
        df |>
            add_undetected_label() |>
            finalize_export_metrics(metrics) |>
            select(any_of(base_cols))
    })
    
    export_summary <- reactive({
        req(n_bio_reps() > 1)
        df <- dCq_summary() |>
            left_join(ddCq_summary(), by = c("Sample", "Target")) |>
            mutate(across(where(is.integer), as.integer)) |>
            select(-matches("^(dCq|ddCq)_(sd|se)_(low|high)$"))

        dCq_metrics <- c("dCq_mean", "exp_dCq_mean")
        ddCq_metrics <- c("ddCq_mean", "exp_ddCq_mean")
        metrics <- c(dCq_metrics, ddCq_metrics)

        df |>
            finalize_export_metrics(metrics) |>
            select(Sample, Target,
                   any_of(c("dCq_n", "dCq_undetected_n")),
                   any_of(dCq_metrics),
                   any_of(c("dCq_sd", "dCq_se")),
                   Reference_Sample,
                   any_of(c("ddCq_n", "ddCq_undetected_n")),
                   any_of(ddCq_metrics),
                   any_of(c("ddCq_sd", "ddCq_se")))
    })
    
    # Data Preview Tables ======================================================
    
    output$preview_raw_cq <- DT::renderDataTable({
        df <- export_raw_cq()
        df |>
            DT::datatable(
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE,
                class = "compact stripe"
            )
    })

    output$preview_technical <- DT::renderDataTable({
        df <- export_technical()
        df |>
            DT::datatable(
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE,
                class = "compact stripe"
            ) |>
            DT::formatRound(columns = names(df)[sapply(df, is.numeric)], digits = 4)
    })
    
    output$preview_bio_rep <- DT::renderDataTable({
        df <- export_bio_rep()
        df |>
            DT::datatable(
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE,
                class = "compact stripe"
            ) |>
            DT::formatRound(columns = names(df)[sapply(df, is.numeric)], digits = 4)
    })
    
    output$preview_summary <- DT::renderDataTable({
        df <- export_summary()
        df |>
            DT::datatable(
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE,
                class = "compact stripe"
            ) |>
            DT::formatRound(columns = names(df)[sapply(df, is.numeric)], digits = 4)
    })
    
    # Download Handlers ========================================================
    
    # Helper: compute total figure size from the built ggplot grob
    export_plot_dims <- reactive({
        get_plot_dims(export_plot_obj())
    })
    
    # Download: Plot as PNG
    output$download_plot_png <- downloadHandler(
        filename = function() {
            paste0("qPCR_plot_", input$select_out_target, "_", Sys.Date(), ".png")
        },
        content = function(file) {
            dims <- export_plot_dims()
            ggsave(
                file,
                plot = export_plot_obj(),
                width = dims$width,
                height = dims$height,
                units = "cm",
                dpi = 300,
                bg = "white"
            )
        }
    )
    
    # Download: Plot as PDF
    output$download_plot_pdf <- downloadHandler(
        filename = function() {
            paste0("qPCR_plot_", input$select_out_target, "_", Sys.Date(), ".pdf")
        },
        content = function(file) {
            dims <- export_plot_dims()
            ggsave(
                file,
                plot = export_plot_obj(),
                width = dims$width,
                height = dims$height,
                units = "cm",
                device = cairo_pdf
            )
        }
    )
    
    # Download: Data as XLSX
    output$download_data_xlsx <- downloadHandler(
        filename = function() {
            paste0("qPCR_data_", Sys.Date(), ".xlsx")
        },
        content = function(file) {
            wb <- createWorkbook()
            
            # Sheet 1: Raw Cq data
            addWorksheet(wb, "Raw Cq")
            writeData(wb, "Raw Cq", export_raw_cq())

            # Sheet 2: Technical Replicates
            addWorksheet(wb, "Technical Replicates")
            writeData(wb, "Technical Replicates", export_technical())
            
            # Sheet 3: Bio Rep Averages
            addWorksheet(wb, "Bio Rep Averages")
            writeData(wb, "Bio Rep Averages", export_bio_rep())
            
            # Sheet 4: Summary (if ≥2 bio reps)
            tryCatch({
                if (n_bio_reps() > 1) {
                    addWorksheet(wb, "Summary")
                    writeData(wb, "Summary", export_summary())
                }
            }, error = function(e) NULL)
            

            
            saveWorkbook(wb, file, overwrite = TRUE)
        }
    )
    
    # Download: Statistics as XLSX (Current Target)
    output$download_stats_xlsx <- downloadHandler(
        filename = function() {
            paste0("qPCR_stats_", input$select_out_target, "_", Sys.Date(), ".xlsx")
        },
        content = function(file) {
            wb <- createWorkbook()
            # Excel sheet name limit is 31 characters
            sheet_name <- substr(input$select_out_target, 1, 31)
            
            addWorksheet(wb, sheet_name)
            
            tryCatch({
                if (!is.null(stats_result()) && is.null(stats_result()$error)) {
                    row <- 1
                    
                    # Omnibus results
                    if (!is.null(stats_result()$omnibus_res)) {
                        omnibus_header <- paste("Omnibus:", stats_result()$omnibus_label %||% "")
                        writeData(wb, sheet_name, data.frame(x = omnibus_header), startRow = row, colNames = FALSE)
                        row <- row + 1
                        writeData(wb, sheet_name, stats_result()$omnibus_res, startRow = row)
                        row <- row + nrow(stats_result()$omnibus_res) + 2
                    }
                    
                    # Extra results (if placed after Omnibus)
                    if (!is.null(stats_result()$extra_res) && isTRUE(stats_result()$extra_position == "omnibus")) {
                         extra_header <- stats_result()$extra_label %||% "Additional Results"
                         writeData(wb, sheet_name, data.frame(x = extra_header), startRow = row, colNames = FALSE)
                         row <- row + 1
                         writeData(wb, sheet_name, stats_result()$extra_res, startRow = row)
                         row <- row + nrow(stats_result()$extra_res) + 2
                    }

                    # Test results
                    if (!is.null(stats_result()$test_res)) {
                        comp_header <- if (!is.null(stats_result()$omnibus_res)) {
                            paste("Post-hoc:", stats_result()$test_label %||% "")
                        } else {
                            stats_result()$test_label %||% "Pairwise Comparisons"
                        }
                        writeData(wb, sheet_name, data.frame(x = comp_header), startRow = row, colNames = FALSE)
                        row <- row + 1
                        writeData(wb, sheet_name, stats_result()$test_res, startRow = row)
                        row <- row + nrow(stats_result()$test_res) + 2
                    }
                    
                    # Extra results (if placed after Comparisons - Default)
                    if (!is.null(stats_result()$extra_res) && !isTRUE(stats_result()$extra_position == "omnibus")) {
                         extra_header <- stats_result()$extra_label %||% "Additional Results"
                         writeData(wb, sheet_name, data.frame(x = extra_header), startRow = row, colNames = FALSE)
                         row <- row + 1
                         writeData(wb, sheet_name, stats_result()$extra_res, startRow = row)
                         row <- row + nrow(stats_result()$extra_res) + 2
                    }
                    
                    # Method
                    if (!is.null(stats_result()$method)) {
                         writeData(wb, sheet_name, data.frame(x = stats_result()$method), startRow = row, colNames = FALSE)
                    }
                } else {
                     writeData(wb, sheet_name, "No valid statistical results available.")
                }
            }, error = function(e) {
                # Fallback to "Error" if sheet name is invalid, though tryCatch handles logic errors
                tryCatch({
                    writeData(wb, sheet_name, paste("Error exporting statistics:", e$message))
                }, error = function(e2) {
                     # specific case where sheet name might be the issue? unlikely if addWorksheet passed
                })
            })
            
            saveWorkbook(wb, file, overwrite = TRUE)
        }
    )

    # Plot Preview Title (Dynamic)
    output$plot_preview_title <- renderUI({
        req(input$select_out_target)
        tagList(
            "Plot Preview for gene ", 
            tags$strong(input$select_out_target)
        )
    })
}

# Run App ======================================================================

shinyApp(ui = ui, server = server)
