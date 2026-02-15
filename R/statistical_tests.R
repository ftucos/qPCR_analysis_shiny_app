library(tidyverse)
library(glue)
library(lme4)
library(broom)
library(lmerTest)
library(emmeans)
library(nlme)
library(PMCMRplus)

format_response <- function(response) {
    case_match(response,
        "dCq"      ~ "ΔCq",
        "ddCq"     ~ "ΔΔCq",
        "exp_ddCq" ~  "2^-ΔΔCq",
        .default = response
    )
}

# Helper: convert p-values to significance symbols
to_significance <- function(p) {
    case_when(
        p <= 0.001 ~ "***",
        p <= 0.01  ~ "**",
        p <= 0.05  ~ "*",
        TRUE       ~ "ns"
    )
}

# Helper: split contrast column to group1 and group2 for ggpubr compatibility
format_emmeans <- function(df) {
    df |>
        mutate(contrast = ifelse(
            str_detect(contrast, "\\(.*\\) - \\(.*\\)"),
            str_remove_all(contrast, "(^\\(|\\)$)") |> str_split("\\) - \\("),
            str_split(contrast, " - ")
        )) |>
        rename(group = contrast) |>
        unnest_wider(group, names_sep = "")
}

# Helper: format PMCMRplus matrix output to tidy format with group columns
format_pmcmr <- function(pmcmr_result) {
    PMCMRplus::toTidy(pmcmr_result) |>
        rename(group1 = group1, group2 = group2) |>
        mutate(Significance = to_significance(p.value))
}

# Helper: format pairwise.*.test output to tidy format
format_pairwise_test <- function(pairwise_result) {
    # Extract p-value matrix
    p_mat <- pairwise_result$p.value

    # Convert matrix to long format
    tibble(
        group1 = rep(rownames(p_mat), each = ncol(p_mat)),
        group2 = rep(colnames(p_mat), times = nrow(p_mat)),
        `p-value` = as.vector(t(p_mat))
    ) |>
        filter(!is.na(`p-value`)) |>
        mutate(Significance = to_significance(`p-value`))
}

# ANCOVA -----------------------------------------------------------------------
# For dCq with biological replicates, adjusting for reference sample variance

run_ancova <- function(x,
                       response = c("dCq"),
                       comparison = c("pairwise", "trt.vs.ctrl")) {

    response <- match.arg(response)
    comparison <- match.arg(comparison)

    test_formula <- dCq ~ Sample + ref_dCq

    # Omnibus: ANCOVA
    omnibus <- aov(test_formula, data = x)
    omnibus_method  <- "ANCOVA"
    omnibus_pvalue  <- summary(omnibus)[[1]]["Sample", "Pr(>F)"]

    omnibus_res <- broom::tidy(omnibus) |>
        rename(
            Term       = term,
            Df         = df,
            `Sum Sq`   = sumsq,
            `Mean Sq`  = meansq,
            `F-value`  = statistic,
            `p-value`  = p.value
        ) |>
        mutate(Term = str_replace(Term, "ref_dCq", "Covariate (Reference dCq)"))

    omnibus_label <- glue(
        "ANCOVA (adj. ref dCq): Sample F({sample_df},{residuals_df}) = {sample_f_value}, p = {prettyNum(signif(omnibus_pvalue, 2))}",
        sample_df      = omnibus_res$Df[1],
        residuals_df   = omnibus_res$Df[nrow(omnibus_res)],
        sample_f_value = round(omnibus_res$`F-value`[1], 2)
    )

    # Post-hoc: estimated marginal means
    emm <- emmeans(omnibus, ~Sample)

    if (comparison == "pairwise") {
        post_hoc <- contrast(emm, method = "pairwise", adjust = "tukey")
        post_hoc_res <- post_hoc |>
            broom::tidy() |>
            format_emmeans() |>
            rename(
                Term           = term,
                Df             = df,
                `q-value`      = statistic,
                `Adj. p-value` = adj.p.value
            ) |>
            mutate(Significance = to_significance(`Adj. p-value`)) |>
            select(-null.value, -estimate, -std.error)
        post_hoc_method <- "Tukey HSD"

    } else {
        set.seed(123)
        post_hoc <- contrast(emm, method = "trt.vs.ctrl", adjust = "mvt")
        post_hoc_res <- post_hoc |>
            broom::tidy() |>
            format_emmeans() |>
            rename(
                Term           = term,
                Df             = df,
                `t-value`      = statistic,
                `Adj. p-value` = adj.p.value
            ) |>
            mutate(Significance = to_significance(`Adj. p-value`)) |>
            select(-null.value, -estimate, -std.error)

        post_hoc_method <- "Dunnett"
    }

    method <- glue("Two-sided ANCOVA on {format_response(response)} with reference sample {format_response(response)} covariate adjustment, followed by {post_hoc_method} post-hoc test.")

    list(
        test_res        = post_hoc_res,
        test_label      = post_hoc_method,
        method          = method,
        omnibus_res     = omnibus_res,
        omnibus_label   = omnibus_label,
        omnibus_pvalue  = omnibus_pvalue
    )
}

run_ancova_2_sample <- function(x, response = c("dCq")) {

    response <- match.arg(response)

    test_formula <- dCq ~ Sample + ref_dCq

    test <- aov(test_formula, data = x)

    test_res <- broom::tidy(test) |>
        rename(
            Term       = term,
            Df         = df,
            `Sum Sq`   = sumsq,
            `Mean Sq`  = meansq,
            `F-value`  = statistic,
            `p-value`  = p.value
        ) |>
        mutate(Term = str_replace(Term, "ref_dCq", "Covariate (Reference dCq)"))

    test_label <- "ANCOVA"

    method <- glue("Two-sided ANCOVA on {format_response(response)} with reference sample {format_response(response)} covariate adjustment.")

    list(
        test_res        = test_res,
        test_label      = test_method,
        method          = method
    )
}

# Mixed Effect Model -----------------------------------------------------------
# For dCq with biological replicates as random effect

run_mixed_effect <- function(x,
                             response = c("dCq"),
                             comparison = c("pairwise", "trt.vs.ctrl"),
                             equal.var  = TRUE) {

    response <- match.arg(response)
    comparison <- match.arg(comparison)
    stopifnot(is.logical(equal.var), length(equal.var) == 1L, !is.na(equal.var))

    if (equal.var) {
        # Equal variance: use lmer with random intercept for Replicate
        test_formula <- dCq ~ Sample + (1 | Replicate)
        omnibus <- lmerTest::lmer(test_formula, data = x)
        omnibus_method <- "Mixed-effect model"

        # Extract F-test for Sample fixed effect
        anova_res  <- anova(omnibus)
        omnibus_pvalue <- anova_res["Sample", "Pr(>F)"]

        omnibus_res <- broom::tidy(anova_res) |>
            mutate(NumDF = round(NumDF, 3),
                   DenDF = round(DenDF, 3)) |>
            rename(
                Term             = term,
                `Sum Sq`         = sumsq,
                `Mean Sq`        = meansq,
                `Numerator Df`   = NumDF,
                `Denominator Df` = DenDF,
                `F-value`        = statistic,
                `p-value`        = p.value
            )

        omnibus_label <- glue(
            "Mixed-effect model: Sample F({numDf},{denDf}) = {sample_f}, p = {prettyNum(signif(omnibus_pvalue, 2))}",
            sample_f  = round(omnibus_res$`F-value`[1], 2),
            numDf       = round(omnibus_res$`Numerator Df`[1], 1),
            denDf       = round(omnibus_res$`Denominator Df`[1], 1)
        )

        # Post-hoc via emmeans
        emm <- emmeans(omnibus, ~Sample)

        if (comparison == "pairwise") {
            post_hoc        <- contrast(emm, method = "pairwise", adjust = "tukey")
            post_hoc_method <- "Tukey HSD"
        } else {
            set.seed(123)
            post_hoc        <- contrast(emm, method = "trt.vs.ctrl", adjust = "mvt")
            post_hoc_method <- "Dunnett"
        }

    } else {
        # Unequal variance: use nlme::lme with varIdent
        test_formula <- dCq ~ Sample
        omnibus <- nlme::lme(
            test_formula,
            random  = ~ 1 | Replicate,
            weights = nlme::varIdent(form = ~ 1 | Sample),
            data    = x
        )
        omnibus_method <- "Mixed-effect model (unequal variances)"

        # Extract F-test
        anova_res      <- anova(omnibus)
        omnibus_pvalue <- anova_res["Sample", "p-value"]

        omnibus_res <- as_tibble(anova_res, rownames = "Term") |>
            mutate(numDF = round(numDF, 3),
                   denDF = round(denDF, 3)) |>
            rename(
                `Numerator Df`   = numDF,
                `Denominator Df` = denDF
            )

        omnibus_label <- glue(
            "Mixed-effect model (unequal var): Sample F({numDf}, {denDf}) = {sample_f}, p = {prettyNum(signif(omnibus_pvalue, 2))}",
            sample_f = round(omnibus_res$`F-value`[omnibus_res$Term == "Sample"], 2),
            numDf    = round(omnibus_res$`Numerator Df`[omnibus_res$Term == "Sample"], 2),
            denDf    = round(omnibus_res$`Denominator Df`[omnibus_res$Term == "Sample"], 2)
        )

        # Post-hoc via emmeans with Satterthwaite df
        emm <- emmeans(omnibus, ~Sample)
        set.seed(123)

        if (comparison == "pairwise") {
            post_hoc        <- contrast(emm, method = "pairwise", adjust = "mvt")
            post_hoc_method <- "Dunnett T3"  # T3 accounts for unequal variances
            # TODO: this is wrong
        } else {
            post_hoc        <- contrast(emm, method = "trt.vs.ctrl", adjust = "mvt")
            post_hoc_method <- "Dunnett adjusted for unequal variances"  # is there a better name?

        }
    }

    # Format post-hoc results
    post_hoc_res <- post_hoc |>
        broom::tidy() |>
        format_emmeans() |>
        rename(
            Term = term,
            Df             = df,
            `Adj. p-value` = adj.p.value
        ) |>
        select(-null.value, -estimate, -std.error) |>
        mutate(Significance = to_significance(`Adj. p-value`))
    
        # rename the statistic column based on the type of post-hoc test
        if(post_hoc_method %in% c("Tukey HSD")) {
            post_hoc_res <- post_hoc_res |>
                rename(`q-value` = statistic)
        } else if (post_hoc_method %in% c("Dunnett", "Dunnett T3", "Dunnett adjusted for unequal variances")) {
            post_hoc_res <- post_hoc_res |>
                rename(`t-value` = statistic)
        } else {
            # throw an error saing that the method is unregonized
            stop("Unexpected post-hoc method")
        }

    method <- glue("Two-sided {omnibus_method} on {format_response(response)} with Replicate as random effect intercept, followed by {post_hoc_method} post-hoc test.")

    list(
        test_res        = post_hoc_res,
        test_label      = post_hoc_method,
        method          = method,
        omnibus_res     = omnibus_res,
        omnibus_label   = omnibus_label,
        omnibus_pvalue  = omnibus_pvalue
    )
}

# One-way ANOVA ----------------------------------------------------------------
# For ddCq or exp_ddCq (independent samples)

run_anova <- function(x,
                      comparison = c("pairwise", "trt.vs.ctrl"),
                      response   = c("ddCq", "exp_ddCq")) {

    comparison <- match.arg(comparison)
    response   <- match.arg(response)

    test_formula <- reformulate("Sample", response = response)

    # Omnibus: one-way ANOVA (assumes equal variances)
    omnibus        <- aov(test_formula, data = x)
    omnibus_method <- "One-way ANOVA"
    omnibus_pvalue <- summary(omnibus)[[1]]["Sample", "Pr(>F)"]

    omnibus_res <- broom::tidy(omnibus) |>
        rename(
            Term       = term,
            Df         = df,
            `Sum Sq`   = sumsq,
            `Mean Sq`  = meansq,
            `F-value`  = statistic,
            `p-value`  = p.value
        )

    omnibus_label <- glue(
        "One-way ANOVA: F({sample_df},{residuals_df}) = {sample_f_value}, p = {prettyNum(signif(omnibus_pvalue, 2))}",
        sample_df      = omnibus_res$Df[1],
        residuals_df   = omnibus_res$Df[2],
        sample_f_value = round(omnibus_res$`F-value`[1], 2)
    )
    
    emm <- emmeans(omnibus, ~Sample)
    
    if (comparison == "pairwise") {
        post_hoc <- contrast(emm, method = "pairwise", adjust = "tukey")
        post_hoc_res <- post_hoc |>
            broom::tidy() |>
            format_emmeans() |>
            rename(
                Term           = term,
                Df             = df,
                `q-value`      = statistic,
                `Adj. p-value` = adj.p.value
            ) |>
            mutate(Significance = to_significance(`Adj. p-value`)) |>
            select(-null.value, -estimate, -std.error)
        
        post_hoc_method <- "Tukey HSD"
        
    } else {
        set.seed(123)
        post_hoc <- contrast(emm, method = "trt.vs.ctrl", adjust = "mvt")
        post_hoc_res <- post_hoc |>
            broom::tidy() |>
            format_emmeans() |>
            rename(
                Term           = term,
                Df             = df,
                `t-value`      = statistic,
                `Adj. p-value` = adj.p.value
            ) |>
            mutate(Significance = to_significance(`Adj. p-value`)) |>
            select(-null.value, -estimate, -std.error)
        
        post_hoc_method <- "Dunnett"
    }

    method <- glue("Two-sided one-way ANOVA on {format_response(response)}, followed by {post_hoc_method} post-hoc test.")

    list(
        test_res        = post_hoc_res,
        test_label      = post_hoc_method,
        method          = method,
        omnibus_res     = omnibus_res,
        omnibus_label   = omnibus_label,
        omnibus_pvalue  = omnibus_pvalue
    )
}

# Kruskal-Wallis ---------------------------------------------------------------
# Non-parametric alternative to one-way ANOVA

run_kruskal <- function(x,
                        comparison = c("pairwise", "trt.vs.ctrl"),
                        response   = c("ddCq", "exp_ddCq"),
                        p.adjust   = c("BH", "holm", "none")) {

    comparison <- match.arg(comparison)
    response   <- match.arg(response)
    
    if (comparison == "pairwise") { # not requiret for many vs one dunn test (single-stpe adjustment)
        p.adjust   <- match.arg(p.adjust)
    }

    test_formula <- reformulate("Sample", response = response)

    # Omnibus: Kruskal-Wallis
    omnibus        <- kruskal.test(test_formula, data = x)
    omnibus_method <- "Kruskal-Wallis"
    omnibus_pvalue <- omnibus$p.value

    omnibus_res <- broom::tidy(omnibus) |>
        rename(
            `Chi-sq`  = statistic,
            Df        = parameter,
            `p-value` = p.value
        ) |>
        mutate(Term = "Sample") |> 
        select(Term, `Chi-sq`, Df, `p-value`)

    omnibus_label <- glue(
        "Kruskal-Wallis: χ²({df}) = {chi_sq}, p = {prettyNum(signif(omnibus_pvalue, 2))}",
        chi_sq = round(omnibus_res$`Chi-sq`, 2),
        df     = omnibus_res$Df
    )

    # Post-hoc: Dunn's test
    if (comparison == "pairwise") {
        post_hoc <- PMCMRplus::kwAllPairsDunnTest(
            test_formula, data = x, p.adjust.method = p.adjust
        )
        
        post_hoc_method <- ifelse(
            p.adjust != "none",
            glue("Dunn's ({str_replace(p.adjust, 'holm', 'Holm')} adjusted)"),
            glue("Dunn's (unadjusted)")
        )
        
    } else {
        post_hoc <- PMCMRplus::kwManyOneDunnTest(
            test_formula, data = x, p.adjust.method = "single-step"
        )
        post_hoc_method <- glue("Dunn's")
    }

    post_hoc_res <- format_pmcmr(post_hoc) |>
        mutate(Significance = to_significance(p.value),
               Term         = "Sample") |>
        rename(
            `z-value` = statistic,
        ) |>
        select(Term, group1, group2, `z-value`, p.value, Significance)
    
    
    if (p.adjust != "none") {
        post_hoc_res <- post_hoc_res |>
            rename(`Adj. p-value` = p.value)
    } else {
        post_hoc_res <- post_hoc_res |>
            rename(`p-value` = p.value)
    }
    
    method <- glue("Two-sided Kruskal-Wallis test on {format_response(response)}, followed by {post_hoc_method}.",
                   # add `post-hoc test` before the adjustment method
                   post_hoc_method = str_replace(post_hoc_method, "Dunn's", "Dunn's post-hoc test") 
    )

    list(
        test_res        = post_hoc_res,
        test_label      = post_hoc_method,
        method          = method,
        omnibus_res     = omnibus_res,
        omnibus_label   = omnibus_label,
        omnibus_pvalue  = omnibus_pvalue
    )
}

# Pairwise t-test --------------------------------------------------------------
# For ddCq or exp_ddCq (independent samples, no omnibus)

run_pairwise_ttest <- function(x,
                               response   = c("ddCq", "exp_ddCq"),
                               comparison = c("pairwise", "trt.vs.ctrl"),
                               equal.var  = TRUE,
                               p.adjust   = c("BH", "holm", "none")) {

    response <- match.arg(response)
    comparison <- match.arg(comparison)
    p.adjust <- match.arg(p.adjust)
    stopifnot(is.logical(equal.var), length(equal.var) == 1L, !is.na(equal.var))
    
    test_formula <- reformulate("Sample", response = response)
    
    if (comparison == "trt.vs.ctrl") {
        reference_sample <- x$Sample |> levels() |> head(n=1)
    } else {
        reference_sample <- NULL # pairwise comparisons
    }
        
    test <- rstatix::pairwise_t_test(
        data = x,
        formula = test_formula,
        var.equal = equal.var,
        ref.group = reference_sample,
        pool.sd = FALSE, # when pooling SD it becomes a Fisher's LSD test
        p.adjust.method = p.adjust
    )
    
    test_res <- test |>
        select(Term = ".y.", group1, group2, Df = df, `t-value` = statistic,
               `p-value` = p, `Adj. p-value` = p.adj) |>
        mutate(Significance = to_significance(`Adj. p-value`))
    
    if (p.adjust == "none") {
        # drop the adj column in case of no adjustment
        test_res <- select(test_res, -`Adj. p-value`) 
    } 
    
    adjust_label <- ifelse(
        p.adjust != "none",
        glue("({str_replace(p.adjust, 'holm', 'Holm')} adjusted)"),
        "(unadjusted)"
    )
    
    base_test <- case_when(
        comparison == "pairwise" & equal.var  ~ "pairwise t-test",
        comparison == "pairwise" & !equal.var ~ "pairwise Welch's t-test",
        comparison == "trt.vs.ctrl" & equal.var  ~ "many-to-one t-test",
        comparison == "trt.vs.ctrl" & !equal.var ~ "many-to-one one-sample t-test"
    )
    
    test_label <- glue("{base_test} {adjust_label}")
    
    
    method <- glue("Two-sided {test_label} on {format_response(response)}.")

    list(
        test_res   = test_res,
        test_label = test_label,
        method     = method
    )
}


# t-test (2 samples) --------------------------------------------------------------
# For ddCq or exp_ddCq (independent samples, no omnibus)
run_ttest <- function(x,
                      response  = c("ddCq", "exp_ddCq"),
                      equal.var = TRUE) {
    
    response <- match.arg(response)
    stopifnot(is.logical(equal.var), length(equal.var) == 1L, !is.na(equal.var))
    
    test_formula <- reformulate("Sample", response = response)
    
    
    test <- rstatix::t_test(
        data = x,
        formula = test_formula,
        var.equal = equal.var
    )
    
    test_res <- test |>
        select(Term = ".y.", group1, group2, Df = df, `t-value` = statistic,
               `p-value` = p) |>
        mutate(Significance = to_significance(`p-value`))
    
    test_label <- ifelse(
        equal.var,
        "t-test",
        "one-sample t-test"
    )
    
    method <- glue("Two-sided {test_label} on {format_response(response)}.")
    
    list(
        test_res   = test_res,
        test_label = test_label,
        method     = method
    )
}

# Pairwise paired t-test -------------------------------------------------------
# For dCq with biological replicates (paired by Replicate)

run_pairwise_paired_ttest <- function(x,
                                      response = c("dCq"),
                                      comparison = c("pairwise", "trt.vs.ctrl"),
                                      p.adjust = c("BH", "holm", "none")) {

    response   <- match.arg(response)
    comparison <- match.arg(comparison)
    p.adjust   <- match.arg(p.adjust)

    test_formula <- reformulate("Sample", response = response)
    
    if (comparison == "trt.vs.ctrl") {
        reference_sample <- x$Sample |> levels() |> head(n=1)
    } else {
        reference_sample <- NULL # pairwise comparisons
    }
    
    test <- rstatix::pairwise_t_test(
        data      = x |> arrange(Replicate), # ensure proper pairing
        formula   = test_formula,
        ref.group = reference_sample,
        paired    = T,
        p.adjust.method = p.adjust
    )
    
    test_res <- test |>
        select(Term = ".y.", group1, group2, Df = df, `t-value` = statistic,
               `p-value` = p, `Adj. p-value` = p.adj) |>
        mutate(Significance = to_significance(`Adj. p-value`))
    
    if (p.adjust == "none") {
        # drop the adj column in case of no adjustment
        test_res <- select(test_res, -`Adj. p-value`) 
    } 
    
    adjust_label <- ifelse(
        p.adjust != "none",
        glue("({str_replace(p.adjust, 'holm', 'Holm')} adjusted)"),
        "(unadjusted)"
    )
    
    base_test <- ifelse(comparison == "pairwise",
        "pairwise paired t-test",
        "many-to-one paired t-test"
    )
    
    test_label <- glue("{base_test} {adjust_label}")
    
    method <- glue("Two-sided {test_label} on {format_response(response)}.")
    
    list(
        test_res   = test_res,
        test_label = test_label,
        method     = method
    )
}


# Pairwise paired t-test -------------------------------------------------------
# For dCq with biological replicates (paired by Replicate)

run_paired_ttest <- function(x,
                             response = c("dCq")) {
    
    response   <- match.arg(response)
    test_formula <- reformulate("Sample", response = response)
    
    
    test <- rstatix::t_test(
        data      = x |> arrange(Replicate), # ensure proper pairing
        formula   = test_formula,
        paired    = T,
    )
    
    test_res <- test |>
        select(Term = ".y.", group1, group2, Df = df, `t-value` = statistic,
               `p-value` = p) |>
        mutate(Significance = to_significance(`p-value`))
    
    test_label <- "paired t-test"
    
    method <- glue("Two-sided paired t-test on {format_response(response)}.")
    
    list(
        test_res   = test_res,
        test_label = test_label,
        method     = method
    )
}


# Pairwise Wilcoxon signed-rank test -------------------------------------------
# Non-parametric paired test for dCq with replicates
run_pairwise_wilcoxon <- function(x,
                                      response = c("dCq"),
                                      comparison = c("pairwise", "trt.vs.ctrl"),
                                      p.adjust = c("BH", "holm", "none")) {
    
    response   <- match.arg(response)
    comparison <- match.arg(comparison)
    p.adjust   <- match.arg(p.adjust)
    
    test_formula <- reformulate("Sample", response = response)
    
    if (comparison == "trt.vs.ctrl") {
        reference_sample <- x$Sample |> levels() |> head(n=1)
    } else {
        reference_sample <- NULL # pairwise comparisons
    }
    
    test <- rstatix::pairwise_wilcox_test(
        data      = x |> arrange(Replicate), # ensure proper pairing
        formula   = test_formula,
        ref.group = reference_sample,
        paired    = T,
        p.adjust.method = p.adjust
    )
    
    test_res <- test |>
        select(Term = ".y.", group1, group2, `V-value` = statistic,
               `p-value` = p, `Adj. p-value` = p.adj) |>
        mutate(Significance = to_significance(`Adj. p-value`))
    
    if (p.adjust == "none") {
        # drop the adj column in case of no adjustment
        test_res <- select(test_res, -`Adj. p-value`) 
    } 
    
    adjust_label <- ifelse(
        p.adjust != "none",
        glue("({str_replace(p.adjust, 'holm', 'Holm')} adjusted)"),
        "(unadjusted)"
    )
    
    base_test <- ifelse(comparison == "pairwise",
                        "pairwise Wilcoxon signed-rank test",
                        "many-to-one Wilcoxon signed-rank test"
    )
    
    test_label <- glue("{base_test} {adjust_label}")
    
    method <- glue("Two-sided {test_label} on {format_response(response)}.")
    
    list(
        test_res   = test_res,
        test_label = test_label,
        method     = method
    )
}

# Wilcoxon signed-rank test -------------------------------------------
# Non-parametric paired test for dCq with replicates
run_wilcoxon <- function(x, response = c("dCq")) {
    
    response   <- match.arg(response)

    test_formula <- reformulate("Sample", response = response)
    
    test <- rstatix::wilcox_test(
        data      = x |> arrange(Replicate), # ensure proper pairing
        formula   = test_formula,
        paired    = T,
    )
    
    test_res <- test |>
        select(Term = ".y.", group1, group2, `V-value` = statistic,
               `p-value` = p) |>
        mutate(Significance = to_significance(`p-value`))

    test_label <- "Wilcoxon signed-rank test"
    
    method <- glue("Two-sided Wilcoxon signed-rank test on {format_response(response)}.")
    
    list(
        test_res   = test_res,
        test_label = test_label,
        method     = method
    )
}

# Pairwise Mann-Whitney U test -------------------------------------------------
# Non-parametric unpaired test for ddCq

run_pairwise_mann_whitney <- function(x,
                                      response = c("ddCq", "exp_ddCq"),
                                      comparison = c("pairwise", "trt.vs.ctrl"),
                                      p.adjust = c("BH", "holm", "none")) {

    response   <- match.arg(response)
    comparison <- match.arg(comparison)
    p.adjust   <- match.arg(p.adjust)

    test_formula <- reformulate("Sample", response = response)
    
    if (comparison == "trt.vs.ctrl") {
        reference_sample <- x$Sample |> levels() |> head(n=1)
    } else {
        reference_sample <- NULL # pairwise comparisons
    }
    
    test <- rstatix::pairwise_wilcox_test(
        data      = x,
        formula   = test_formula,
        ref.group = reference_sample,
        paired    = F,
        p.adjust.method = p.adjust
    )

    test_res <- test |>
        select(Term = ".y.", group1, group2, `U-value` = statistic,
               `p-value` = p, `Adj. p-value` = p.adj) |>
        mutate(Significance = to_significance(`Adj. p-value`))
    
    if (p.adjust == "none") {
        # drop the adj column in case of no adjustment
        test_res <- select(test_res, -`Adj. p-value`) 
    } 
    
    adjust_label <- ifelse(
        p.adjust != "none",
        glue("({str_replace(p.adjust, 'holm', 'Holm')} adjusted)"),
        "(unadjusted)"
    )
    
    base_test <- ifelse(comparison == "pairwise",
                        "pairwise Wilcoxon-Mann–Whitney test",
                        "many-to-one Wilcoxon-Mann–Whitney test"
    )
    
    test_label <- glue("{base_test} {adjust_label}")
    
    method <- glue("Two-sided {test_label} on {format_response(response)}.")
    
    list(
        test_res   = test_res,
        test_label = test_label,
        method     = method
    )
}

# Mann-Whitney U test -------------------------------------------------
# Non-parametric unpaired test for ddCq

run_mann_whitney <- function(x, response = c("ddCq", "exp_ddCq")) {
    
    response   <- match.arg(response)

    test_formula <- reformulate("Sample", response = response)
    
    test <- rstatix::wilcox_test(
        data      = x,
        formula   = test_formula,
        paired    = F
    )
    
    test_res <- test |>
        select(Term = ".y.", group1, group2, `U-value` = statistic,
               `p-value` = p) |>
        mutate(Significance = to_significance(`p-value`))
    
    test_label <- "Wilcoxon-Mann–Whitney test"
    
    method <- glue("Two-sided Wilcoxon-Mann–Whitney test on {format_response(response)}.")
    
    list(
        test_res   = test_res,
        test_label = test_label,
        method     = method
    )
}


# test the functions ---------
x <- read_csv("data/simulated_qPCR_data.csv") |>
    mutate(Sample = factor(Sample))
run_ancova(x, comparison = "trt.vs.ctrl")
# run_mixed_effect(x, comparison = "trt.vs.ctrl", equal.var = T)
# run_anova(x, comparison = "pairwise", response = "ddCq")
# run_kruskal(x, comparison = "pairwise", response = "exp_ddCq", p.adjust = "holm")
# run_pairwise_ttest(x, comparison = "trt.vs.ctrl", response = "ddCq", p.adjust = "BH")
# run_pairwise_wilcoxon(x, comparison = "trt.vs.ctrl", response = "dCq", p.adjust = "BH")
# run_pairwise_mann_whitney(x, comparison = "trt.vs.ctrl", response = "ddCq", p.adjust = "BH")

# TODO:
# protect post-hoc test
# when only 2 samples available, run it with adjustment = "none"
