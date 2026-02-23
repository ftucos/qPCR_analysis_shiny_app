library(tidyverse)
library(glue)
library(lme4)
library(broom)
library(lmerTest)
library(emmeans)
library(nlme)
library(PMCMRplus)
library(rstatix)

# change defaults of add_signif to mark up to 3 stars
add_signif <- partial(add_significance,
  cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  symbols = c("***", "**", "*", "ns")
)

format_response <- function(response) {
    case_match(response,
        "dCq"      ~ "ΔCq",
        "ddCq"     ~ "ΔΔCq",
        "exp_ddCq" ~  "2^-ΔΔCq",
        .default = response
    )
}

# Helper: split contrast column to group1 and group2 for ggpubr compatibility
format_emmeans <- function(df, sample_sizes) {
    df |>
        mutate(contrast = ifelse(
            str_detect(contrast, "\\(.*\\) - \\(.*\\)"),
            str_remove_all(contrast, "(^\\(|\\)$)") |> str_split("\\) - \\("),
            str_split(contrast, " - ")
        )) |>
        rename(group = contrast) |>
        unnest_wider(group, names_sep = "") |>
        left_join(sample_sizes |> rename("n1" = "n"), by = c("group1" = "Sample")) |>
        left_join(sample_sizes |> rename("n2" = "n"), by = c("group2" = "Sample")) |>
        relocate(c("n1", "n2"), .after = "group2")
}

# Helper: format PMCMRplus matrix output to tidy format with group columns
format_pmcmr <- function(pmcmr_result, sample_sizes) {
    PMCMRplus::toTidy(pmcmr_result) |>
        rename(group1 = group1, group2 = group2) |>
        add_signif(p.col = "p.value", output.col = "Significance") |>
        left_join(sample_sizes |> rename("n1" = "n"), by = c("group1" = "Sample")) |>
        left_join(sample_sizes |> rename("n2" = "n"), by = c("group2" = "Sample")) |>
        relocate(c("n1", "n2"), .after = "group2")
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
        add_signif(p.col = "p-value", output.col = "Significance")
    }

# ANCOVA -----------------------------------------------------------------------
# For dCq with biological replicates, adjusting for reference sample variance

run_ancova <- function(x,
                       response = c("dCq"),
                       comparison = c("pairwise", "trt.vs.ctrl")) {

    response <- match.arg(response)
    comparison <- match.arg(comparison)

    test_formula <- dCq ~ Sample + ref_dCq
    
    sample_sizes <- x |>
      tidyr::drop_na(dCq, ref_dCq) |>
      count(Sample)

    # Omnibus: ANCOVA
    omnibus <- aov(test_formula, data = x)
    omnibus_method  <- "ANCOVA"
    omnibus_pvalue  <- summary(omnibus)[[1]]["Sample", "Pr(>F)"]

    omnibus_res <- broom::tidy(omnibus) |>
        rename(
            Term      = term,
            Df        = df,
            `Sum Sq`  = sumsq,
            `Mean Sq` = meansq,
            `F-value` = statistic,
            `p-value` = p.value
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
            format_emmeans(sample_sizes = sample_sizes) |>
            rename(
                Term           = term,
                Df             = df,
                `q-value`      = statistic,
                `Adj. p-value` = adj.p.value
            ) |>
            add_signif(p.col = "Adj. p-value", output.col = "Significance") |>
            select(-null.value, -estimate, -std.error)
        
        post_hoc_method <- "Tukey HSD"

    } else {
        set.seed(123)
        post_hoc <- contrast(emm, method = "trt.vs.ctrl", adjust = "mvt")
        post_hoc_res <- post_hoc |>
            broom::tidy() |>
            format_emmeans(sample_sizes = sample_sizes) |>
            rename(
                Term           = term,
                Df             = df,
                `t-value`      = statistic,
                `Adj. p-value` = adj.p.value
            ) |>
            add_signif(p.col = "Adj. p-value", output.col = "Significance") |>
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
    
    sample_sizes <- x |>
      tidyr::drop_na(dCq, ref_dCq) |>
      count(Sample)

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
        mutate(Term = str_replace(Term, "ref_dCq", "Covariate (Reference dCq)"),
               group1 = ifelse(Term != "Residuals", levels(x$Sample)[1], NA),
               group2 = ifelse(Term != "Residuals", levels(x$Sample)[2], NA)) |>
        left_join(sample_sizes |> rename("n1" = "n"), by = c("group1" = "Sample")) |>
        left_join(sample_sizes |> rename("n2" = "n"), by = c("group2" = "Sample")) |>
        add_signif(p.col = "p-value", output.col = "Significance") |>
        relocate(c("group1", "group2", "n1", "n2"), .after = Term)
    
    test_label <- "ANCOVA"

    method <- glue("Two-sided ANCOVA on {format_response(response)} with reference sample {format_response(response)} covariate adjustment.")

    list(
        test_res       = test_res,
        test_label     = test_label,
        method         = method
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
    
    sample_sizes <- x |>
      tidyr::drop_na(all_of(response)) |>
      count(Sample)

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
                Term           = term,
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
        
        rand_effect_res <- as.data.frame(VarCorr(omnibus)) |>
          mutate(grp = str_replace(grp, "Replicate", "Replicate (Intercept)"),
                 vcov  = signif(vcov, 3),
                 sdcor = signif(sdcor, 3)) |>
          # rename to SD and variance given is a random intercept only
          select(Term = grp, Variance = vcov, `St. Dev.` = sdcor) |>
          as_tibble()

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
        
        rand_effect_res <- VarCorr(omnibus)[1:2, 1:2] |> 
          as.data.frame() |>
          rownames_to_column("Term") |>
          mutate(Term = str_replace(Term, "\\(Intercept\\)", "Replicate (Intercept)"),
                 Variance  = signif(as.numeric(Variance), 3),
                 StdDev = signif(as.numeric(StdDev), 3)) |>
          # rename to SD and variance given is a random intercept only
          select(Term, Variance, `St. Dev.` = StdDev) |>
          as_tibble()

        # Post-hoc via emmeans with Satterthwaite df
        emm <- emmeans(omnibus, ~Sample)
        set.seed(123)

        if (comparison == "pairwise") {
            post_hoc        <- contrast(emm, method = "pairwise", adjust = "mvt")
            post_hoc_method <- "Dunnett T3"  # T3 accounts for unequal variances
            # TODO: double check this
        } else {
            post_hoc        <- contrast(emm, method = "trt.vs.ctrl", adjust = "mvt")
            post_hoc_method <- "Dunnett adjusted for unequal variances"  # is there a better name?

        }
    }

    # Format post-hoc results
    post_hoc_res <- post_hoc |>
        broom::tidy() |>
        format_emmeans(sample_sizes = sample_sizes) |>
        rename(
            Term           = term,
            Df             = df,
            `Adj. p-value` = adj.p.value
        ) |>
        select(-null.value, -estimate, -std.error) |>
        add_signif(p.col = "Adj. p-value", output.col = "Significance")
        
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
        extra_res       = rand_effect_res,
        extra_label     = "Random Effect:",
        extra_position  = "omnibus",
        omnibus_label   = omnibus_label,
        omnibus_pvalue  = omnibus_pvalue
    )
}

# Mixed Effect Model (2 samples) -----------------------------------------------
run_mixed_effect_2_sample <- function(x,
                             response = c("dCq"),
                             equal.var  = TRUE) {
  
  response <- match.arg(response)
  stopifnot(is.logical(equal.var), length(equal.var) == 1L, !is.na(equal.var))
  
  sample_sizes <- x |>
    tidyr::drop_na(all_of(response)) |>
    count(Sample)
  
  if (equal.var) {
    # Equal variance: use lmer with random intercept for Replicate
    test_formula <- dCq ~ Sample + (1 | Replicate)
    test <- lmerTest::lmer(test_formula, data = x)
    test_label <- "Mixed-effect model"
    
    # Extract F-test for Sample fixed effect

    test_res <- anova(test) |>
      broom::tidy() |>
      mutate(NumDF = round(NumDF, 3),
             DenDF = round(DenDF, 3)
             ) |>
      rename(
        Term             = term,
        `Sum Sq`         = sumsq,
        `Mean Sq`        = meansq,
        `Numerator Df`   = NumDF,
        `Denominator Df` = DenDF,
        `F-value`        = statistic,
        `p-value`        = p.value
      ) 
    
    rand_effect_res <- as.data.frame(VarCorr(test)) |>
      mutate(grp = str_replace(grp, "Replicate", "Replicate (Intercept)"),
             vcov  = signif(vcov, 3),
             sdcor = signif(sdcor, 3)) |>
      # rename to SD and variance given is a random intercept only
      select(Term = grp, Variance = vcov, `St. Dev.` = sdcor) |>
      as_tibble()
    
  } else {
    # Unequal variance: use nlme::lme with varIdent
    test_formula <- dCq ~ Sample
    test <- nlme::lme(
      test_formula,
      random  = ~ 1 | Replicate,
      weights = nlme::varIdent(form = ~ 1 | Sample),
      data    = x
    )
    test_label <- "Mixed-effect model (unequal variances)"
    
    # Extract F-test
    test_res <- anova(test) |>
      as_tibble(rownames = "Term") |>
      mutate(numDF = round(numDF, 3),
             denDF = round(denDF, 3)) |>
      rename(
        `Numerator Df`   = numDF,
        `Denominator Df` = denDF
      )
    
    rand_effect_res <- VarCorr(test)[1:2, 1:2] |> 
      as.data.frame() |>
      rownames_to_column("Term") |>
      mutate(Term = str_replace(Term, "\\(Intercept\\)", "Replicate (Intercept)"),
             Variance  = signif(as.numeric(Variance), 3),
             StdDev = signif(as.numeric(StdDev), 3)) |>
      # rename to SD and variance given is a random intercept only
      select(Term, Variance, `St. Dev.` = StdDev) |>
      as_tibble()
  }
  
  # add groups and number of
  test_res <- test_res |>
    mutate(group1 = levels(x$Sample)[1],
           group2 = levels(x$Sample)[2]) |>
   left_join(sample_sizes |> rename("n1" = "n"), by = c("group1" = "Sample")) |>
   left_join(sample_sizes |> rename("n2" = "n"), by = c("group2" = "Sample")) |>
   add_signif(p.col = "p-value", output.col = "Significance") |>
   relocate(c("group1", "group2", "n1", "n2"), .after = "Term")
  
  method <- glue("Two-sided {test_label} on {format_response(response)} with Replicate as random effect intercept.")
  
  list(
    test_res        = test_res,
    extra_res       = rand_effect_res,
    extra_label     = "Random Effect:",
    extra_position  = "comparisons",
    test_label      = test_label,
    method          = method
  )
}

# One-way ANOVA ----------------------------------------------------------------
# For ddCq or exp_ddCq (independent samples)

run_anova <- function(x,
                      comparison = c("pairwise", "trt.vs.ctrl"),
                      response   = c("ddCq", "exp_ddCq")) {

    comparison <- match.arg(comparison)
    response   <- match.arg(response)
    
    sample_sizes <- x |>
      tidyr::drop_na(all_of(response)) |>
      count(Sample)

    test_formula <- reformulate("Sample", response = response)

    # Omnibus: one-way ANOVA (assumes equal variances)
    omnibus        <- aov(test_formula, data = x)
    omnibus_method <- "One-way ANOVA"
    omnibus_pvalue <- summary(omnibus)[[1]]["Sample", "Pr(>F)"]

    omnibus_res <- broom::tidy(omnibus) |>
        rename(
            Term     = term,
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
            format_emmeans(sample_sizes = sample_sizes) |>
            rename(
                Term           = term,
                Df             = df,
                `q-value`      = statistic,
                `Adj. p-value` = adj.p.value
            ) |>
            add_signif(p.col = "Adj. p-value", output.col = "Significance") |>
            select(-null.value, -estimate, -std.error)
        
        post_hoc_method <- "Tukey HSD"
        
    } else {
        set.seed(123)
        post_hoc <- contrast(emm, method = "trt.vs.ctrl", adjust = "mvt")
        post_hoc_res <- post_hoc |>
            broom::tidy() |>
            format_emmeans(sample_sizes = sample_sizes) |>
            rename(
                Term           = term,
                Df             = df,
                `t-value`      = statistic,
                `Adj. p-value` = adj.p.value
            ) |>
            add_signif(p.col = "Adj. p-value", output.col = "Significance") |>
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
                        p_adjust_method = c("BH", "holm", "none")) {

    comparison <- match.arg(comparison)
    response   <- match.arg(response)
    
    sample_sizes <- x |>
      tidyr::drop_na(all_of(response)) |>
      count(Sample)
    
    if (comparison == "pairwise") { # not requiret for many vs one dunn test (single-stpe adjustment)
        p_adjust_method   <- match.arg(p_adjust_method)
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
            test_formula, data = x, p.adjust.method = p_adjust_method
        )
        
        post_hoc_method <- ifelse(
            p_adjust_method != "none",
            glue("Dunn's test ({str_replace(p_adjust_method, 'holm', 'Holm')} adjusted)"),
            "Dunn's test (unadjusted)"
        )
        
    } else {
        post_hoc <- PMCMRplus::kwManyOneDunnTest(
            test_formula, data = x, p.adjust.method = "single-step"
        )
        post_hoc_method <- "Dunn's test"
    }

    post_hoc_res <- post_hoc |>
        format_pmcmr(sample_sizes = sample_sizes) |>
        add_signif(p.col = "p.value", output.col = "Significance") |>
        mutate(Term = "Sample") |>
        rename(
            `z-value` = statistic,
        ) |>
        select(Term, group1, group2, n1, n2, `z-value`, p.value, Significance)
    
    
    if (p_adjust_method != "none") {
        post_hoc_res <- post_hoc_res |>
            rename(`Adj. p-value` = p.value)
    } else {
        post_hoc_res <- post_hoc_res |>
            rename(`p-value` = p.value)
    }
    
    method <- glue("Two-sided Kruskal-Wallis test on {format_response(response)}, followed by {post_hoc_method}.",
                   # add `post-hoc test` before the adjustment method
                   post_hoc_method = str_replace(post_hoc_method, "Dunn's test", "Dunn's post-hoc test") 
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

run_repeated_ttest <- function(x,
                               response   = c("ddCq", "exp_ddCq"),
                               comparison = c("pairwise", "trt.vs.ctrl"),
                               equal.var  = TRUE,
                               p_adjust_method   = c("BH", "holm", "none")) {

    response <- match.arg(response)
    comparison <- match.arg(comparison)
    p_adjust_method <- match.arg(p_adjust_method)
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
        p.adjust.method = p_adjust_method
    )
    
    test_res <- test |>
        mutate(Term = "Sample") |>
        select(Term, group1, group2, n1, n2, Df = df, `t-value` = statistic,
               `p-value` = p, `Adj. p-value` = p.adj) |>
        add_signif(p.col = "Adj. p-value", output.col = "Significance")
        
    if (p_adjust_method == "none") {
        # drop the adj column in case of no adjustment
        test_res <- select(test_res, -`Adj. p-value`) 
    } 
    
    adjust_label <- ifelse(
        p_adjust_method != "none",
        glue("({str_replace(p_adjust_method, 'holm', 'Holm')} adjusted)"),
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
        mutate(Term = "Sample") |>
        select(Term, group1, group2, n1, n2, Df = df, `t-value` = statistic,
               `p-value` = p) |>
        add_signif(p.col = "p-value", output.col = "Significance")
        
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

run_repeated_paired_ttest <- function(x,
                                      response = c("dCq"),
                                      comparison = c("pairwise", "trt.vs.ctrl"),
                                      p_adjust_method = c("BH", "holm", "none")) {

    response   <- match.arg(response)
    comparison <- match.arg(comparison)
    p_adjust_method   <- match.arg(p_adjust_method)

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
        p.adjust.method = p_adjust_method
    )
    
    test_res <- test |>
        mutate(Term = "Sample") |>
        select(Term, group1, group2, n1, n2, Df = df, `t-value` = statistic,
               `p-value` = p, `Adj. p-value` = p.adj) |>
        add_signif(p.col = "Adj. p-value", output.col = "Significance")
        
    if (p_adjust_method == "none") {
        # drop the adj column in case of no adjustment
        test_res <- select(test_res, -`Adj. p-value`) 
    } 
    
    adjust_label <- ifelse(
        p_adjust_method != "none",
        glue("({str_replace(p_adjust_method, 'holm', 'Holm')} adjusted)"),
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
        mutate(Term = "Sample") |>
        select(Term, group1, group2, n1, n2, Df = df, `t-value` = statistic,
               `p-value` = p) |>
        add_signif(p.col = "p-value", output.col = "Significance")
        
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
run_repeated_wilcoxon <- function(x,
                                      response = c("dCq"),
                                      comparison = c("pairwise", "trt.vs.ctrl"),
                                      p_adjust_method = c("BH", "holm", "none")) {
    
    response   <- match.arg(response)
    comparison <- match.arg(comparison)
    p_adjust_method   <- match.arg(p_adjust_method)
    
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
        p.adjust.method = p_adjust_method
    )
    
    test_res <- test |>
        mutate(Term = "Sample") |>
        select(Term, group1, group2, n1, n2, `V-value` = statistic,
               `p-value` = p, `Adj. p-value` = p.adj) |>
        add_signif(p.col = "Adj. p-value", output.col = "Significance")
        
    if (p_adjust_method == "none") {
        # drop the adj column in case of no adjustment
        test_res <- select(test_res, -`Adj. p-value`) 
    } 
    
    adjust_label <- ifelse(
        p_adjust_method != "none",
        glue("({str_replace(p_adjust_method, 'holm', 'Holm')} adjusted)"),
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
        mutate(Term = "Sample") |>
        select(Term, group1, group2, n1, n2, `V-value` = statistic,
               `p-value` = p) |>
        add_signif(p.col = "p-value", output.col = "Significance")
        
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

run_repeated_mann_whitney <- function(x,
                                      response = c("ddCq", "exp_ddCq"),
                                      comparison = c("pairwise", "trt.vs.ctrl"),
                                      p_adjust_method = c("BH", "holm", "none")) {

    response   <- match.arg(response)
    comparison <- match.arg(comparison)
    p_adjust_method   <- match.arg(p_adjust_method)

    sample_sizes <- x |>
      tidyr::drop_na(all_of(response)) |>
      count(Sample)
    
    test_formula <- reformulate("Sample", response = response)
    
    # Identify groups
    groups <- x$Sample |> droplevels() |> levels()
    if (is.null(groups)) groups <- unique(x$Sample) # fallback in case Sample is not a factor
    
    # Define comparisons
    comparisons_list <- list()
    
    if (comparison == "pairwise") {
        # All pairwise combinations
        comparisons_list <- combn(groups, 2, simplify = FALSE)
    } else {
        # Treatment vs Control (Reference)
        ref_group <- groups[1]
        other_groups <- groups[-1]
        comparisons_list <- lapply(other_groups, function(g) c(ref_group, g))
    }

    # Run tests manually
    test_res <- map(comparisons_list, function(sample_pair) {
        subset_data <- x |> filter(Sample %in% sample_pair)
        
        wilcox.test(test_formula, data = subset_data) |>
            broom::tidy() |>
            mutate(Term = "Sample",
                   group1 = sample_pair[[1]],
                   group2 = sample_pair[[2]],
                   method = str_extract(method, "(exact|with continuity correction)") |>
                                        str_to_sentence()
            ) |>
            select(Term, group1, group2, `U-value` = statistic, `p-value` = p.value, Details = method)
    }) |>
        bind_rows() |>
        left_join(sample_sizes |> rename("n1" = "n"), by = c("group1" = "Sample")) |>
        left_join(sample_sizes |> rename("n2" = "n"), by = c("group2" = "Sample")) |>
        relocate(c("n1", "n2"), .after = "group2")
    
    if (p_adjust_method == "none") {
        # drop the adj column in case of no adjustment
        test_res <- test_res |>
            add_signif(p.col = "p-value", output.col = "Significance")
    } else {
        test_res <- test_res |>
            mutate("Adj. p-value" = p.adjust(`p-value`, method = p_adjust_method)) |>
            add_signif(p.col = "Adj. p-value", output.col = "Significance")
    
    }
    
    adjust_label <- ifelse(
        p_adjust_method != "none",
        glue("({str_replace(p_adjust_method, 'holm', 'Holm')} adjusted)"),
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
    
    sample_sizes <- x |>
      tidyr::drop_na(all_of(response)) |>
      count(Sample)

    test_formula <- reformulate("Sample", response = response)
    
    # use stats::wilcox.test in place of rstatix::wilcox_test because it fails when 
    test <- wilcox.test(formula = test_formula, data = x)
    
    test_res <- test |>
        broom::tidy() |>
        mutate(Term = "Sample",
               group1 = x$Sample[[1]],
               group2 = x$Sample[[2]],
               ) |>
        add_signif(p.col = "p.value", output.col = "Significance") |>
        select(Term, group1, group2, `U-value` = statistic, `p-value` = p.value, Significance) |>
        left_join(sample_sizes |> rename("n1" = "n"), by = c("group1" = "Sample")) |>
        left_join(sample_sizes |> rename("n2" = "n"), by = c("group2" = "Sample")) |>
        relocate(c("n1", "n2"), .after = "group2")

    test_label <- case_match(test$method,
                        "Wilcoxon rank sum test with continuity correction" ~ "Wilcoxon-Mann–Whitney test with continuity correction",
                        "Wilcoxon rank sum exact test" ~ "Wilcoxon-Mann–Whitney test")
    
    method <- glue("Two-sided {test_label} on {format_response(response)}.")
    
    list(
        test_res   = test_res,
        test_label = test_label,
        method     = method
    )
}


# test the functions ---------
# x <- read_csv("data/simulated_qPCR_data.csv") |>
#     mutate(Sample = factor(Sample))
# run_ancova(x, comparison = "pairwise")
# run_mixed_effect(x, comparison = "trt.vs.ctrl", equal.var = T)
# run_anova(x, comparison = "pairwise", response = "ddCq")
# run_kruskal(x, comparison = "pairwise", response = "exp_ddCq", p_adjust_method = "holm")
# run_repeated_paired_ttest(x, comparison = "trt.vs.ctrl", response = "dCq", p_adjust_method = "BH")
# run_repeated_ttest(x, comparison = "trt.vs.ctrl", response = "ddCq", p_adjust_method = "BH")
# run_repeated_wilcoxon(x, comparison = "trt.vs.ctrl", response = "dCq", p_adjust_method = "BH")
# run_repeated_mann_whitney(x, comparison = "trt.vs.ctrl", response = "ddCq", p_adjust_method = "BH")
# x_2 <- x |> filter(Sample %in% c("Ctrl", "TrtA"))
# run_ancova_2_sample(x_2, response = "dCq")
# run_mixed_effect_2_sample(x_2, response = "dCq")
# run_paired_ttest(x_2, response = "dCq")
# run_ttest(x_2, response = "ddCq")
# run_wilcoxon(x_2, response = "dCq")
# run_mann_whitney(x_2, response = "ddCq")

# TODO:
# protect post-hoc test
