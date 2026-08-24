# qPCR Statistical Tests Reference

Overview of all supported statistical tests in the Shiny qPCR app for the different response metrics.

---

## 1 — $ΔCq$ tests

This app assumes that dCq tests operate on **non-independent (paired) data** (there is a batch effect to be accounted for between the different biological replicates), and all the statistical tests try to correct for that in different ways.

- The **ANCOVA** uses the **reference sample's dCq** as a covariate to adjust for variability between the different runs (this is conceptually equivalent to what you do with $ΔΔCq$ calculation, but you let the statistical method to apprirpately correct for that.
- **Mixed Effect Models (MEM)** and **paired t-tests** correct for the average dCq value for each run. MEMs first correct for the average of each run and then performs pairwise comparisons on estimated marginal means. Repeated paired t-tests instead correct for the average dCq of each pair of comparison. A secondary difference is that repeated paired t-tests only use complete observations while MEMs can handle missing values. 

### ≥ 3 samples (parametric)

| Test | Recomended | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Reference | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment |
|---|:---:|---|---|---|---|---|---|
| **ANCOVA** | ✅ | `aov(dCq ~ Sample + ref_dCq)` | Tukey HSD | Dunnett | `stats::aov` | `emmeans::contrast` | Tukey / Multivariate t distribution |
| **Mixed Effect Model** (equal var) | ✅ | `lmerTest::lmer(dCq ~ Sample + (1 | Replicate))` | Tukey HSD | Dunnett | `lmerTest::lmer` | `emmeans::contrast` | Tukey / Multivariate t distribution |
| **Mixed Effect Model** (unequal var) | | `nlme::lme(dCq ~ Sample, random = ~1|Replicate, weights = varIdent)` | Dunnett T3 | Dunnett (uneq. var) | `nlme::lme` | `emmeans::contrast` | Multivariate t distribution |
| **Repeated paired t-test** | | — | — | — | — | `rstatix::pairwise_t_test(paired=T)` | BH / Holm / none |

### ≥ 3 samples (non-parametric, ≥ 5 bio reps)

| Test | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Reference | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment |
|---|---|---|---|---|---|---|
| **Repeated Wilcoxon signed-rank** | — | — | — | — | `rstatix::pairwise_wilcox_test(paired=T)` | BH / Holm / none |

### = 2 samples (parametric)

In the presence of complete data, the paired t-test and Mixed Effect Model matematically converge to the same result for 2 samples.

| Test | Recomended | Package::function |
|---|:---:|---|
| **ANCOVA** (2 sample) | ✅ | `stats::aov` |
| **Mixed Effect Model** (2 sample, equal var) | ✅ | `lmerTest::lmer` |
| **Mixed Effect Model** (2 sample, unequal var) | | `nlme::lme` + `nlme::varIdent` |
| **Paired t-test** | | `rstatix::t_test(paired=T)` |

### = 2 samples (non-parametric, ≥ 5 bio reps)

| Test | Package::function |
|---|---|
| **Wilcoxon signed-rank** | `rstatix::wilcox_test(paired=T)` |

---

## 2 — $ΔΔCq$ and $2^{-ΔΔCq}$ tests

These tests operate on **independent** (unpaired) data because they assume that the ΔΔCq normalisation already corrected for any batch variability in biological replicates. Both $ΔΔCq$ and $2^{-ΔΔCq}$ share the same set of tests.

It is recommended to test in the log space ($ΔCq$/$ΔΔCq$) because qPCR variance is usually more homogeneous (homoscedastic) there. If testing on $2^{-ΔΔCq}$ data, we recommend turning on the "Handle unequal variance" toggle in tests that support it.

> [!NOTE]
> Non-parametric tests additionally require ≥ 5 biological replicates.

### ≥ 3 samples (parametric)

| Test | Recomended | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Reference | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment options |
|---|:---:|---|---|---|---|---|---|
| **One-way ANOVA** | ✅ | `stats::aov` | Tukey HSD | Dunnett | `stats::aov` | `emmeans::contrast` | Tukey / Multivariate t distribution |
| **Repeated t-test** (equal var) | | — | — | — | — | `rstatix::pairwise_t_test` | BH / Holm / none |
| **Repeated Welch's t-test** (unequal var) | | — | — | — | — | `rstatix::pairwise_t_test(var.equal=F)` | BH / Holm / none |

### ≥ 3 samples (non-parametric, ≥ 5 bio reps)

| Test | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Reference | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment options |
|---|---|---|---|---|---|---|
| **Kruskal-Wallis** | `stats::kruskal.test` | Dunn's | Dunn's | `stats::kruskal.test` | `PMCMRplus::kwAllPairsDunnTest` / `PMCMRplus::kwManyOneDunnTest` | BH / Holm / none (pairwise) · single-step (many-to-one) |
| **Repeated Wilcoxon-Mann–Whitney** | — | — | — | — | `stats::wilcox.test` | BH / Holm / none |

### = 2 samples (parametric)

| Test | Default | Package::function | Notes |
|---|:---:|---|---|
| **Welch's t-test** (unequal var) | ✅ | `rstatix::t_test(var.equal=F)` | Single comparison |
| **Student's t-test** (equal var) | | `rstatix::t_test(var.equal=T)` | Single comparison |

### = 2 samples (non-parametric, ≥ 5 bio reps)

| Test | Package::function | Notes |
|---|---|---|
| **Wilcoxon-Mann–Whitney** | `stats::wilcox.test` | Performs the exact test by default but may apply "continuity correction" in cases of samples with all the same value (e.g. all undetected) |

---

## 3 — Handling of undetected/censored values

- The default undetected replacement cycle is 40. If a detected Cq is greater than the configured replacement, the app raises the replacement to the first integer above that detected value.
- Within a technical-replicate group, undetected measurements are excluded whenever at least one technical replicate is detected. If every technical replicate is undetected, their replacement values are retained and labeled `Undetected` in both the Technical Replicates and Bio Rep Averages exports. 
`Cq_n` records the measurements actually used for the calculation.
- The app retains numeric replacements and censoring flags internally for Cq, ΔCq, ΔΔCq, and exponentiated results. Exports omit the internal `*_numeric` and `*_censored` columns and apply `>`/`<` labels directly to the displayed metric (for example `>40`, `>20.2`, or `<-20.2`). Summary tables report `dCq_undetected_n` and `ddCq_undetected_n`.
- Any sample can be selected as the reference, including one with an undetected target value, a missing target value, or a biological replicate excluded because its housekeeping gene was undetected. The app recommends a complete alternative but does not force the change. Affected replicate-level ΔΔCq values become `NA` and are omitted from ΔΔCq statistical tests and ANCOVA; unaffected replicates remain available.

> [!IMPORTANT]
> Replacing undetected Cq values with a fixed maximum-cycle value can bias expression estimates and the resulting statistical inference. See McCall et al. (2014), [*On non-detects in qPCR data*](https://doi.org/10.1093/bioinformatics/btu239).

---

## 4 — Recommendations

### Recommended tests

- **ANCOVA** is recommended when you have a clear reference/control sample (e.g. an untreated condition). It adjusts for inter-replicate variability across runs via the reference sample's dCq covariate.
- **Mixed Effect Model** is better suited when the reference sample is arbitrary across replicates (e.g. comparing expression between different patients or cell lines). It models replicate-level variance as a random intercept, making it more appropriate for designs without a fixed control.
- Other tests (repeated paired t-tests, etc.) are available for completeness but are generally not the first choice.

### Variance assumptions

- On the **ΔCq / ΔΔCq** (log) scale, data is usually homoscedastic, so the default equal-variance assumption is generally appropriate. However, it is still worth inspecting the data.
- On the **2⁻ΔΔCq** (linear) scale, variance heterogeneity is much more common because qPCR measurement noise is multiplicative (proportional to expression level). The log scale (ΔCq/ΔΔCq) converts this multiplicative noise into additive noise by the property of logarithms, yielding homogeneous variance. If testing on 2⁻ΔΔCq, we recommend using a test that models unequal variances (Welch's t-test).

### Non-parametric tests

Non-parametric tests are provided for completeness but are generally not recommended for qPCR data. They require at least 5 biological replicates and tend to have lower statistical power compared to their parametric counterparts. They may still be useful when the sample size is sufficient.
