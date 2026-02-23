# qPCR Statistical Tests Reference

Overview of all supported statistical tests in the Shiny qPCR app for the different response metrics.

> [!NOTE]
> Non-parametric tests additionally require ≥ 5 biological replicates.
> The sample count thresholds refer to the number of *finite* samples for parametric tests and total samples for non-parametric tests.

---

## 1 — ΔCq (dCq) tests

dCq tests operate on **non-independent (paired) data** (different samples' data points are linked by the biological replicate), and all the statistical tests try to correct for that in different ways.

- The **ANCOVA** uses the **reference sample's dCq** as a covariate to adjust for variability between the different runs.
- **Mixed Effect Models (MEM)** and **paired t-tests** correct for the average dCq value for each run. MEMs first correct for the average of each run and then performs pairwise comparisons on estimated marginal means. Repeated paired t-tests instead correct for the average dCq of each pair of comparison. A secondary difference is that repeated paired t-tests only use complete observations while MEMs can handle missing values. 

### ≥ 3 samples (parametric)

| Test | Recomended | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Control | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment |
|---|:---:|---|---|---|---|---|---|
| **ANCOVA** | ✅ | `aov(dCq ~ Sample + ref_dCq)` | Tukey HSD | Dunnett | `stats::aov` | `emmeans::contrast` | Tukey / Multivariate t distribution |
| **Mixed Effect Model** (equal var) | ✅ | `lmerTest::lmer(dCq ~ Sample + (1 | Replicate))` | Tukey HSD | Dunnett | `lmerTest::lmer` | `emmeans::contrast` | Tukey / Multivariate t distribution |
| **Mixed Effect Model** (unequal var) | | `nlme::lme(dCq ~ Sample, random = ~1|Replicate, weights = varIdent)` | Dunnett T3 | Dunnett (uneq. var) | `nlme::lme` | `emmeans::contrast` | Multivariate t distribution |
| **Repeated paired t-test** | | — | — | — | — | `rstatix::pairwise_t_test(paired=T)` | BH / Holm / none |

### ≥ 3 samples (non-parametric, ≥ 5 bio reps)

| Test | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Control | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment |
|---|---|---|---|---|---|---|
| **Repeated Wilcoxon signed-rank** | — | — | — | — | `rstatix::pairwise_wilcox_test(paired=T)` | BH / Holm / none |

### = 2 samples (parametric)

In the presence of complete data, the paired t-test and Mixed Effect Model converge to the same result for 2 samples.

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

## 2 — ΔΔCq (ddCq) and 2⁻ΔΔCq (exp_ddCq) tests

These tests operate on **independent** (unpaired) data because they assume that the ΔΔCq normalisation already corrected for any batch variability in biological replicates. Both ΔΔCq and 2⁻ΔΔCq share the same set of tests.

It is recommended to always test in the log space (ΔCq/ΔΔCq) for qPCR data because variance is usually more homogeneous (homoscedastic). The only limitation is that parametric tests cannot handle undetected (`Inf`) values, while in the linear space (2⁻ΔΔCq) those data points assume a value of zero and can be therefore modelled. If testing on 2⁻ΔΔCq data, we recommend turning on the "Handle unequal variance" toggle.

### ≥ 3 samples (parametric)

| Test | Recomended | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Control | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment options |
|---|:---:|---|---|---|---|---|---|
| **One-way ANOVA** | ✅ | `stats::aov` | Tukey HSD | Dunnett | `stats::aov` | `emmeans::contrast` | Tukey / Multivariate t distribution |
| **Repeated t-test** (equal var) | | — | — | — | — | `rstatix::pairwise_t_test` | BH / Holm / none |
| **Repeated Welch's t-test** (unequal var) | | — | — | — | — | `rstatix::pairwise_t_test(var.equal=F)` | BH / Holm / none |

### ≥ 3 samples (non-parametric, ≥ 5 bio reps)

| Test | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Control | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment options |
|---|---|---|---|---|---|---|
| **Kruskal-Wallis** | `stats::kruskal.test` | Dunn's | Dunn's | `stats::kruskal.test` | `PMCMRplus::kwAllPairsDunnTest` / `PMCMRplus::kwManyOneDunnTest` | BH / Holm / none (pairwise) · single-step (many-to-one) |
| **Repeated Wilcoxon-Mann–Whitney** | — | — | — | — | `stats::wilcox.test` | BH / Holm / none |

### = 2 samples (parametric)

| Test | Default | Package::function | Notes |
|---|:---:|---|---|
| **t-test** (equal var) | ✅ | `rstatix::t_test` | Single comparison |
| **Welch's t-test** (unequal var) | | `rstatix::t_test(var.equal=F)` | Single comparison |

### = 2 samples (non-parametric, ≥ 5 bio reps)

| Test | Package::function | Notes |
|---|---|---|
| **Wilcoxon-Mann–Whitney** | `stats::wilcox.test` | Performs the exact test by default but may apply "continuity correction" in cases of samples with all the same value (e.g. all undetected) |

---

## 3 — Handling of undetected/undefined values 

- Undetected values are modelled as `Inf` values by this tool.

- Undefined values arise from the ΔΔCq calculation when the reference sample is undetected. For genes with this behaviour, statistical testing is allowed only for the dCq metric.

Different tests handle undetected Cq values (which propagate as `Inf` in dCq and ddCq) in different ways:

| Test category | Handling strategy | Consequence |
|---|---|---|
| **ANCOVA** (dCq) | Drops entire replicate runs where reference dCq is `Inf`; then drops individual `Inf` response values. | May lose substantial data if the reference sample has undetected values. Reduces effective sample size. |
| **Paired t-test** (dCq) | For each comparison, uses only complete observations. | E.g. for TreatA vs. Ctrl comparison, if TreatA R3 is undetected, then Ctrl R3 is also ignored while it is preserved when testing against other samples. |
| **Other parametric** **tests** (ΔCq / ΔΔCq) | Drops individual undetected/undefined dCq values | Only individual undetected points are excluded. |
| **Parametric tests** (2⁻ΔΔCq) | When exponentiating, undetected values acquire a value of 0. | All undetected data points are preserved. Replicates for which the reference samples is undetected remains as undefined and are discarded. |
| **Non-parametric** (all metrics) | Replaces `Inf` with `999` to preserve rank information. | Undetected values are kept as the highest rank; no data loss. |

> [!IMPORTANT]
> Based on the selected test, a warning is displayed for discarded samples/replicates.

---

## 4 — Recommendations

### Recommended tests

- **ANCOVA** is recommended when you have a clear reference/control sample (e.g. an untreated condition). It adjusts for inter-replicate variability across runs via the reference sample's dCq covariate.
- **Mixed Effect Model** is better suited when the reference sample is arbitrary across replicates (e.g. comparing expression between different patients or cell lines). It models replicate-level variance as a random intercept, making it more appropriate for designs without a fixed control.
- Other tests (repeated paired t-tests, etc.) are available for completeness but are generally not the first choice.

### When to test on 2⁻ΔΔCq

Testing on the linear scale (2⁻ΔΔCq) is recommended **only** when the biological question specifically concerns genes that become completely undetectable (i.e. turned off). In this scenario, undetected values acquire a value of 0 in the linear space and can be included in the analysis, whereas on the log scale they are `Inf` and must be dropped. For all other cases, testing on ΔCq or ΔΔCq (log scale) is preferred because variance is typically more homoscedastic.

### Variance assumptions

- On the **ΔCq / ΔΔCq** (log) scale, data is usually homoscedastic, so the default equal-variance assumption is generally appropriate. However, it is still worth inspecting the data.
- On the **2⁻ΔΔCq** (linear) scale, variance heterogeneity is much more common because qPCR measurement noise is multiplicative (proportional to expression level). The log scale (ΔCq/ΔΔCq) converts this multiplicative noise into additive noise by the property of logarithms, yielding homogeneous variance. If testing on 2⁻ΔΔCq, we recommend using a test that models unequal variances (e.g. Mixed Effect Model with unequal variance toggle, or Welch's t-test).

### Non-parametric tests

Non-parametric tests are provided for completeness but are generally not recommended for qPCR data. They require at least 5 biological replicates and tend to have lower statistical power compared to their parametric counterparts. They may still be useful when the sample size is sufficient.
