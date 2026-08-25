# qPCR Statistical Tests Reference

Overview of all supported statistical tests in the Shiny qPCR app for the different response metrics.

---

## 1 — Tests on Log₂ normalized expression ($ΔCq$)

This app assumes that dCq tests operate on **non-independent (paired) data** (there is a batch effect to be accounted for between the different biological replicates), and all the statistical tests try to correct for that in different ways. Analysis in Cq space is preferred because Cq is the observed experimental value ([Yuan et al., 2006](https://doi.org/10.1186/1471-2105-7-85)).

- The **ANCOVA** uses the **reference sample's dCq** as a covariate to adjust for variability between runs. Its estimated sample effects correspond to $ΔΔCq$, allowing the model to derive the normalized effect directly ([Yuan et al., 2006](https://doi.org/10.1186/1471-2105-7-85)).
- **Mixed Effect Models (MEM)** and **paired t-tests** correct for the average dCq value for each run. MEMs first correct for the average of each run and then performs pairwise comparisons on estimated marginal means. Repeated paired t-tests instead correct for the average dCq of each pair of comparison. A secondary difference is that repeated paired t-tests only use complete observations while MEMs can handle missing values. 

### ≥ 3 samples (parametric)

| Test | Recommended | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Reference | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment |
|---|:---:|---|---|---|---|---|---|
| **ANCOVA** | ✅ | `aov(dCq ~ Sample + ref_dCq)` | Tukey HSD | Dunnett | `stats::aov` | `emmeans::contrast` | Tukey / Multivariate t distribution |
| **Mixed Effect Model** (equal var) | ✅ | `lmerTest::lmer(dCq ~ Sample + (1 | Replicate))` | Tukey HSD | Dunnett | `lmerTest::lmer` | `emmeans::contrast` | Tukey / Multivariate t distribution |
| **Mixed Effect Model** (unequal var) | | `nlme::lme(dCq ~ Sample, random = ~1|Replicate, weights = varIdent)` | Dunnett T3 | Dunnett (uneq. var) | `nlme::lme` | `emmeans::contrast` | Multivariate t distribution |
| **Repeated paired t-test** | | — | — | — | — | `rstatix::pairwise_t_test(paired=T)` | BH / Holm / none |

### ≥ 3 samples (non-parametric, ≥ 5 bio reps)

| Test | Recommended | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Reference | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment |
|---|:---:|---|---|---|---|---|---|
| **Repeated Wilcoxon signed-rank** | | — | — | — | — | `rstatix::pairwise_wilcox_test(paired=T)` | BH / Holm / none |

### = 2 samples (parametric)

In the presence of complete data, the paired t-test and Mixed Effect Model matematically converge to the same result for 2 samples.

| Test | Recommended | Package::function |
|---|:---:|---|
| **ANCOVA** (2 sample) | ✅ | `stats::aov` |
| **Mixed Effect Model** (2 sample, equal var) | ✅ | `lmerTest::lmer` |
| **Mixed Effect Model** (2 sample, unequal var) | | `nlme::lme` + `nlme::varIdent` |
| **Paired t-test** | | `rstatix::t_test(paired=T)` |

### = 2 samples (non-parametric, ≥ 5 bio reps)

| Test | Recommended | Package::function |
|---|:---:|---|
| **Wilcoxon signed-rank** | | `rstatix::wilcox_test(paired=T)` |

---

## 2 — Tests on (Log₂) fold change ($ΔΔCq$ and $2^{-ΔΔCq}$)

Testing on the **ΔCq scale is preferred**. Tests on $ΔΔCq$ or $2^{-ΔΔCq}$ are provided as a fallback for users who specifically want to analyze reference-normalized fold changes. These tests operate on **independent** (unpaired) data because they assume that the ΔΔCq normalization already corrected for batch variability among biological replicates.

By construction, every reference-sample value is fixed at $0$ on the $ΔΔCq$ scale and at $1$ on the $2^{-ΔΔCq}$ scale. The reference group therefore has zero variance. Welch's ANOVA cannot accommodate a group with zero variance and is not offered for these metrics. Pairwise Welch's t-tests remain usable: a comparison against the fixed reference group reduces mathematically to a one-sample t-test against $0$ for $ΔΔCq$, or against $1$ for $2^{-ΔΔCq}$. Consequently, **repeated Welch's t-test is the default when there are three or more samples**, and **Welch's t-test is the default when there are two samples**.

The $ΔCq$/$ΔΔCq$ scales are logarithmic and generally have more homogeneous variance than the exponentiated $2^{-ΔΔCq}$ scale. Statistical analysis should be performed in Cq space, with the linear transformation used for reporting expression levels ([Taylor et al., 2019](https://doi.org/10.1016/j.tibtech.2018.12.002)).

> [!NOTE]
> Non-parametric tests additionally require ≥ 5 biological replicates.

### ≥ 3 samples (parametric)

| Test | Recommended | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Reference | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment options | Notes |
|---|:---:|---|---|---|---|---|---|---|
| **Repeated Welch's t-test** (unequal var; default) | ✅ | — | — | — | — | `rstatix::pairwise_t_test(var.equal=F)` | BH / Holm / none | Comparisons vs. reference reduce to a one-sample t-test against 0 (ΔΔCq) or 1 (2⁻ΔΔCq) |
| **Repeated t-test** (equal var) | | — | — | — | — | `rstatix::pairwise_t_test` | BH / Holm / none ||
| **One-way ANOVA** |  | `stats::aov` | Tukey HSD | Dunnett | `stats::aov` | `emmeans::contrast` | Tukey / Multivariate t distribution ||

### ≥ 3 samples (non-parametric, ≥ 5 bio reps)

| Test | Recommended | Omnibus | Post-hoc: Pairwise | Post-hoc: All vs Reference | Package::function (omnibus) | Package::function (post-hoc) | p-adjustment options |
|---|:---:|---|---|---|---|---|---|
| **Kruskal-Wallis** || `stats::kruskal.test` | Dunn's | Dunn's | `stats::kruskal.test` | `PMCMRplus::kwAllPairsDunnTest` / `PMCMRplus::kwManyOneDunnTest` | BH / Holm / none (pairwise) · single-step (many-to-one) |
| **Repeated Wilcoxon-Mann–Whitney** |  | — | — | — | `stats::wilcox.test` | BH / Holm / none ||

### = 2 samples (parametric)

| Test | Recommended | Package::function | Notes |
|---|:---:|---|---|
| **Welch's t-test** (unequal var; default) | ✅ | `rstatix::t_test(var.equal=F)` | Reduces to a one-sample t-test against 0 (ΔΔCq) or 1 (2⁻ΔΔCq) |
| **Student's t-test** (equal var) | | `rstatix::t_test(var.equal=T)` |  |

### = 2 samples (non-parametric, ≥ 5 bio reps)

| Test | Recommended | Package::function | Notes |
|---|:---:|---|---|
| **Wilcoxon-Mann–Whitney** | `stats::wilcox.test` | Performs the exact test by default but may apply "continuity correction" in cases of samples with all the same value (e.g. all undetected) |

---

## 3 — Handling of undetected/censored values

- The default undetected replacement cycle is 40. If a detected Cq is greater than the configured replacement, the app raises the replacement to the first integer above that detected value.
- Within a technical-replicate group, undetected measurements are excluded whenever at least one technical replicate is detected. If every technical replicate is undetected, their replacement values are retained and labeled `Undetected` in both the Technical Replicates and Bio Rep Averages exports. 
`Cq_n` records the number of measurements actually used for the calculation.
- The app retains numeric replacements and censoring flags internally for Cq, ΔCq, ΔΔCq, and exponentiated results. Exports omit the internal `*_numeric` and `*_censored` columns and apply `>`/`<` labels directly to the displayed metric (e.g. `>40`, `>20.2`, or `<-20.2`).
- For ΔΔCq and ANCOVA, we recommend selecting reference samples with complete observations. You can still choose samples with missing or undetected replicates, the app will suggest a complete alternative but will not block your choice. Any affected replicate-level ΔΔCq values will be set to NA and excluded from statistical tests, while unaffected replicates remain available.

> [!IMPORTANT]
> Replacing undetected Cq values with a fixed maximum-cycle value can bias expression estimates and the resulting statistical inference. See McCall et al. (2014), [*On non-detects in qPCR data*](https://doi.org/10.1093/bioinformatics/btu239).

---

## 4 — Recommendations

### Recommended tests

- **ANCOVA** is recommended when you have a clear reference/control sample (e.g. an untreated condition). It adjusts for inter-replicate variability across runs via the reference sample's dCq covariate.
- **Mixed Effect Model** is better suited when the reference sample is arbitrary across replicates (e.g. comparing expression between different patients or cell lines). It models replicate-level variance as a random intercept, making it more appropriate for designs without a fixed control.
- Testing on **ΔCq is preferred**. If a user instead tests $ΔΔCq$ or $2^{-ΔΔCq}$, repeated Welch's t-tests are the default because the fixed reference group has zero variance.
- Other tests (repeated paired t-tests, etc.) are available for completeness but are generally not the first choice.

### Variance assumptions

- On the **ΔCq** scale, data is often reasonably homoscedastic, although the assumption should still be checked. On the **ΔΔCq** scale, the reference group is fixed at zero and therefore has zero variance; repeated Welch's t-tests are used by default instead of an omnibus Welch's ANOVA.
- On the **2⁻ΔΔCq** (linear) scale, variance heterogeneity is much more common because qPCR measurement noise is multiplicative (proportional to expression level). The log scale (ΔCq/ΔΔCq) converts this multiplicative noise into additive noise by the property of logarithms, yielding homogeneous variance. If testing on 2⁻ΔΔCq, we recommend using a test that models unequal variances (Welch's t-test).

### Non-parametric tests

Non-parametric tests are provided for completeness but are generally not recommended for qPCR data. They require at least 5 biological replicates and tend to have lower statistical power compared to their parametric counterparts. They may still be useful when the sample size is sufficient.
