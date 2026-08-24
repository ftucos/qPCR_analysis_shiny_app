# Example qPCR data scenarios

Load `example_qPCR_data.csv` with the app's **Load Example** button. The data
contain four conditions, five biological replicates (`R1`–`R5`), and three
technical measurements per biological replicate. Select both `ACTB` and `TBP`
as housekeeping genes; the app should select them automatically.

| Target or HK | Scenario | Expected live behavior |
|---|---|---|
| `Target_Stable` | No modulation in any condition | Similar normalized expression in all four conditions |
| `Target_Induced_A` | Only `Treatment_A` is induced | Approximately 8-fold induction relative to `Control`; the other conditions remain near control |
| `Target_Repressed_B` | Only `Treatment_B` is repressed | Approximately 8-fold repression relative to `Control` |
| `Target_Graded` | Ordered modulation | Expression increases from `Control` through `Treatment_C` |
| `Target_TechCensor` | One of three technical measurements is undetected in `Treatment_A / R2` | The undetected well is shown, but the replicate mean uses the two detected wells and is not censored |
| `Target_BioCensor` | All technical measurements are undetected in `Treatment_B / R3` | That biological replicate is censored; the other `Treatment_B` replicates remain detected |
| `Target_ReferenceCensor` | All technical measurements are undetected in `Control / R4` | `Control` is marked as containing undetected data and cannot be used as the ΔΔCq/ANCOVA reference for this target; another condition can be selected |
| `Target_SampleAbsent` | `Treatment_C` is undetected in every biological replicate | The complete sample-target combination is censored |
| `Target_TwoSamplesAbsent` | `Treatment_B` and `Treatment_C` are always undetected | The app warns that the two fully censored samples cannot be ranked or meaningfully compared |
| `Target_AllAbsent` | Every measurement is undetected | All samples are censored and no valid reference exists for this target |
| `ACTB` | One technical measurement is undetected in `Treatment_A / R2` | The biological sample remains HK-valid because detected ACTB measurements exist |
| `TBP` | All technical measurements are undetected in `Treatment_B / R4` | Only `Treatment_B / R4` fails HK validation and is excluded from ΔCq processing; the rest of the run remains usable |
| `Target_HighCq` | One detected value is `40.40` | The maximum-cycle replacement automatically increases from 40 to 41 |

The censoring examples deliberately use `Undetermined`, `No Ct`, and `>40` to
exercise the accepted non-detect input formats. They should all be normalized
to the same explicit censored display using the current replacement cycle.

Regenerate the CSV from the repository root with:

```r
source("data/generate_example_qPCR_data.R")
```

