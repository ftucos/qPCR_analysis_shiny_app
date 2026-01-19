data.frame(value = c(1, 2, NA, 1, 2, 3),
           sample = c("A", "A", "A", "B", "B", "B")) |> 
    group_by(sample) |>
    summarize(count = length(na.omit(value)))

# compare qPCR SD and SE propagation

n_samples = 2
n_genes = 3
n_rep = 4
raw_data <- data.frame(
    Sample = rep(paste0("Sample", seq(1:n_samples)), each = n_rep * n_genes),
    Target = rep(rep(c("HK1", "HK2", "TG"), each = n_rep), n_samples),
    Cq     = round(rnorm(n_samples * n_genes * n_rep, mean = 23.5, sd = 1.2), 2)
) 

raw_data <- raw_data[2:nrow(raw_data), ] # drop first row to unbalance n

hk_summary_int <- raw_data |>
    filter(Target %in% c("HK1", "HK2")) |>
    group_by(Sample, Target) |>
    summarize(
        HK_mean = mean(Cq, na.rm = TRUE),
        HK_sd   = sd(Cq, na.rm = TRUE),
        HK_n    = length(na.omit(Cq)), # sample size for each gene
        HK_se   = HK_sd/sqrt(HK_n),
        .groups = "drop"
    )

hk_summary <- hk_summary_int |>
    group_by(Sample) |>
    summarize(
        HK_mean = mean(HK_mean, na.rm = TRUE),
        n_HK_genes = sum(HK_n > 0),
        # pooled SD and SE formula: http://stats.libretexts.org/Bookshelves/Applied_Statistics/An_Introduction_to_Psychological_Statistics_(Foster_et_al.)/10%3A__Independent_Samples/10.05%3A_Standard_Error_and_Pooled_Variance
        HK_sd_pool   = sqrt(
            sum(HK_sd^2*(HK_n-1))/
                (sum(HK_n)-n_HK_genes)
        ), 
        HK_avg_sd = (1/2)*sqrt(sum(HK_sd^2)),
        HK_se   = (1/n_HK_genes)*sqrt(sum(HK_sd_pool^2/HK_n)),
        HK_se2  = (1/n_HK_genes)*sqrt(sum(HK_sd^2/HK_n)),
        .groups = "drop"
    )
