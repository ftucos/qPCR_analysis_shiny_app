library(readr)
library(dplyr)
library(stringr)

parse_Cq <- function(x) {
    original <- str_trim(as.character(x))

    # Treat explicit non-detect labels and censored values (for example,
    # ">40") as censored. Missing and otherwise invalid values are not censored.
    Cq_censored <- str_detect(original, "[A-Za-z]{2,}") |
        str_detect(original, fixed(">"))
    Cq_censored <- replace_na(Cq_censored, FALSE)

    uses_decimal_comma <- str_count(original, ",") == 1 &
        str_count(original, "\\.") == 0

    Cq <- case_when(
        Cq_censored ~ NA,
        uses_decimal_comma ~ suppressWarnings(
            parse_number(original, locale = locale(decimal_mark = ","))
        ),
        TRUE ~ suppressWarnings(parse_number(original))
    )

    tibble(Cq = Cq, Cq_censored = Cq_censored)
}

parse_Cq_data <- function(data) {
    if (!"Cq" %in% names(data)) {
        stop("`data` must contain a `Cq` column.")
    }

    parsed_cq <- parse_Cq(data$Cq)

    data |>
        mutate(
            Cq = parsed_cq$Cq,
            Cq_censored = parsed_cq$Cq_censored
        )
}
