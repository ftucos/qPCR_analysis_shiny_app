library(readr)
library(dplyr)
library(stringr)

parse_Cq <- function(x) {
    x <- str_trim(as.character(x))
    x <- case_when(
        # undetermined or other non-numeric labels
        str_detect(x, "[A-Za-z]{2,}") ~ "Inf",
        # higher than
        str_detect(x, ">") ~ "Inf",
        # , as decimal mark
        str_count(x, ",") == 1 & str_count(x, "\\.") == 0 ~ suppressWarnings(parse_number(x, locale = locale(decimal_mark = ","))) |> as.character(),
        TRUE ~ suppressWarnings(parse_number(x)) |> as.character()
    )
}