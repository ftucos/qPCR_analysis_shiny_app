library(stringr)

# try to automatically pick HK genes
is_HK <- function(x) {
    HK_regex <- "(ACT|B2M|EF1|GAPDH|GUSB|BGUS|HPRT|HMBS|PGK|PPIA|RPL|RPS|SDH|TBP|TUBULIN|TUB(1)?(A)?$|UBC|YWHAZ|18S)"
    
    # make case insensitive
    x <- str_to_upper(x) |>
        # remove ws, dashes etc
        str_remove_all("( |-|_)")
    
    str_detect(x, HK_regex)
}
