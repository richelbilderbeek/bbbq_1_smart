# Create 'table_f_tmh.latex'
#
# Usage:
#
#   Rscript create_table_f_tmh.R
#
# For example:
#
#   Rscript create_table_f_tmh.R
#
# Output:
#
#  * File named 'table_f_tmh.latex'
#
library(testthat, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)

target_latex_filename <- "table_f_tmh.latex"
message("'target_latex_filename': '", target_latex_filename, "'")

haplotypes_filename <- "haplotypes_lut.csv"
message("'haplotypes_filename': '", haplotypes_filename, "'")
testthat::expect_true(file.exists(haplotypes_filename))
t_haplotypes <- readr::read_csv(
  haplotypes_filename,
  col_types = readr::cols(
    haplotype = readr::col_character(),
    mhc_class = readr::col_double(),
    haplotype_id = readr::col_character()
  )
)

t <- readr::read_csv("counts.csv",
  col_types = readr::cols(
    target = readr::col_character(),
    haplotype_id = readr::col_character(),
    protein_id = readr::col_character(),
    n_binders = readr::col_skip(),
    n_binders_tmh = readr::col_skip(),
    n_spots = readr::col_double(),
    n_spots_tmh = readr::col_double()
  )
)

# Group by proteins
t <- t %>% dplyr::group_by(target, haplotype_id) %>%
  dplyr::summarise(
    n_spots = sum(n_spots),
    n_spots_tmh = sum(n_spots_tmh),
    .groups = "keep"
  )

# Add mhc_class
t$mhc_class <- NA
for (i in seq_len(nrow(t))) {
  t$mhc_class[i] <- t_haplotypes$mhc_class[t_haplotypes$haplotype_id == t$haplotype_id[i] ]
}

# Group by proteins
t <- t %>% dplyr::group_by(target, mhc_class) %>%
  dplyr::summarise(
    n_spots = mean(n_spots),
    n_spots_tmh = mean(n_spots_tmh),
    .groups = "keep"
  )

t$f_tmh <- format(100.0 * (t$n_spots_tmh / t$n_spots), digits = 3)

knitr::kable(
  t, "latex",
  caption = paste0(
    "Percentage of spots and spots that overlap with a TMH"
  ),
  label = "f_tmh"
) %>% cat(., file = target_latex_filename)
