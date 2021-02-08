# Create 'table_f_tmh.latex'
#
# Usage:
#
#   Rscript create_table_f_tmh.R [percentage]
#
# For example:
#
#   Rscript create_table_f_tmh.R 2
#
# Output:
#
#  * File named 'table_f_tmh_[percentage].latex',
#    for example 'table_f_tmh_2.latex
#
args <- commandArgs(trailingOnly = TRUE)
if (1 == 2) {
  args <- c("2")
}
testthat::expect_equal(length(args), 1)
message("Running with arguments {", paste0(args, collapse = ", "), "}")
percentage <- as.numeric(args[1])
message("percentage: ", percentage)

suppressPackageStartupMessages({
  library(testthat, warn.conflicts = FALSE)
  library(dplyr, warn.conflicts = FALSE)
})

counts_csv_filename <- paste0("counts_", percentage, ".csv")
message("counts_csv_filename: ", counts_csv_filename)
testthat::expect_true(file.exists(counts_csv_filename))

target_latex_filename <- paste0("table_f_tmh_", percentage, ".latex")
message("target_latex_filename: ", target_latex_filename, " (the output file)")

t <- readr::read_csv(counts_csv_filename,
  col_types = readr::cols(
    target = readr::col_character(),
    haplotype = readr::col_character(),
    name = readr::col_character(),
    n_binders = readr::col_skip(),
    n_binders_tmh = readr::col_skip(),
    n_spots = readr::col_double(),
    n_spots_tmh = readr::col_double()
  )
)

# Group by proteins
t <- t %>% dplyr::group_by(target, haplotype) %>%
  dplyr::summarise(
    n_spots = sum(n_spots),
    n_spots_tmh = sum(n_spots_tmh),
    .groups = "keep"
  )

# Add mhc_class
t$mhc_class <- NA
for (i in seq_len(nrow(t))) {
  if (t$haplotype[i] %in% bbbq::get_mhc1_haplotypes()) {
    t$mhc_class[i] <- 1
  } else {
    testthat::expect_true(t$haplotype[i] %in% bbbq::get_mhc2_haplotypes())
    t$mhc_class[i] <- 2
  }
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

testthat::expect_true(file.exists(target_latex_filename))
