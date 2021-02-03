# Merges all tables for MHC-I or MHC-II,
# that have measured the number of binders and the
# number of binders that are TMH
#
# Usage:
#
#  Rscript create_table_tmh_binders_mhc.R [MHC]
#
#  * [MHC] is either 'mhc1' or 'mhc2'
#

library(dplyr, warn.conflicts = FALSE)
library(testthat, warn.conflicts = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args <- "mhc1"
}
expect_equal(length(args), 1)
message("Running with argument '", args[1], "'")
mhc <- args[1]
message("mhc: '", mhc, "'")
expect_equal(4, stringr::str_length(mhc))
mhc_class <- stringr::str_sub(mhc, 4, 4)
message("mhc_class: '", mhc_class, "'")
the_mhc_class <- mhc_class # Needed for filtering later

target_csv_filename <- paste0("table_tmh_binders_mhc", mhc_class, ".csv")
message("target_csv_filename: '", target_csv_filename, "'")

target_latex_filename <- paste0("table_tmh_binders_mhc", mhc_class, ".latex")
message("target_latex_filename: '", target_latex_filename, "'")

percentile <- bbbq::get_ic50_percentile_binder()
message("'percentile': '", percentile, "' (as hard-coded by BBBQ)")


raw_table_filename <- "counts.csv"
testthat::expect_true(file.exists(raw_table_filename))
t_raw <- readr::read_csv(
  raw_table_filename,
  col_types = readr::cols(
    target = readr::col_character(),
    haplotype_id = readr::col_character(),
    protein_id = readr::col_character(),
    n_binders = readr::col_double(),
    n_binders_tmh = readr::col_double(),
    n_spots = readr::col_double(),
    n_spots_tmh = readr::col_double()
  )
)
# Create the BBBQ haplotype LUT
haplotypes_filename <- "haplotypes_lut.csv"
message("haplotypes_filename: '", haplotypes_filename, "'")
expect_true(file.exists(haplotypes_filename))
t_haplotypes <- readr::read_csv(
  haplotypes_filename,
  col_types = readr::cols(
    haplotype = readr::col_character(),
    mhc_class = readr::col_double(),
    haplotype_id = readr::col_character()
  )
)
t_haplotypes$name <- mhcnuggetsr::to_mhcnuggets_names(t_haplotypes$haplotype)

# Only keep the desired MHC class
t_haplotypes <- t_haplotypes %>% filter(mhc_class == the_mhc_class)

t_long <- t_raw %>% dplyr::filter(haplotype_id %in% t_haplotypes$haplotype_id)

t_long$target <- as.factor(t_long$target)
t_long$haplotype_id <- as.factor(t_long$haplotype_id)
t_long$protein_id <- as.factor(t_long$protein_id)

# Group all proteins
t_long <- t_long %>% dplyr::group_by(target, haplotype_id) %>%
    dplyr::summarize(
      n_binders = sum(n_binders, na.rm = TRUE),
      n_binders_tmh = sum(n_binders_tmh, na.rm = TRUE),
      n_spots = sum(n_spots, na.rm = TRUE),
      n_spots_tmh = sum(n_spots_tmh, na.rm = TRUE),
      .groups = "keep"
    ) %>% dplyr::ungroup()

t_long$f <- 100.0 * t_long$n_binders_tmh / t_long$n_binders
t_long$f <- paste0(
  format(t_long$f, digits = 4), " ",
  "(", t_long$n_binders_tmh, "/", t_long$n_binders, ")"
)
t_long$haplotype <- NA
for (i in seq_len(nrow(t_long))) {
  id <- t_long$haplotype_id[i]
  t_long$haplotype[i] <- t_haplotypes$haplotype[id == t_haplotypes$haplotype_id]
}
t_long <- t_long %>% dplyr::select(target, haplotype, f)
t_long$target <- as.factor(t_long$target)

# Wide form
t_wide <- tidyr::pivot_wider(
  t_long,
  names_from = "target",
  values_from = "f"
)
readr::write_csv(t_wide, target_csv_filename)


roman_mhc_class <- NA
if (mhc_class == 1) roman_mhc_class <- "I"
if (mhc_class == 2) roman_mhc_class <- "II"

knitr::kable(
  t_wide, "latex",
  caption = paste0(
    "Percentage of MHC-", roman_mhc_class, " ",
    bbbq::get_mhc_peptide_length(mhc_class) ,"-mers overlapping with TMH. ",
    "Values in brackets show the number of binders ",
    "that have at least one residue overlapping with a TMH (first value)",
    "as well as the number of binders (second value). ",
    "Percentile used: ", percentile
  ),
  label = paste0("tmh_binders_mhc", mhc_class)
) %>% cat(., file = target_latex_filename)
