# Create a table that shows the IC50 values below
# which a peptide is considered a binder
#
# Usage:
#
#  Rscript create_table_ic50_binders.R
#
library(dplyr, warn.conflicts = FALSE)
library(knitr, warn.conflicts = FALSE)
library(testthat, warn.conflicts = FALSE)

csv_target_filename <- "table_ic50_binders.csv"
message("'csv_target_filename': '", csv_target_filename, "'")
latex_target_filename <- "table_ic50_binders.latex"
message("'latex_target_filename': '", latex_target_filename, "'")

percentile <- bbbq::get_ic50_percentile_binder()
message("'percentile': '", percentile, "' (as hard-coded by BBBQ)")

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


t <- t_haplotypes
t$haplotype_id <- NULL
t$ic50 <- NA


for (i in seq_len(nrow(t))) {
  t$ic50[i] <- mhcnpreds::get_ic50_threshold(
    peptide_length = bbbq::get_mhc_peptide_length(t$mhc_class[i]),
    mhc_haplotype = mhcnuggetsr::to_mhcnuggets_name(t$haplotype[i]),
    percentile = percentile
  )
}
t$mhc_class <- NULL

readr::write_csv(t, csv_target_filename)

knitr::kable(
  t, "latex",
  caption = paste0(
    "IC50 values (in nM) per haplotype ",
    "below which a peptide is considered a binder. ",
    "Percentile used: ", percentile
  ),
  label = "tab:ic50_binders"
) %>% cat(., file = latex_target_filename)
