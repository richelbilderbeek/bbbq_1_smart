# Create a table that shows the IC50 values below
# which a peptide is considered a binder
#
# Usage:
#
#  Rscript create_table_ic50_binders.R [mhc class] [percentage]
#
# Example:
#
#  Rscript create_table_ic50_binders.R mhc1 2%
#
args <- commandArgs(trailingOnly = TRUE)
if (1 == 2) {
  args <- c("mhc1", "2%")
}
testthat::expect_equal(length(args), 2)
message("Running with arguments {", paste0(args, collapse = ", "), "}")

testthat::expect_equal(4, nchar(args[1]))
mhc_class <- as.numeric(stringr::str_sub(args[1], 4, 4))
message("mhc_class: ", mhc_class)

testthat::expect_true(nchar(args[2]) > 1)
testthat::expect_equal("%", stringr::str_sub(args[2], nchar(args[2]), nchar(args[2])))
percentage <- as.numeric(stringr::str_sub(args[2], 1, nchar(args[2]) - 1))
message("percentage: ", percentage)

suppressPackageStartupMessages({
  library(dplyr, warn.conflicts = FALSE)
  library(knitr, warn.conflicts = FALSE)
  library(testthat, warn.conflicts = FALSE)
})

csv_target_filename <- paste0("table_ic50_binders_mhc", mhc_class, "_", percentage, ".csv")
message("'csv_target_filename': '", csv_target_filename, "'")
latex_target_filename <- paste0("table_ic50_binders_mhc", mhc_class, "_", percentage, ".latex")
message("'latex_target_filename': '", latex_target_filename, "'")

if (mhc_class == 1) {
  haplotypes <- bbbq::get_mhc1_haplotypes()
  peptide_length <- 9
  ic50_prediction_tool <- "EpitopePrediction"
} else {
  testthat::expect_equal(2, mhc_class)
  haplotypes <- bbbq::get_mhc2_haplotypes()
  peptide_length <- 15
  ic50_prediction_tool <- "mhcnuggetsr"
}
message("haplotypes: {", paste0(haplotypes, collapse = ", "), "}")
message("peptide_length: ", peptide_length)
message("ic50_prediction_tool: ", ic50_prediction_tool)

t <- tidyr::expand_grid(
  target = c("covid", "human", "myco"),
  haplotype = haplotypes,
  ic50 = NA
)

for (i in seq_len(nrow(t))) {
  target_name <- t$target[i]
  haplotype <- t$haplotype[i]
  t_ic50s <- bbbq::get_ic50s_lut(
    target_name = target_name,
    haplotype = haplotype,
    peptide_length = peptide_length,
    ic50_prediction_tool = ic50_prediction_tool
  )
  t$ic50[i] <- as.numeric(quantile(x = t_ic50s$ic50, probs = percentage / 100))
}

readr::write_csv(t, csv_target_filename)

t_wide <- tidyr::pivot_wider(
  t,
  names_from = "target",
  values_from = "ic50"
)

knitr::kable(
  t_wide, "latex",
  caption = paste0(
    "IC50 values (in nM) per haplotype ",
    "below which a peptide is considered a binder. ",
    "percentage used: ", percentage
  ),
  label = paste0("ic50_binders_mhc", mhc_class)
) %>% cat(., file = latex_target_filename)
