# Create a table that shows the IC50 values below
# which a peptide is considered a binder
#
# Usage:
#
#  Rscript create_table_ic50_binders.R [percentage]
#
# Example:
#
#  Rscript create_table_ic50_binders.R 2
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
  library(dplyr, warn.conflicts = FALSE)
  library(knitr, warn.conflicts = FALSE)
  library(testthat, warn.conflicts = FALSE)
})

csv_target_filename <- paste0("table_ic50_binders_", percentage, ".csv")
message("'csv_target_filename': '", csv_target_filename, "'")
latex_target_filename <- paste0("table_ic50_binders_", percentage, ".latex")
message("'latex_target_filename': '", latex_target_filename, "'")

t <- tidyr::expand_grid(
  target = c("covid", "human", "myco"),
  haplotype = bbbq::get_mhc_haplotypes(),
  mhc_class = NA,
  peptide_length = NA,
  ic50 = NA,
  ic50_prediction_tool = ""
)

for (i in seq_len(nrow(t))) {
  if (t$haplotype[i] %in% bbbq::get_mhc1_haplotypes()) {
    t$mhc_class[i] <- 1
    t$peptide_length[i] <- 9
    t$ic50_prediction_tool[i] <- "EpitopePrediction"
  } else {
    testthat::expect_true(t$haplotype[i] %in% bbbq::get_mhc2_haplotypes())
    t$mhc_class[i] <- 2
    t$peptide_length[i] <- 15
    t$ic50_prediction_tool[i] <- "mhcnuggetsr"
  }
}
t
for (i in seq_len(nrow(t))) {
  target_name <- t$target[i]
  haplotype <- t$haplotype[i]
  peptide_length <- t$peptide_length[i]
  ic50_prediction_tool <- t$ic50_prediction_tool[i]
  t_ic50s <- bbbq::get_ic50s_lut(
    target_name = target_name,
    haplotype = haplotype,
    peptide_length = peptide_length,
    ic50_prediction_tool = ic50_prediction_tool
  )
  t$ic50[i] <- as.numeric(quantile(x = t_ic50s$ic50, probs = percentage / 100))
}

readr::write_csv(t, csv_target_filename)

knitr::kable(
  t, "latex",
  caption = paste0(
    "IC50 values (in nM) per haplotype ",
    "below which a peptide is considered a binder. ",
    "percentage used: ", percentage
  ),
  label = "tab:ic50_binders"
) %>% cat(., file = latex_target_filename)
