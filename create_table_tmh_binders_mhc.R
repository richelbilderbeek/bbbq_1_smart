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
args <- commandArgs(trailingOnly = TRUE)
if (1 == 2) {
  args <- c("mhc1", "2")
}
message("Running with arguments {", paste0(args, collapse = ", "), "}")
testthat::expect_equal(length(args), 2)

mhc <- args[1]
message("mhc: '", mhc, "'")
testthat::expect_equal(4, stringr::str_length(mhc))

percentage <- as.numeric(args[2])
message("percentage: ", percentage)


mhc_class <- as.numeric(stringr::str_sub(mhc, 4, 4))
message("mhc_class: '", mhc_class, "'")
the_mhc_class <- mhc_class # Needed for filtering later

haplotypes <- NA
if (mhc_class == 1) {
  haplotypes <- bbbq::get_mhc1_haplotypes()
  peptide_length <- 9
} else {
  testthat::expect_equal(2, mhc_class)
  haplotypes <- bbbq::get_mhc2_haplotypes()
  peptide_length <- 15
}
message("haplotypes: {", paste0(haplotypes, collapse = ", "), "}")
message("peptide_length: ", peptide_length)

target_csv_filename <- paste0("table_tmh_binders_mhc", mhc_class, "_", percentage, ".csv")
message("target_csv_filename: '", target_csv_filename, "'")

target_latex_filename <- paste0("table_tmh_binders_mhc", mhc_class, "_", percentage, ".latex")
message("target_latex_filename: '", target_latex_filename, "'")

suppressPackageStartupMessages({
  library(dplyr, warn.conflicts = FALSE)
  library(testthat, warn.conflicts = FALSE)
})


raw_table_filename <- paste0("counts_", percentage, ".csv")
testthat::expect_true(file.exists(raw_table_filename))
t_haplotypes <- readr::read_csv(
  raw_table_filename,
  col_types = readr::cols(
    target = readr::col_character(),
    haplotype = readr::col_character(),
    name = readr::col_character(),
    n_binders = readr::col_double(),
    n_binders_tmh = readr::col_double(),
    n_spots = readr::col_double(),
    n_spots_tmh = readr::col_double()
  )
)
nrow(t_haplotypes)

# Create the BBBQ haplotype LUT

# Only keep the desired MHC class

t_long <- t_haplotypes %>% filter(haplotype %in% haplotypes)
nrow(t_long)


#t_long$target <- as.factor(t_long$target)
#t_long$haplotype_id <- as.factor(t_long$haplotype_id)
#t_long$protein_id <- as.factor(t_long$protein_id)

# Group all proteins
t_long <- t_long %>% dplyr::group_by(target, haplotype) %>%
    dplyr::summarize(
      n_binders = sum(n_binders),
      n_binders_tmh = sum(n_binders_tmh),
      n_spots = sum(n_spots),
      n_spots_tmh = sum(n_spots_tmh),
      .groups = "keep"
    ) %>% dplyr::ungroup()

t_long$f <- 100.0 * t_long$n_binders_tmh / t_long$n_binders
t_long$f <- paste0(
  format(t_long$f, digits = 4), " ",
  "(", t_long$n_binders_tmh, "/", t_long$n_binders, ")"
)
t_long <- t_long %>% dplyr::select(target, haplotype, f)
t_long$target <- as.factor(t_long$target)

# Wide form
t_wide <- tidyr::pivot_wider(
  t_long,
  names_from = "target",
  values_from = "f"
)
knitr::kable(head(t_wide))

readr::write_csv(t_wide, target_csv_filename)


knitr::kable(
  t_wide, "latex",
  caption = paste0(
    "Percentage of MHC-", as.roman(mhc_class), " ",
    peptide_length ,"-mers overlapping with TMH. ",
    "Values in brackets show the number of binders ",
    "that have at least one residue overlapping with a TMH (first value)",
    "as well as the number of binders (second value). ",
    "percentage used: ", percentage
  ),
  label = paste0("tmh_binders_mhc", mhc_class)
) %>% cat(., file = target_latex_filename)
