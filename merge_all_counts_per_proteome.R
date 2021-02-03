# Merge all counts

# Merge files named `[target]_[percentage]_counts.csv` to `counts.csv`.
#
# `[target]_[percentage]_counts.csv`:
#
# name | n_binders | n_binders_tmh | n_spots | n_spots_tmh| haplotype
# [1]  |       221 |            51 |    4397 |        416 |	HLA-A*01:01
#
# * [1] Protein name, e.g. sp|P0DTC1|R1A_SARS2 Replicase polyprotein 1a OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 PE=1 SV=1
#
# `counts_[percentage].csv` must look like this:
#
# target | haplotype 	 | name | n_binders | n_binders_tmh | n_spots | n_spots_tmh
# covid  | HLA-A*01:01 | [1]  |        11 |             5 |     100 |          20
# covid  |         ... |  ... |       ... |           ... |     ... |         ...
# covid  | HLA-A*02:01 |  ... |       ... |           ... |     ... |         ...
# covid  |         ... |  ... |       ... |           ... |     ... |         ...
# human  |         ... |  ... |       ... |           ... |     ... |         ...
#   ...  |         ... |  ... |       ... |           ... |     ... |         ...
# myco   |         ... |  ... |       ... |           ... |     ... |         ...
#   ...  |         ... |  ... |       ... |           ... |     ... |         ...
#

args <- commandArgs(trailingOnly = TRUE)
if (1 == 2) {
  args <- c("5")
}
message("Running with arguments {", paste0(args, collapse = ", "), "}")
testthat::expect_equal(length(args), 1)
percentage <- as.numeric(args[1])
message("percentage: ", percentage)
testthat::expect_true(percentage >= 0)
testthat::expect_true(percentage <= 100)

target_filename <- paste0("counts_", percentage, ".csv")
message("target_filename: ", target_filename)

targets <- bbbq::get_target_names()
targets <- c("covid", "human", "myco")
message("Running for targets {", paste0(targets, collapse = ", "), "}")

# All little tibbles
tibbles <- list()

for (i in seq_along(targets)) {
  target <- targets[i]
  target_counts_filename <- paste0(
    target, "_", percentage, "_counts.csv"
  )
  testthat::expect_true(file.exists(target_counts_filename))

  t <- readr::read_csv(
    target_counts_filename,
    col_types = readr::cols(
      haplotype = readr::col_character(),
      name = readr::col_character(),
      n_binders = readr::col_double(),
      n_binders_tmh = readr::col_double(),
      n_spots = readr::col_double(),
      n_spots_tmh = readr::col_double()
    )
  )
  t <- dplyr::relocate(t, haplotype)
  t$target <- target
  t <- dplyr::relocate(t, target)
  tibbles[[i]] <- t
}

t <- dplyr::bind_rows(tibbles)

testthat::expect_equal(
  names(t),
  c("target", "haplotype", "name", "n_binders", "n_binders_tmh", "n_spots", "n_spots_tmh")
)

readr::write_csv(t, target_filename)

