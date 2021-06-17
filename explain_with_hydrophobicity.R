# Usage:
#
# Rscript explain_with_hydrophobicity.R [target] [percentage] [MHC class]
#
# Rscript explain_with_hydrophobicity.R covid 2% MHC-I
#
# Results in a file called `[target]_[percentage]_counts.csv`,
# with a table like this:
#
#`haplotype`   |`binder`      |hydrophobicity
#--------------|--------------|--------------
# [h1]         |FAMILYVWW     |0.123
# [h1]         |VWFAMILYY     |-12.3
#
#
# Create all counts
args <- commandArgs(trailingOnly = TRUE)
if (1 == 2) {
  args <- c("covid", "2%", "MHC-I")
  args <- c("myco", "2%", "MHC-I")
  args <- c("human", "2%", "MHC-I")
  args <- c("human", "2%", "MHC-II")
}
message("args: {", paste0(args, collapse = ", "), "}")

testthat::expect_equal(length(args), 3)
target_name <- args[1]
message("target_name: ", target_name)
bbbq::check_target_name(target_name)

testthat::expect_true(nchar(args[2]) > 1)
testthat::expect_equal("%", stringr::str_sub(args[2], nchar(args[2]), nchar(args[2])))
percentage <- as.numeric(stringr::str_sub(args[2], 1, nchar(args[2]) - 1))
testthat::expect_true(percentage >= 0)
testthat::expect_true(percentage <= 100)
testthat::expect_true(percentage / 100 >= 0.0)
testthat::expect_true(percentage / 100 <= 1.0)
message("percentage: ", percentage)


testthat::expect_equal("MHC-", stringr::str_sub(args[3], 1, 4))
mhc_class_roman <- stringr::str_sub(args[3], 5)
message("mhc_class_roman: ", mhc_class_roman)

mhc_class <- NA
if (mhc_class_roman == "I") mhc_class <- 1
if (mhc_class_roman == "II") mhc_class <- 2
testthat::expect_false(is.na(mhc_class))

library(dplyr, warn.conflicts = FALSE, quietly = TRUE)

topology_prediction_tool <- "tmhmm"
message("topology_prediction_tool: ", topology_prediction_tool)
bbbq::check_topology_prediction_tool(topology_prediction_tool)

mhc1_ic50_prediction_tool <- "EpitopePrediction"
message("mhc1_ic50_prediction_tool: ", mhc1_ic50_prediction_tool)
bbbq::check_ic50_prediction_tool(mhc1_ic50_prediction_tool)

mhc2_ic50_prediction_tool <- "mhcnuggetsr"
message("mhc2_ic50_prediction_tool: ", mhc2_ic50_prediction_tool)
bbbq::check_ic50_prediction_tool(mhc2_ic50_prediction_tool)

haplotypes <- NA
if (mhc_class == 1) haplotypes <- bbbq::get_mhc1_haplotypes()
if (mhc_class == 2) haplotypes <- bbbq::get_mhc2_haplotypes()
testthat::expect_false(is.na(haplotypes[1]))
message("haplotypes: ", paste0(haplotypes, collapse = ", "))

if (target_name == "human") {
  proteome_type <- "representative"
} else {
  proteome_type <- "full"
}
message("proteome_type: ", proteome_type)

t_proteome <- bbbq::get_proteome(
  target_name = target_name,
  proteome_type = proteome_type,
  keep_selenoproteins = FALSE
)

t_topology <- bbbq::get_topology(
  target_name = target_name,
  proteome_type = proteome_type,
  keep_selenoproteins = FALSE,
  topology_prediction_tool = topology_prediction_tool
)
testthat::expect_true(all(t_proteome$name == t_topology$name))

tibbles <- list()

for (haplotype_index in seq_along(haplotypes)) {
  haplotype <- haplotypes[haplotype_index]
  message("haplotype: ", haplotype)

  peptide_length <- NA
  if (haplotype %in% bbbq::get_mhc1_haplotypes()) {
    message("Haplotype is MHC-I")
    peptide_length <- 9
    ic50_prediction_tool <- mhc1_ic50_prediction_tool
  } else {
    message("Haplotype is MHC-II")
    peptide_length <- 14
    ic50_prediction_tool <- mhc2_ic50_prediction_tool
  }
  message("peptide_length: ", peptide_length)
  message("ic50_prediction_tool: ", ic50_prediction_tool)

  # Remove the sequences and topologies that are too short
  testthat::expect_true(all(t_proteome$name == t_topology$name))
  testthat::expect_equal(nrow(t_proteome), nrow(t_topology))
  t_proteome <- t_proteome %>% filter(nchar(sequence) >= peptide_length)
  t_topology <- t_topology %>% filter(nchar(sequence) >= peptide_length)
  testthat::expect_equal(nrow(t_proteome), nrow(t_topology))
  testthat::expect_true(all(t_proteome$name == t_topology$name))

  # Keep only the sequences with the 20 standard amino acids
  regexp <- paste0("^[", paste0(Peptides::aaList(), collapse = ""), "]+$")
  keep_indices <- stringr::str_which(string = t_proteome$sequence, pattern = regexp)
  t_proteome <- t_proteome[keep_indices, ]
  t_topology <- t_topology[keep_indices, ]
  testthat::expect_equal(nrow(t_proteome), nrow(t_topology))
  testthat::expect_true(all(t_proteome$name == t_topology$name))

  t_proteome_n_mers <- bbbq::create_n_mers_tibble(
    strings = t_proteome$sequence,
    n = peptide_length
  )
  nrow(t_proteome_n_mers)
  t_proteome_n_mers$string <- NULL
  t_proteome_n_mers$name <- rep(t_proteome$name, times = nchar(t_proteome$sequence) - peptide_length + 1)

  t_topology_n_mers <- bbbq::create_n_mers_tibble(
    strings = t_topology$sequence,
    n = peptide_length
  )
  t_topology_n_mers$string <- NULL
  t_topology_n_mers$name <- rep(t_topology$name, times = nchar(t_topology$sequence) - peptide_length + 1)
  testthat::expect_equal(nrow(t_proteome_n_mers), nrow(t_topology_n_mers))
  testthat::expect_true(all(t_proteome_n_mers$name == t_topology_n_mers$name))
  t_topology_n_mers$overlap_with_tmh <- stringr::str_detect(
    string = t_topology_n_mers$n_mer,
    pattern = "[mM]"
  )

  t_ic50s <- bbbq::get_ic50s_lut(
    target_name = target_name,
    haplotype = haplotype,
    peptide_length = peptide_length,
    ic50_prediction_tool = ic50_prediction_tool
  )

  # Use all IC50 values in all proteomes' epitopes
  ic50_threshold <- as.numeric(quantile(x = t_ic50s$ic50, probs = percentage / 100))

  t_ic50s$is_binder <- t_ic50s$ic50 <= ic50_threshold
  testthat::expect_equal(0, sum(is.na(t_ic50s$is_binder)))

  testthat::expect_true("n_mer" %in% names(t_proteome_n_mers))
  testthat::expect_true("peptide" %in% names(t_ic50s))
  t <- dplyr::left_join(x = t_proteome_n_mers, y = t_ic50s, by = c("n_mer" = "peptide"))

  testthat::expect_equal(nrow(t), nrow(t_proteome_n_mers))
  testthat::expect_equal(nrow(t), nrow(t_topology_n_mers))
  testthat::expect_equal(t$name, t_topology_n_mers$name)
  t$overlap_with_tmh <- t_topology_n_mers$overlap_with_tmh


  f_tmh_epitopes <- mean(t$is_binder & t$overlap_with_tmh) / mean(t$overlap_with_tmh)
  hydrophobicity <- mean(Peptides::hydrophobicity(dplyr::filter(t, is_binder)$n_mer))

  # Show counts per protein
  t_here <- tibble::tibble(
    haplotype = haplotype,
    f_tmh_epitopes = f_tmh_epitopes,
    hydrophobicity = hydrophobicity
  )
  message(knitr::kable(head(t_here)))
  tibbles[[haplotype_index]] <- t_here

}


t <- dplyr::bind_rows(tibbles)

t$simple_haplotype <- bbbq::simplify_haplotype_names(t$haplotype)

p <- ggplot2::ggplot(t, ggplot2::aes(x = f_tmh_epitopes, y = hydrophobicity)) +
  ggplot2::geom_point() +
  ggplot2::geom_text(
    ggplot2::aes(label = simple_haplotype),
    size = 5,
    nudge_y = 0.04,
    check_overlap = TRUE
  ) +
  ggplot2::scale_y_continuous("Hydrophobicity preference score") +
  ggplot2::scale_x_continuous("% of TMH epitopes") +
  ggpmisc::stat_poly_eq(
    formula = y ~ x,
    ggplot2::aes(label = paste(..rr.label.., sep = "~~~")),
    parse = TRUE,
    size = 8
  ) + bbbq::get_bbbq_theme()
p
p + ggplot2::geom_smooth(method = "lm", formula = y ~ x, col = "red", se = FALSE); ggplot2::ggsave(
    paste0("~/fig_hydrophobicity_mhc", mhc_class,".png"),
    width = 7,
    height = 7
  )
p + ggplot2::geom_smooth(method = "lm", formula = y ~ x, col = "red", se = FALSE); ggplot2::ggsave(
    paste0("~/fig_hydrophobicity_mhc", mhc_class,".tiff"),
    width = 7,
    height = 7
  )

p + ggplot2::geom_smooth(method = "lm", formula = y ~ x, col = "black", se = FALSE); ggplot2::ggsave(
    paste0("~/fig_hydrophobicity_mhc", mhc_class,"_bw.png"),
    width = 7,
    height = 7
  )
p + ggplot2::geom_smooth(method = "lm", formula = y ~ x, col = "black", se = FALSE); ggplot2::ggsave(
    paste0("~/fig_hydrophobicity_mhc", mhc_class,"_bw.tiff"),
    width = 7,
    height = 7
  )
