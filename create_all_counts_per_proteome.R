# Usage:
#
# Rscript create_all_counts_per_proteome.R [target] [percentage]
#
# Rscript create_all_counts_per_proteome.R human 5
#
# Results in a file called `[target]_[percentage]_counts.csv`,
# with a table like this:
#
#`haplotype`   |`protein_name`|`n_binders`|`n_binders_tmh`|`n_spots`|`n_spots_tmh`
#--------------|--------------|-----------|---------------|---------|-------------
# [h1]         |[p1]          |11         |5              |100      |20
# [h1]         |[p2]          |12         |4              |10       |2
#
# [h1] = HLA-A*01:01
# [p1] = tr|A0A024R1R8|A0A024R1R8_HUMAN Coiled-coil domain-containing protein 72 OS=Homo sapiens (Human) OX=9606 GN=ENSG00000225528 PE=3 SV=1 # nolint indeed a long line
#
# Create all counts
args <- commandArgs(trailingOnly = TRUE)
message("args: {", paste0(args, collapse = ", "), "}")

if (1 == 2) {
  args <- c("human", "5")
}
testthat::expect_equal(length(args), 2)
target_name <- args[1]
message("target_name: ", target_name)
bbbq::check_target_name(target_name)
percentage <- as.numeric(args[2]) # use percentage values for filename
message("percentage: ", percentage)
testthat::expect_true(percentage >= 0)
testthat::expect_true(percentage <= 100)
testthat::expect_true(percentage / 100 >= 0.0)
testthat::expect_true(percentage / 100 <= 1.0)

target_filename <- paste0(target_name, "_", percentage, "_counts.csv")
message("target_filename: ", target_filename)

topology_prediction_tool <- "tmhmm"
message("topology_prediction_tool: ", topology_prediction_tool)
bbbq::check_topology_prediction_tool(topology_prediction_tool)
ic50_prediction_tool <- "EpitopePrediction"
message("ic50_prediction_tool: ", ic50_prediction_tool)
bbbq::check_ic50_prediction_tool(ic50_prediction_tool)

haplotypes <- NA
if (ic50_prediction_tool == "EpitopePrediction") {
  message("Use MHC-I haplotypes")
  haplotypes <- bbbq::get_mhc1_haplotypes()
} else {
  stop("Not implemented yet")
}

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
  } else {
    message("Haplotype is MHC-II")
    peptide_length <- 15
  }
  message("peptide_length: ", peptide_length)

  # Remove the sequences and topologies that are too short
  testthat::expect_true(all(t_proteome$name == t_topology$name))
  testthat::expect_equal(nrow(t_proteome), nrow(t_topology))
  library(dplyr)
  t_proteome <- t_proteome %>% filter(nchar(sequence) >= peptide_length)
  t_topology <- t_topology %>% filter(nchar(sequence) >= peptide_length)
  testthat::expect_equal(nrow(t_proteome), nrow(t_topology))
  testthat::expect_true(all(t_proteome$name == t_topology$name))

  # Keep only the sequences with a TMHs
  keep_indices <- stringr::str_which(string = t_topology$sequence, pattern = "[mM]")
  t_proteome <- t_proteome[keep_indices, ]
  t_topology <- t_topology[keep_indices, ]
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
    haplotype = haplotype,
    peptide_length = peptide_length,
    ic50_prediction_tool = ic50_prediction_tool
  )

  # Use all IC50 values in all proteomes' epitopes
  ic50_threshold <- as.numeric(quantile(x = t_ic50s$ic50, probs = percentage / 100))

  t_ic50s$is_binder <- t_ic50s$ic50 <= ic50_threshold
  sum(is.na(t_ic50s$is_binder))
  which(is.na(t_ic50s$is_binder))
  which(is.na(as.numeric(t_ic50s$is_binder)))

  t <- dplyr::left_join(x = t_proteome_n_mers, y = t_ic50s, by = c("n_mer" = "peptide"))

  testthat::expect_equal(nrow(t), nrow(t_proteome_n_mers))
  testthat::expect_equal(nrow(t), nrow(t_topology_n_mers))
  testthat::expect_equal(t$name, t_topology_n_mers$name)
  t$overlap_with_tmh <- t_topology_n_mers$overlap_with_tmh

  t_counts <- t %>% group_by(name) %>% summarise(
    n_binders = sum(is_binder),
    n_binders_tmh = sum(is_binder & overlap_with_tmh),
    n_spots = n(),
    n_spots_tmh = sum(overlap_with_tmh),
    .groups = "drop"
  )
  # Show counts per protein
  message(knitr::kable(head(t_counts)))
  t_counts$haplotype <- haplotype
  tibbles[[haplotype_index]] <- t_counts

  # Other debug info
  f_n_mers_overlap_with_tmh <- sum(t$overlap_with_tmh) / nrow(t)
  message("f_n_mers_overlap_with_tmh: ", f_n_mers_overlap_with_tmh)
  sum(is.na(t$is_binder))
  which(is.na(t$is_binder))
  which(is.na(as.numeric(t$is_binder)))
  f_binder <- sum(as.numeric(t$is_binder)) / nrow(t)
  message("f_binder: ", f_binder)
  f_binders_that_overlap_with_tmh <- sum(t$overlap_with_tmh & t$is_binder) / sum(t$is_binder)
  message("f_binders_that_overlap_with_tmh: ", f_binders_that_overlap_with_tmh)
}

readr::write_csv(
  dplyr::bind_rows(tibbles),
  target_filename
)
