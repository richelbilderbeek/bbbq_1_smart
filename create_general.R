# Create the 'general.csv' file, containing general info

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
message("args: {", paste0(args, collapse = ", "), "}")

if (1 == 2) {
  args <- c("1AA")
  args <- c("2AA")
}

testthat::expect_equal(length(args), 1)
testthat::expect_true(nchar(args[1]) > 1)
testthat::expect_equal("AA", stringr::str_sub(args[1], 2, nchar(args[1])))
n_aas_overlap <- as.numeric(stringr::str_sub(args[1], 1, nchar(args[1]) - 2))
message("number of AAs overlap: ", n_aas_overlap)



targets <- bbbq::get_target_names()
targets <- c("covid", "human", "myco")
n_targets <- length(targets)



t_general <- tibble::tibble(
  target = targets,
  english_name = NA,
  n_proteins = NA,
  n_soluble_proteins = NA,
  n_tmh_proteins = NA,
  n_tmh = NA,
  n_aas_tmh = NA,
  n_aas_non_tmh = NA,
  n_aas = NA
)

# Names
for (i in seq_len(nrow(t_general))) {
  t_general$english_name[i] <- bbbq::get_target_english_name(t_general$target[i])
}

# Number of AAs
for (i in seq_len(nrow(t_general))) {
  target_name <- t_general$target[i]
  message("target_name: ", target_name)
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
  t_general$n_proteins[i] <- nrow(t_proteome)
  t_general$n_aas[i] <- sum(
    stringr::str_length(
      t_proteome$sequence
    )
  )
}

# Number of TMHs according to TMHMM
for (i in seq_len(nrow(t_general))) {
  target_name <- t_general$target[i]
  message("target_name: ", target_name)
  if (target_name == "human") {
    proteome_type <- "representative"
  } else {
    proteome_type <- "full"
  }
  message("proteome_type: ", proteome_type)
  topology <- bbbq::get_topology(
    target_name = target_name,
    proteome_type = proteome_type,
    keep_selenoproteins = FALSE,
    topology_prediction_tool = "tmhmm"
  )
  names(topology) <- c("name", "topology")
  testthat::expect_equal(names(topology), c("name", "topology"))
  t_general$n_soluble_proteins[i] <- sum(
    stringr::str_detect(
      topology$topology, "^[iIoO]+$"
    )
  )
  t_general$n_tmh_proteins[i] <- sum(
    stringr::str_detect(
      topology$topology, paste0("[mM]{", n_aas_overlap, "}"
    )
  )
  t_general$n_aas_tmh[i] <- sum(
    stringr::str_count(topology$topology, paste0("[mM]{", n_aas_overlap, "}")
  )
  t_general$n_aas_non_tmh[i] <- sum(
    stringr::str_count(topology$topology, "[iIoO]")
  )

  t_tmhs <- tmhmm::tally_tmhs(topology)
  t_general$n_tmh[i] <- sum(t_tmhs$n_tmhs)
}

t_general

testthat::expect_equal(
  t_general$n_proteins,
  t_general$n_soluble_proteins + t_general$n_tmh_proteins
)
testthat::expect_equal(
  t_general$n_aas,
  t_general$n_aas_non_tmh + t_general$n_aas_tmh
)
readr::write_csv(x = t_general, file = "general.csv")
