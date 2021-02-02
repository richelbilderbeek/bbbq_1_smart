# Create the 'general.csv' file, containing general info
targets <- bbbq::get_target_names()
targets <- c("covid", "human", "myco")
n_targets <- length(targets)

t_general <- tibble::tibble(
  target = targets,
  english_name = NA,
  n_tmh_tmhmm = seq(1, n_targets),
  # n_tmh_pureseqtm = seq(1, n_targets),
  n_aas_tmh_tmhmm = seq(1, n_targets),
  # n_aas_tmh_pureseqtm = seq(1, n_targets),
  n_aas_non_tmh_tmhmm = seq(1, n_targets),
  # n_aas_non_tmh_pureseqtm = seq(1, n_targets),
  n_aas = seq(1000, n_targets * 1000,by = 1000)
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
  # Sum the non-protein lines
  t_general$n_aas[i] <- sum(
    stringr::str_length(
      t_proteome$sequence
    )
  )
}

shortest_index <- which(t_general$n_aas == min(t_general$n_aas))

# Number of TMHs according to PureseqTM
if (1 == 2) {
  for (i in seq_len(nrow(t_general))) {
    filename <- paste0(t_general$target[i], ".fasta")
    topology <- bbbq::get_topology(
      target_name = t_general$target[i],
      proteome_type = "full",
      keep_selenoproteins = FALSE,
      topology_prediction_tool = "pureseqtmr"
    )
    names(topology) <- c("name", "topology")
    testthat::expect_equal(names(topology), c("name", "topology"))
    t_general$n_aas_tmh_pureseqtm[i] <- sum(
      stringr::str_count(topology$topology, "1")
    )
    t_general$n_aas_non_tmh_pureseqtm[i] <- sum(
      stringr::str_count(topology$topology, "0")
    )
  }
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
  t_general$n_aas_tmh_tmhmm[i] <- sum(
    stringr::str_count(topology$topology, "[mM]")
  )
  t_general$n_aas_non_tmh_tmhmm[i] <- sum(
    stringr::str_count(topology$topology, "[iIoO]")
  )

  t_tmhs <- tmhmm::tally_tmhs(topology)
  t_general$n_tmh_tmhmm[i] <- sum(t_tmhs$n_tmhs)
}
t_general
readr::write_csv(x = t_general, file = "general.csv")
