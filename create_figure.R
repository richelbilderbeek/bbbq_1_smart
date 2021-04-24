# Creates a figure from 'table_1.csv' or 'table_2.csv'
# to show the measured the number of binders and the
# number of binders that are TMH
#
# Usage:
#
#  Rscript create_figure.R [MHC] [percentage]
#
#  * [MHC] is either 'mhc1' or 'mhc2'
#  * [percentage] is a value
#
#  Rscript create_figure.R mhc1 5
#
#
#
#
suppressPackageStartupMessages({
  library(dplyr, warn.conflicts = FALSE)
  library(testthat, warn.conflicts = FALSE)
  library(ggplot2, quietly = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)
if (1 == 2) {
  args <- c("mhc1", "2")
  args <- c("mhc2", "2")
}
testthat::expect_equal(length(args), 2)
message("Running with arguments {", paste0(args, collapse = ", "), "}")
mhc <- args[1]
message("mhc: ", mhc)
testthat::expect_equal(4, stringr::str_length(mhc))
percentage <- as.numeric(args[2])
message("percentage: ", percentage)

mhc_class <- as.numeric(stringr::str_sub(mhc, 4, 4))
message("mhc_class: ", mhc_class)
haplotypes <- NA
if (mhc_class == 1) {
  haplotypes <- bbbq::get_mhc1_haplotypes()
  peptide_length <- 9
} else {
  testthat::expect_equal(2, mhc_class)
  haplotypes <- bbbq::get_mhc2_haplotypes()
  peptide_length <- 14
}
message("haplotypes: ", paste0(haplotypes, collapse = ", "))
message("peptide_length: ", peptide_length)


targets <- c("covid", "human", "myco")
message("targets: ", paste0(targets, collapse = ", "))

general_filename <- "general.csv"
message("general_filename: '", general_filename, "'")
testthat::expect_true(file.exists(general_filename))
t_general <- readr::read_csv(
  general_filename,
  col_types = readr::cols(
    target = readr::col_character(),
    english_name = readr::col_character(),
    n_tmh = readr::col_double(),
    n_aas = readr::col_double()
  )
)

table_filename <- paste0("counts_", percentage, ".csv")
message("table_filename: '", table_filename, "'")
testthat::expect_true(file.exists(table_filename))
t_tmh_binders_all <- readr::read_csv(
  table_filename,
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

# Only keep the desired MHC class
t_tmh_binders_mhc <- t_tmh_binders_all %>% filter(haplotype %in% haplotypes)

# Group all proteins
t_tmh_binders <- t_tmh_binders_mhc %>% dplyr::group_by(target, haplotype) %>%
    dplyr::summarize(
      n_binders = sum(n_binders, na.rm = TRUE),
      n_binders_tmh = sum(n_binders_tmh, na.rm = TRUE),
      n_spots = sum(n_spots, na.rm = TRUE),
      n_spots_tmh = sum(n_spots_tmh, na.rm = TRUE),
      .groups = "drop"
    )


t_tmh_binders$f_tmh <- NA
t_tmh_binders$f_tmh <- t_tmh_binders$n_binders_tmh / t_tmh_binders$n_binders
t_tmh_binders$haplotype <- as.factor(t_tmh_binders$haplotype)

t_coincidence <- t_tmh_binders %>% dplyr::group_by(target) %>%
    dplyr::summarize(
      n_spots = mean(n_spots),
      n_spots_tmh = mean(n_spots_tmh),
      .groups = "drop"
    )
t_coincidence$f_tmh <- t_coincidence$n_spots_tmh / t_coincidence$n_spots

f_covid <- t_coincidence$f_tmh[t_coincidence$target == "covid"]
f_human <- t_coincidence$f_tmh[t_coincidence$target == "human"]
f_myco  <- t_coincidence$f_tmh[t_coincidence$target == "myco"]

roman_mhc_class <- as.character(as.roman(mhc_class))
testthat::expect_true(
  mhc_class == 1 && roman_mhc_class == "I" ||
  mhc_class == 2 && roman_mhc_class == "II"
)

# Humans-only, to compare with other studies
p <- ggplot(
  t_tmh_binders %>% dplyr::filter(target == "human"),
  aes(x = haplotype, y = f_tmh)
) +
  geom_col(position = position_dodge(), color = "#000000") +
  xlab(paste0("MHC-", roman_mhc_class, " HLA haplotype")) +
  ylab("Epitopes overlapping \nwith transmembrane helix") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 2),
    breaks = seq(0.0, 1.0, by = 0.1),
    minor_breaks = seq(0.0, 1.0, by = 0.1)
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(aes(yintercept = f_human), lty = "dashed") +
  labs(
    title = "% epitopes that overlap with TMH per haplotype",
    caption = paste0(
      "Dashed horizontal line: % ", peptide_length,
      "-mers that overlaps with TMH: ",
      formatC(100.0 * mean(f_human), digits = 3),"%"
    )
  )
p
p + ggsave(
    paste0("fig_f_tmh_mhc", mhc_class, "_", percentage, "_human.png"),
    width = 7,
    height = 7
)
p + ggsave(
    paste0("fig_f_tmh_mhc", mhc_class, "_", percentage, "_human.tiff"),
    width = 7,
    height = 7
)


caption_text <- paste0(
  "Horizontal lines: % ", peptide_length, "-mers that overlaps with TMH in ",
  "SARS-Cov2 (",     formatC(100.0 * mean(f_covid), digits = 3),"%), ",
  "humans (",        formatC(100.0 * mean(f_human), digits = 3),"%), ",
  "Mycobacterium (", formatC(100.0 * mean(f_myco ), digits = 3),"%), "
)
p <- ggplot(t_tmh_binders, aes(x = haplotype, y = f_tmh, fill = target)) +
  geom_col(position = position_dodge(), color = "#000000") +
  xlab(paste0("MHC-", roman_mhc_class, " HLA haplotype")) +
  ylab("Epitopes overlapping \nwith transmembrane helix") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 2),
    breaks = seq(0.0, 1.0, by = 0.1),
    minor_breaks = seq(0.0, 1.0, by = 0.1)
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(data = t_coincidence, aes(yintercept = f_tmh, lty = target)) +
  labs(
    title = "% epitopes that overlap with TMH per haplotype",
    caption = caption_text
  )

p + ggsave(
  paste0("fig_f_tmh_mhc", mhc_class, "_", percentage, ".png"),
  width = 7,
  height = 7
)
p + ggsave(
  paste0("fig_f_tmh_mhc", mhc_class, "_", percentage, ".tiff"),
  width = 7,
  height = 7
)


# Facet labels
facet_labels <- paste0(
  t_general$english_name, "\n",
  "TMHs: ", t_general$n_tmh, "\n",
  t_general$n_aas, " AAs"
)
names(facet_labels) <- t_general$target


p_facet <- p + facet_grid(
  target ~ ., scales = "free",
  labeller = ggplot2::as_labeller(facet_labels)
) + ggplot2::theme(strip.text.y.right = ggplot2::element_text(angle = 0)) +
  ggplot2::theme(legend.position = "none")

p_facet + ggsave(
  paste0("fig_f_tmh_mhc", mhc_class, "_", percentage, "_grid.png"),
  width = 7, height = 14
)
p_facet + ggsave(
  paste0("fig_f_tmh_mhc", mhc_class, "_", percentage, "_grid.tiff"),
  width = 7, height = 14
)


# Normalize
t_tmh_binders$coincidence <- NA
for (i in seq_len(nrow(t_tmh_binders))) {
  target <- t_tmh_binders$target[i]
  coincidence <- NA
  if (target == "covid") coincidence <- f_covid
  #else if (target == "flua") coincidence <- f_flua
  #else if (target == "hepa") coincidence <- f_hepa
  #else if (target == "hiv") coincidence <- f_hiv
  else if (target == "human") coincidence <- f_human
  else if (target == "myco") coincidence <- f_myco
  #else if (target == "polio") coincidence <- f_polio
  #else if (target == "rhino") coincidence <- f_rhino
  else stop("?")
  t_tmh_binders$coincidence[i] <- coincidence

}
t_tmh_binders$normalized_f_tmh <- t_tmh_binders$f_tmh / t_tmh_binders$coincidence

p <- ggplot(t_tmh_binders, aes(x = haplotype, y = normalized_f_tmh, fill = target)) +
  geom_col(position = position_dodge(), color = "#000000") +
  xlab(paste0("MHC-", roman_mhc_class, " HLA haplotype")) +
  ylab("Normalized epitopes overlapping \nwith transmembrane helix") +
  scale_y_continuous() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(aes(yintercept = 1), lty = "dashed") +
  labs(
    title = "Normalized % epitopes that overlap with TMH per haplotype",
    caption = glue::glue(
      "Dashed line: normalized expected percentage of epitopes ",
      "that have one residue overlapping with a TMH"
    )
  )

p + ggsave(
  paste0("fig_f_tmh_mhc", mhc_class, "_", percentage, "_normalized.png"),
  width = 7, height = 7
)
p + ggsave(
  paste0("fig_f_tmh_mhc", mhc_class, "_", percentage, "_normalized.tiff"),
  width = 7, height = 7
)


