# Creates a figure from 'table_1.csv' and 'table_2.csv'
# to show the measured the number of binders and the
# number of binders that are TMH
#
# Usage:
#
#  Rscript create_pabel_figure.R [percentage]
#
#  * [percentage] is a value
#
#  Rscript create_panel_figure.R 2
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
args <- c("2")
testthat::expect_equal(length(args), 1)
message("Running with arguments {", paste0(args, collapse = ", "), "}")
percentage <- as.numeric(args[1])
message("percentage: ", percentage)

targets <- c("covid", "human", "myco")
message("targets: ", paste0(targets, collapse = ", "))

target_filename <- paste0("fig_f_tmh_", percentage, ".png")
message("target_filename: '", target_filename, "'")
target_filename_grid <- paste0("fig_f_tmh_", percentage, "_grid.png")
message("target_filename_grid: '", target_filename_grid, "'")
target_filename_normalized <- paste0("fig_f_tmh_", percentage, "_normalized.png")
message("target_filename_normalized: '", target_filename_normalized, "'")

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

# Group all proteins
t_tmh_binders <- t_tmh_binders_all %>% dplyr::group_by(target, haplotype) %>%
    dplyr::summarize(
      n_binders = sum(n_binders, na.rm = TRUE),
      n_binders_tmh = sum(n_binders_tmh, na.rm = TRUE),
      n_spots = sum(n_spots, na.rm = TRUE),
      n_spots_tmh = sum(n_spots_tmh, na.rm = TRUE),
      .groups = "drop"
    )


t_tmh_binders$f_tmh <- NA
t_tmh_binders$f_tmh <- t_tmh_binders$n_binders_tmh / t_tmh_binders$n_binders

# Add the MHC class
t_tmh_binders$mhc_class <- NA
t_tmh_binders$mhc_class[
  t_tmh_binders$haplotype %in% bbbq::get_mhc1_haplotypes()
] <- 1
t_tmh_binders$mhc_class[
  t_tmh_binders$haplotype %in% bbbq::get_mhc2_haplotypes()
] <- 2
testthat::expect_equal(0, sum(is.na(t_tmh_binders$mhc_class)))
t_tmh_binders$mhc_class <- as.factor(t_tmh_binders$mhc_class)
t_tmh_binders$mhc_class <- forcats::fct_recode(
  t_tmh_binders$mhc_class,
  "I" = "1",
  "II" = "2"
)

t_tmh_binders$haplotype <- as.factor(t_tmh_binders$haplotype)



t_coincidence <- t_tmh_binders %>% dplyr::group_by(target, mhc_class) %>%
    dplyr::summarize(
      n_spots = mean(n_spots),
      n_spots_tmh = mean(n_spots_tmh),
      .groups = "drop"
    )
t_coincidence$f_tmh <- t_coincidence$n_spots_tmh / t_coincidence$n_spots

f_covid <- t_coincidence$f_tmh[t_coincidence$target == "covid"]
#f_flua  <- t_coincidence$f_tmh[t_coincidence$target == "flua"]
#f_hepa  <- t_coincidence$f_tmh[t_coincidence$target == "hepa"]
#f_hiv   <- t_coincidence$f_tmh[t_coincidence$target == "hiv"]
f_human <- t_coincidence$f_tmh[t_coincidence$target == "human"]
f_myco  <- t_coincidence$f_tmh[t_coincidence$target == "myco"]
#f_polio <- t_coincidence$f_tmh[t_coincidence$target == "polio"]
#f_rhino <- t_coincidence$f_tmh[t_coincidence$target == "rhino"]

t_intercepts <- dplyr::select(t_coincidence, target, mhc_class, f_tmh)

# Humans-only, to compare with other studies
ggplot(
  t_tmh_binders,
  aes(x = haplotype, y = f_tmh)
) +
  ggplot2::geom_col(position = position_dodge(), fill = "#BBBBBB", width = 0.9) +
  ggplot2::facet_grid(
    target ~ mhc_class,
    scales = "free_x",
    space = "free_x",
    labeller = ggplot2::labeller(
      mhc_class = c(I = "MHC-I", II = "MHC-II"),
      target = c(covid = "SARS-CoV-2", human = "Human", myco = "MTb")
    )
  ) +
  ggplot2::xlab("Haplotype") +
  ggplot2::ylab("Epitopes overlapping with TMH") +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 2),
    breaks = seq(0.0, 1.0, by = 0.1),
    minor_breaks = seq(0.0, 1.0, by = 0.1)
  ) +
  bbbq::geom_hline(data = t_intercepts, aes(yintercept = f_tmh), color = "red") +
  bbbq::get_bbbq_theme() +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) + ggplot2::theme(text = element_text(size = 17)) +
  ggsave(
    paste0("fig_f_tmh_", percentage, "_panel.png"),
    width = 7,
    height = 7
  )
