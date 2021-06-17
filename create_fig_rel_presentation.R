# Plot the relative over-presentation of TMH-derived epitopes
# between MCH-I and MHC-II for the different targets
#
# Usage:
#
#  Rscript create_fig_rel_presentation.R
#
#  * [MHC] is both 'mhc1' and 'mhc2'
#  * [percentage] is 2
#
#  Rscript create_fig_rel_presentation.R
#
suppressPackageStartupMessages({
  library(dplyr, warn.conflicts = FALSE)
  library(testthat, warn.conflicts = FALSE)
  library(ggplot2, quietly = TRUE)
})

percentage <- 2

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

# Add the MHC class
t_tmh_binders_all$mhc_class <- NA
t_tmh_binders_all$mhc_class[
  t_tmh_binders_all$haplotype %in% bbbq::get_mhc1_haplotypes()
] <- 1
t_tmh_binders_all$mhc_class[
  t_tmh_binders_all$haplotype %in% bbbq::get_mhc2_haplotypes()
] <- 2
testthat::expect_equal(0, sum(is.na(t_tmh_binders_all$mhc_class)))
t_tmh_binders_all$mhc_class <- as.factor(t_tmh_binders_all$mhc_class)
t_tmh_binders_all$mhc_class <- forcats::fct_recode(
  t_tmh_binders_all$mhc_class,
  "I" = "1",
  "II" = "2"
)


# Group all proteins
t_tmh_binders <- t_tmh_binders_all %>%
    dplyr::group_by(target, mhc_class, haplotype) %>%
    dplyr::summarize(
      n_binders = sum(n_binders, na.rm = TRUE),
      n_binders_tmh = sum(n_binders_tmh, na.rm = TRUE),
      n_spots = sum(n_spots, na.rm = TRUE),
      n_spots_tmh = sum(n_spots_tmh, na.rm = TRUE),
      .groups = "drop"
    )
t_tmh_binders$f_tmh_observed <- t_tmh_binders$n_binders_tmh / t_tmh_binders$n_binders
t_tmh_binders$f_tmh_chance <- t_tmh_binders$n_spots_tmh / t_tmh_binders$n_spots
t_tmh_binders$haplotype <- as.factor(t_tmh_binders$haplotype)
t_tmh_binders$normalized_f_tmh <- t_tmh_binders$f_tmh_observed / t_tmh_binders$f_tmh_chance

target_to_english_lut <- c(
  covid = "SARS-CoV-2",
  human = "Human",
  myco = "MTb"
)

p1 <- ggplot(t_tmh_binders, aes(x = haplotype, y = normalized_f_tmh, fill = mhc_class)) +
  geom_col(position = position_dodge(), color = "#000000") +
  ggplot2::facet_grid(
    target ~ . ,
    labeller = ggplot2::labeller(target = target_to_english_lut)
  ) +
  xlab(paste0("Haplotype")) +
  ylab("Normalized epitopes overlapping \nwith transmembrane helix") +
  scale_y_continuous() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(aes(yintercept = 1), lty = "dashed") +
  labs(
    title = "Normalized % epitopes that overlap with TMH per haplotype",
    caption = glue::glue(
      "Dashed line: normalized expected percentage of epitopes ",
      "that have one residue overlapping with a TMH"
    ),
    fill = "MHC class"
  )
p1

p1 + ggplot2::ggsave(
  filename = "fig_rel_presentation_per_haplotype.png",
  width = 7,
  height = 7
)

p1 + ggplot2::ggsave(
  filename = "fig_rel_presentation_per_haplotype.tiff",
  width = 7,
  height = 7
)

p1 + ggplot2::scale_color_brewer(palette = "Greys") +
  ggplot2::scale_fill_brewer(palette = "Greys") ; ggplot2::ggsave(
  filename = "fig_rel_presentation_per_haplotype_bw.png",
  width = 7,
  height = 7
)

p1 + ggplot2::scale_color_brewer(palette = "Greys") +
  ggplot2::scale_fill_brewer(palette = "Greys") ; ggplot2::ggsave(
  filename = "fig_rel_presentation_per_haplotype_bw.tiff",
  width = 7,
  height = 7
)



t_per_mhc_class <- t_tmh_binders %>%
  dplyr::group_by(target, mhc_class) %>%
  dplyr::summarise(
    mean_f_tmh_observed = mean(normalized_f_tmh),
    f_tmh_observed_se = stats::sd(normalized_f_tmh)/sqrt(dplyr::n()),
    .groups = "drop"
  )

p2 <- ggplot(t_per_mhc_class, aes(x = mhc_class, y = mean_f_tmh_observed)) +
  geom_col(position = position_dodge(), fill = "#BBBBBB") +
  geom_errorbar(
    aes(
      x = mhc_class,
      ymin = mean_f_tmh_observed - f_tmh_observed_se,
      ymax = mean_f_tmh_observed + f_tmh_observed_se
    ),
    width = 0.4
  ) +
  # ggplot2::geom_text(
  #   ggplot2::aes(
  #     y = 0.0,
  #     label = format(mean_f_tmh_observed, digits = 3)
  #   ),
  #   vjust = -0.5
  # ) +
  ggplot2::facet_grid(
    target ~ .,
    labeller = ggplot2::labeller(target = target_to_english_lut)
  ) +
  xlab("MHC class") +
  ylab("Normalized epitopes overlapping with TMH") +
  scale_y_continuous() +
  geom_hline(aes(yintercept = 1), lty = "dashed") +
  labs(
    title = "Normalized % epitopes that overlap with TMH per MHC class",
    fill = "MHC class"
  ) + ggplot2::theme_bw() +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    strip.background = element_rect(colour="white", fill="#FFFFFF")
  )
p2

p2 + ggplot2::ggsave(
  filename = "fig_rel_presentation.png",
  width = 7,
  height = 7
)

p2 + ggplot2::ggsave(
  filename = "fig_rel_presentation.tiff",
  width = 7,
  height = 7
)

