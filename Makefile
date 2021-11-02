#
# Usage:
#
# Create the data on Peregrine
#
#   make peregrine
#
# Create the results locally, assume data is there
#
#   make results
#
all: \
  covid_2_counts.csv \
  myco_2_counts.csv \
  human_2_counts.csv \
  general.csv \
  counts_2.csv \
  table_f_tmh_2.latex \
  fig_f_tmh_mhc1_2.tiff \
  fig_f_tmh_mhc2_2.tiff \
  fig_rel_presentation.tiff \
  fig_f_tmh_2_human_mhc1.tiff \
  table_ic50_binders_mhc1_2.latex \
  table_ic50_binders_mhc2_2.latex \
  table_tmh_binders_mhc1_2.latex \
  table_tmh_binders_mhc2_2.latex

################################################################################
#
# 1. COUNTS PER PROTEIN
#
################################################################################

covid_2_counts.csv: create_all_counts_per_proteome.R
	Rscript -e 'install.packages("forcats")'
	Rscript -e 'install.packages("tidyr")'
	Rscript -e 'remotes::install_github("richelbilderbeek/mhcnuggetsr")'
	Rscript -e 'remotes::install_github("richelbilderbeek/mhcnuggetsrinstall")'
	Rscript -e 'if (!mhcnuggetsr::is_mhcnuggets_installed()) mhcnuggetsrinstall::install_mhcnuggets()'
	Rscript -e 'remotes::install_github("richelbilderbeek/mhcnpreds")'
	Rscript -e 'remotes::install_github("richelbilderbeek/nmhc2ppreds")'
	Rscript -e 'remotes::install_github("richelbilderbeek/tmhmm")'
	Rscript -e 'remotes::install_github("richelbilderbeek/epiprepreds")'
	Rscript -e 'remotes::install_github("richelbilderbeek/pureseqtmr")'
	Rscript -e 'remotes::install_bioc("Biostrings")'
	Rscript -e 'remotes::install_github("richelbilderbeek/bbbq")'
	Rscript create_all_counts_per_proteome.R covid 2% 2AA

human_2_counts.csv: create_all_counts_per_proteome.R
	Rscript create_all_counts_per_proteome.R human 2% 2AA

myco_2_counts.csv: create_all_counts_per_proteome.R
	Rscript create_all_counts_per_proteome.R myco 2% 2AA

################################################################################
#
# 2. GENERAL
#
################################################################################

general.csv: create_general.R
	Rscript create_general.R 2AA

counts_2.csv: merge_all_counts_per_proteome.R
	Rscript merge_all_counts_per_proteome.R 2

################################################################################
# Create the CSV tables for the binders
################################################################################

table_tmh_binders_mhc1_2.latex: counts_2.csv create_table_tmh_binders_mhc.R
	Rscript create_table_tmh_binders_mhc.R mhc1 2

table_tmh_binders_mhc2_2.latex: counts_2.csv create_table_tmh_binders_mhc.R
	Rscript create_table_tmh_binders_mhc.R mhc2 2

################################################################################
# Create all LaTeX tables
################################################################################

# Easy and general table
table_ic50_binders_mhc1_2.latex: create_table_ic50_binders.R
	Rscript create_table_ic50_binders.R mhc1 2%

table_ic50_binders_mhc2_2.latex: create_table_ic50_binders.R
	Rscript create_table_ic50_binders.R mhc2 2%

table_f_tmh_2.latex: create_table_f_tmh.R
	Rscript create_table_f_tmh.R 2

################################################################################
# Create the figures
################################################################################

fig_f_tmh_mhc1_2.tiff: counts_2.csv general.csv create_figure.R
	Rscript create_figure.R mhc1 2

fig_f_tmh_mhc2_2.tiff: counts_2.csv general.csv create_figure.R
	Rscript create_figure.R mhc2 2

fig_rel_presentation.tiff: counts_2.csv general.csv create_fig_rel_presentation.R
	Rscript create_fig_rel_presentation.R

fig_f_tmh_2_human_mhc1.tiff: counts_2.csv general.csv create_panel_figure.R
	Rscript create_panel_figure.R

################################################################################
# Misc
################################################################################

update_packages:
	Rscript -e 'remotes::install_github("richelbilderbeek/peregrine")'
	Rscript -e 'remotes::install_github("richelbilderbeek/mhcnuggetsr")'
	Rscript -e 'remotes::install_github("richelbilderbeek/mhcnpreds")'
	Rscript -e 'remotes::install_github("jtextor/epitope-prediction")'
	Rscript -e 'remotes::install_github("richelbilderbeek/epiprepreds")'
	Rscript -e 'remotes::install_github("richelbilderbeek/netmhc2pan")'
	Rscript -e 'remotes::install_github("richelbilderbeek/nmhc2ppreds")'
	Rscript -e 'remotes::install_github("richelbilderbeek/pureseqtmr")'
	Rscript -e 'remotes::install_github("richelbilderbeek/bbbq", ref = "develop")'

clean:
	rm -f *.png *.latex *.pdf *.fasta *.csv

