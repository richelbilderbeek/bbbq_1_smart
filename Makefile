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
  fig_f_tmh_mhc1_2.png \
  fig_f_tmh_mhc2_2.png \
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
	Rscript create_all_counts_per_proteome.R covid 2%

human_2_counts.csv: create_all_counts_per_proteome.R
	Rscript create_all_counts_per_proteome.R human 2%

myco_2_counts.csv: create_all_counts_per_proteome.R
	Rscript create_all_counts_per_proteome.R myco 2%

################################################################################
#
# 2. GENERAL
#
################################################################################

general.csv: create_general.R
	Rscript create_general.R

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

fig_f_tmh_mhc1_2.png: counts_2.csv general.csv create_figure.R
	Rscript create_figure.R mhc1 2

fig_f_tmh_mhc2_2.png: counts_2.csv general.csv create_figure.R
	Rscript create_figure.R mhc2 2

update_packages:
	Rscript -e 'remotes::install_github("richelbilderbeek/peregrine")'
	Rscript -e 'remotes::install_github("richelbilderbeek/mhcnuggetsr")'
	Rscript -e 'remotes::install_github("richelbilderbeek/mhcnpreds")'
	Rscript -e 'remotes::install_github("jtextor/epitope-prediction")'
	Rscript -e 'remotes::install_github("richelbilderbeek/epiprepreds")'
	Rscript -e 'remotes::install_github("richelbilderbeek/netmhc2pan")'
	Rscript -e 'remotes::install_github("richelbilderbeek/nmhc2ppreds")'
	Rscript -e 'remotes::install_github("richelbilderbeek/bbbq", ref = "develop")'

clean:
	rm -f *.png *.latex *.pdf *.fasta *.csv

