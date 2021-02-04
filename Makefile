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
  human_2_counts.csv \
  myco_2_counts.csv \
  general.csv \
  counts_2.csv \
  fig_f_tmh_mhc1_2.png \
  fig_f_tmh_mhc2_2.png

#haplotypes_lut.csv \
#     covid_proteins_lut.csv flua_proteins_lut.csv hepa_proteins_lut.csv \
#     hiv_proteins_lut.csv human_proteins_lut.csv \
#     myco_proteins_lut.csv polio_proteins_lut.csv rhino_proteins_lut.csv \
#     flua_h1_counts.csv \
#     hiv_h1_counts.csv \
#     human_h1_counts.csv \
#     rhino_h1_counts.csv \
#     counts.csv \
#     table_tmh_binders_mhc1.latex table_tmh_binders_mhc2.latex \
#     table_ic50_binders.latex \
#     table_f_tmh.latex \
#     fig_f_tmh_mhc1.png fig_f_tmh_mhc2.png general.csv create_figure.R

################################################################################
#
# 1. COUNTS PER PROTEIN
#
################################################################################

covid_2_counts.csv:
	Rscript create_all_counts_per_proteome.R covid 2

human_2_counts.csv:
	Rscript create_all_counts_per_proteome.R human 2

myco_2_counts.csv:
	Rscript create_all_counts_per_proteome.R myco 2

################################################################################
#
# 2. GENERAL
#
################################################################################

general.csv:
	Rscript create_general.R

counts_2.csv:
	Rscript merge_all_counts_per_proteome.R 2

################################################################################
# Create the CSV tables for the binders
################################################################################

table_tmh_binders_mhc1.latex: counts.csv
	Rscript create_table_tmh_binders_mhc.R mhc1

table_tmh_binders_mhc2.latex: counts.csv
	Rscript create_table_tmh_binders_mhc.R mhc2

################################################################################
# Create all LaTeX tables
################################################################################

# Easy and general table
table_ic50_binders.latex: haplotypes_lut.csv
	Rscript create_table_ic50_binders.R

table_f_tmh.latex:
	Rscript create_table_f_tmh.R

################################################################################
# Create the figures
################################################################################

fig_f_tmh_mhc1_2.png: counts_2.csv general.csv
	Rscript create_figure.R mhc1 2

fig_f_tmh_mhc2_2.png: counts_2.csv general.csv
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

