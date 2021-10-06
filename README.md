# bbbq_1_smart

Branch   |[![GitHub Actions logo](pics/GitHubActions.png)](https://github.com/richelbilderbeek/bbbq_1_smart/actions)|[![Codecov logo](pics/Codecov.png)](https://www.codecov.io)
---------|-----------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------
`master` |![make](https://github.com/richelbilderbeek/bbbq_1_smart/workflows/make/badge.svg?branch=master)   |[![codecov.io](https://codecov.io/github/richelbilderbeek/bbbq_1_smart/coverage.svg?branch=master)](https://codecov.io/github/richelbilderbeek/bbbq_1_smart/branch/master)
`develop`|![make](https://github.com/richelbilderbeek/bbbq_1_smart/workflows/make/badge.svg?branch=develop)  |[![codecov.io](https://codecov.io/github/richelbilderbeek/bbbq_1_smart/coverage.svg?branch=develop)](https://codecov.io/github/richelbilderbeek/bbbq_1_smart/branch/develop)

A pipeline to answer the first sub-question of the 
Bianchi, Bilderbeek and Bogaart Question.

 * [Full article](https://github.com/richelbilderbeek/bbbq_article)

## Build

```
make
```

## File structure

### `[target]_[percentage]_counts.csv`

`haplotype`   |`name`      |`n_binders`|`n_binders_tmh`|`n_spots`|`n_spots_tmh`
--------------|------------|-----------|---------------|---------|-------------
h1            |p1          |11         |5              |100      |20
h1            |p2          |12         |4              |10       |2

Note that:

 * `n_spots` and `n_spots_tmh` can vary, due to MHC class-dependent epitope lengths.
 * we only keep track of the membrane proteins, `n_spots_tmh` will thus be
   always greater than zero

```
Rscript create_all_counts_per_proteome.R [target] [percentage]
Rscript create_all_counts_per_proteome.R human 5
```

### `counts_[percentage].csv`

`target`|`haplotype_id`|`protein_id`|`n_binders`|`n_binders_tmh`|`n_spots`|`n_spots_tmh`
--------|--------------|------------|-----------|---------------|---------|-------------
covid   |h1            |p1          |11         |5              |100      |20
covid   |h1            |p2          |12         |6              |101      |20

Note that:

 * we only keep track of the membrane proteins, `n_spots_tmh` will thus be
   always greater than zero


```
Rscript merge_all_counts.R [percentage]
Rscript merge_all_counts.R 5
```

### `table_tmh_binders_mhc[mhc_class].csv`

Pretty-printed version

`haplotype`|`covid`      |`human`
-----------|-------------|-------------
HLA-A*01:01| 38.46 (5/13)| 25.00 (5/20)
HLA-B*39:01| 100.00 (2/2)|58.33 (14/24)
HLA-B*40:02|  55.56 (5/9)| 29.17 (7/24)

Note that:

 * up so far, we only kept track of the membrane proteins. 
   Here we correct the percentages to include the full proteome again

```
Rscript create_table_tmh_binders_mhc.R mhc1
Rscript create_table_tmh_binders_mhc.R mhc2
```

## Figures

Note that:

 * up so far, we only kept track of the membrane proteins. 
   Here we correct the percentages to include the full proteome again

![](fig_f_tmh_mhc1_grid.png)

![](fig_f_tmh_mhc1_normalized.png)

![](fig_f_tmh_mhc1.png)

![](fig_f_tmh_mhc2_grid.png)

![](fig_f_tmh_mhc2_normalized.png)

![](fig_f_tmh_mhc2.png)

