# pozniak_melanoma
This repository contains an R Shiny app and Rmarkdown workflow generating
interactive visualizations from the melanoma data from Pozniak, et al., 2022
paper. Both the R Shiny app and Rmarkdown workflow apply GSVA to create a gene
set score for each gene sets, which is used to produce visualizations including:
* Violin plots
* Correlation plots
* Heatmaps

This workflow is primarily designed with the intention to correlate two gene
sets. Here we load such a data file with two toy gene sets. It is possible to 
upload a custom user-generated gene set as long as the file is a CSV or TSV file
and the gene sets are the first and second columns. This workflow ignores
additional columns, and these columns may be named whatever the end user wishes.
