# pozniak_melanoma
This repository contains an R Shiny app and Rmarkdown workflow generating
interactive visualizations from the human and melanoma data from Pozniak, et al.,
2022 paper. Both the R Shiny app and Rmarkdown workflow apply GSVA to create a gene
set score for two gene sets of the users choice, which is used to produce
visualizations including:
* Violin plots
* Correlation plots
* Heatmaps

This workflow is primarily designed with the intention to correlate two gene
sets. Here we load such a data file with two toy gene sets. It is possible to 
upload a custom user-generated gene set as long as the file is a CSV or TSV file
and the gene sets are the first and second columns. This workflow ignores
additional columns, and these columns may be named whatever the end user wishes.

Two R Shiny apps are available, one Human and one Mouse. Due to the sample
composition, fewer violin and correlation plots are available for the Mouse data.

If users wish to run these Shiny apps on their own local machines, they are
advised to change the path of the data files to where this repo resides on their
machine.
