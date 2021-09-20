# scViewer-Lite
Light version of [SeuratViewer.](https://github.com/Morriseylab/SeuratViewer). This is a R Shiny website for viewing single cell RNA-Seq data analysed using [Seurat.](https://satijalab.org/seurat/) (version 3). You can view the Seurat repository [here](https://github.com/satijalab/seurat)

## Introduction
scViewer Lite reads in the seurat object as an RDS file and enables users to view and download UMAPs and Gene Expression plots of their single cell RNAseq data

## Requirements
- R (version > 3.5)
- RStudio Server
- Shiny Server (if you need to host it online)

If you need help installing the above or getting started, refer to [this](https://deanattali.com/2015/05/09/setup-rstudio-shiny-server-digital-ocean/#install-r)
 
## Installation
For Linux, run the following commands in terminal 
```
sudo apt-get install libcurl4-openssl-dev libssl-dev
sudo apt-get install xorg libx11-dev mesa-common-dev libglu1-mesa-dev
sudo apt-get install libxml2-dev
sudo apt-get install libftgl2 freetype2-demos libfreetype6-dev
sudo apt-get install libhdf5-dev
sudo apt-get install r-cran-rcppeigen
```
Run the following commands in R to install all required packages
```
install.packages(c("devtools","shiny","shinydashboard","shinyjs","shinyBS","RColorBrewer","reshape2","ggplot2",
                   "dplyr","tidyr","openssl","httr","plotly","htmlwidgets","DT","Seurat","cowplot",
                    "data.table","NMF","tibble"))

#Install packages from bioconductor
install.packages("BiocManager")
BiocManager::install(c("biomaRt","Biobase"))

##This package contains helper functions 
require(devtools)
install_github("Morriseylab/scExtras")
```
For linux users, other R dependencies include
- RcppEigen
- lme4
- flexmix


## Creating Input data

The input data is created using this [script](https://github.com/bapoorva/scripts/blob/master/RunSeurat_WT.R). 

### Adding your dataset

Add your data to the param.csv file and move it to the data directory. Please note that the data directory must be in the same location as your server.R and ui.R files. The param.csv file should also be saved in the data directory as the RDS files.

### Shiny App


<img width="1821" alt="Screen Shot 2021-09-20 at 9 19 42 AM" src="https://user-images.githubusercontent.com/43073258/134009248-198f58f7-2964-4edb-b88b-5220bc45008e.png">
