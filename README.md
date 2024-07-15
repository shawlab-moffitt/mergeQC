# mergeQC

http://shawlab-moffitt.shinyapps.io/mergeqc

# Installation

## R Dependencies

* `R` - https://cran.r-project.org/src/base/R-4/

|  |  |  |  |
| --- | --- | --- | --- |
| shiny_1.8.1.1 | shinycssloaders_1.0.0 | shinyjqui_0.4.1 | svglite_2.1.3 |
| ggplot2_3.5.1 | dplyr_1.1.4 | reshape2_1.4.4 | plotly_4.10.4 |
| RecordLinkage_0.4-12.4 | DT_0.33 | data.table_1.15.4 |  |

## Deploy
After the required pacakges are installed and the app.R file is in a local repository or your preference, you can run the app using the command `runApp("[path to app.R file]")`


# Input Data

## Expression Data

*  Tab-delimited data table
  *  First column consists of feature names
  *  Following columns are of numeric expression values, with the column names referring to the sample name

## Meta Data

* Tab-delimited data table
  *  The first column should consist of sample names that are found in the expression data column names
  *  The following columns can be additional information on these samples 


