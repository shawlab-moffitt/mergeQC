# mergeQC

We present an R Shiny application that allows for users to merge multiple datasets of expression and meta data, as well as the ability to perform a brief check of the quality between the datasets.

This application can be downloaded and ran locally or through out public mergeQC page below.

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

* Example data is available to download through the GitHub repo or through a link on the app interface.

## Expression Data

*  Tab-delimited data table
   *  First column consists of feature names
   *  Following columns are of numeric expression values, with the column names referring to the sample name

## Meta Data

* Tab-delimited data table
  *  The first column should consist of sample names that are found in the expression data column names
  *  The following columns can be additional information on these samples 

# Load data into the app

Use CTRL + left click to select the proper files for each input section

![alt text](https://github.com/shawlab-moffitt/mergeQC/blob/main/Example_Images/Load_data.PNG?raw=true)

![alt text](https://github.com/shawlab-moffitt/mergeQC/blob/main/Example_Images/Loaded_data.PNG?raw=true)

## Data Processing

### Expression Data

* Expression matrices are merged based on similar features between input data. Users can view and download both an inner merged matrix or full merged matrix. The inner merged matrix will only contain features that were similar throughout all of the matrices input. The full merged matrix will contain all genes from the input matrices, and if there is a gene not found in one matrix, that value will be NA.

### Meta Data

* The meta data will be merged based on similar sample names found in the first column and similar column names. If similarities are not found the data will be added as additional rows and NA values in columns that do not match.

### Notes

* At the top of the data preview page a summary of the merged data will be displayed.


# Merged Data QC

## Average Expression Density

* Users can view the density of the average expression over all genes in each data set, as well as overlap datasets by selecting multiple rows.
* This can help users check their data and see if any transformations need to be done to ensure it is comparable to the other datasets.
* Users can select datasets on the sidebar panel and the datasets selected will be log2 transformed and the data will be updated on all tabs.


![alt text](https://github.com/shawlab-moffitt/mergeQC/blob/main/Example_Images/dataqc_log_trans.png?raw=true)

![alt text](https://github.com/shawlab-moffitt/mergeQC/blob/main/Example_Images/density.PNG?raw=true)

## Average Expression Comparison

* Users can select two datasets to compare expression data through an interactive scatter plot.
* Upon hover users can see the gene and expression value from each dataset.
* The red points have a higher expression in the y-axis dataset, and the bluse have a higher expression int eh x-axis dataset.
* User can also sear for and highlight specifc genes on the sidebar.

![alt text](https://github.com/shawlab-moffitt/mergeQC/blob/main/Example_Images/avg_expr.PNG?raw=true)



















