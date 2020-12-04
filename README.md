# mrbrowse

## DESCRIPTION

Type: R shiny application
Title: Interactive Browser for High Dimensional Mendelian Randomization Studies
Version: 1.3
Author: Michael C Sachs <michael.sachs@ki.se>
Maintainer: Michael C Sachs <michael.sachs@ki.se>
Description: In the SCALLOP CVD-I project researchers at the Karolinska Institute 
    have performed Mendelian Randomization (MR) analyses of 100s of proteins against 
    100s of clinical outcomes using genetic data (attached manuscript). This work 
    will continue with analysis in other SCALLOP projects as well as the KARMA 
    project. To supplement publication and enable easier browsing of MR results 
    there is a need for a graphical interface to easily search and visualise 
    results from MR analyses.
License: MIT + file LICENSE.md
Encoding: UTF-8
Depends: 
    shiny,
    ggplot2,
    data.table,
	openxlsx,
    R (>= 4.0)

	
## Usage

### Run the shiny app

Clone and unzip the repository, and run the following in R from the mrbrowse working directory: 

    library(shiny)
	runApp("shinyapp")
	
	
### Update or change the data

1. The raw datasets reside in the txtfiles subdirectory
2. Add or remove files in the relevant subdirectory of that (pQTLs or CVDI), ensuring that you use the same file format and naming structure.
3. Run the `parse_datasets.R` script
4. Run the shiny app