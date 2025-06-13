#' ---
#' title: "Puma space use and dispersal in the tropics: Create R data package"
#' author: Bernardo Niebuhr
#' date: ""
#' output: 
#'   rmarkdown::github_document: default
#'   pdf_document: default
#' geometry: margin=2cm
#' ---

# --------------- label=setup, warning=FALSE, message=FALSE, echo=TRUE

# Load packages
library(knitr)
library(devtools)

# Print options for this document
options(width = 165)
opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # root project folder
opts_chunk$set(error = F, message = F, warning = F, cache = F, eval = F, echo = T)

# ----------- label=initialize

#' # 1) Create the package within the folder

# create the package structure
usethis::create_package("PardasIPC")

#' Now that the structure is created, we should re-open this script within the package project and folder,
#' then it is possible to continue.
#' 
#' # 2) Create README
#' 
#' Here we can create a README.md file and edit what to put in this file. This information will be shown
#' in the first page in the Github repository where the package will be stored.

# Begin readme
usethis::use_readme_md()

#' We could also define a license for the data, but in this context of private data it is better to define
#' such a license by hand and add it to the repository.
#' 
#' We should also edit the `DESCRIPTION` file.
#' 
#' # 3) Add data
#' 
#'  Now we can copy our already prepared `.rda` file to the data folder. If it does not exists, we can create it
#'  by hand. We can also store the raw version of the datasets.

# add raw data
usethis::use_data_raw()

#' This function creates the structure of the data-raw folder. I have copied there the raw data files, and
#' added the cleaning script to prepare the actual data file in there.
#' 
#' Now we can create the folder for the data that will be used.

# create data for use in the Rpackage
load("../data/movement_data_pumas_tiete_legado.rda")
pumas <- mov_data
class(pumas)

usethis::use_data(pumas, overwrite = TRUE)

#' Now we must document our data by hand and put info in the `R/` folder.
#' 
#' # 4) Add packages upon which the package depends on

# Necessary packages
usethis::use_package("dplyr", type = "Imports")
# Suggested packages
usethis::use_package("amt", type = "Suggests")

#' # 5) Build the package
#' 
#' Now we can compile/build the package. 