#' ---
#' title: "Puma space use and dispersal in the tropics: Load and organize data"
#' author: Bernardo Niebuhr
#' date: ""
#' geometry: margin=2cm
#' output:
#'   pdf_document: default
#'   github_document: default
#' ---

#' In this document, we provide the code to organize and prepare a single 
#' dataset with GPS data from 14 puma individuals monitored in Southern Brazil.
#' We also show how the dataset might be exported from R into different formats,
#' to be visualized and analyzed in different GIS and statistical software.
#' Furthermore, we provide information about the availability of the GPS data
#' in MoveBank studies and in a private R package called `PardasIPC`. This 
#' document supplements the paper: Niebuhr et al. *Puma space use and dispersal 
#' in the tropics: filling a gap to connect individuals to populations*.

# --------------- label=setup, warning=FALSE, message=FALSE, error=FALSE, echo=FALSE

# Load packages
if(!require(install.load)) install.packages("install.load"); library(install.load)

install.load::install_load("knitr", "rprojroot", "NinaR")
# Print options for this document
options(width = 165)
opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # root project folder
opts_chunk$set(error = F, 
               message = F, 
               warning = F, 
               cache = F, 
               echo = T, 
               eval = T, 
               results = T,
               tidy = "styler",
               dev = c("png", "pdf"),
               dpi = 600,
               fig.path = "figure/",
               fig.align="center",
               fig.pos = 'H')

#' # Load packages and setup

library(install.load)
install.load::install_load("dplyr", "purrr", "lubridate")
install.load::install_load("sf")

# Clean everything before beginning
rm(list = ls())

#' # Loading and organizing data
#' 
#' First of all, we loaded and organized the two sources of movement data: (i) cougar data from
#' the Pardas do Tietê project (n = 13) and (ii) cougar data from Legado das Aguas project (n = 1).
#' Our final dataset then comprises 14 puma individuals and about 56,000 GPS fixes. 
#' Below we present the code for loading and organizing the data.

# --------------- label=load_data

#' ## Loading data
#' 
#' First, we loaded the data and selected which individuals were considered for the analyses. This included
#' all _Puma concolor_ individuals. The dispersal status of the individuals was assessed previously,
#' by visual inspection. Although all individuals were classified as if they dispersed or not 
#' (dispersal, residency) _a priori_, this status was checked with analytical tools in the 
#' following scripts.

# --------------- label=load_data2, eval=TRUE

# load meta data
meta_data <- readr::read_csv("data/meta_data.csv")

# Load data
mov_data_tiete <- readr::read_csv("data/raw_movement_data_Pardas_Tiete_all_individuals.csv") %>% 
  dplyr::select(tagID:tempC)

# mov_data_tiete
# ncol(mov_data_tiete)
str(mov_data_tiete)
unique(mov_data_tiete$id)
head(mov_data_tiete$timestamp)

# left join with metadata, filter by release date
mov_data_tiete <- mov_data_tiete %>%
  # sf::st_transform(crs = 4326) %>%
  # dplyr::bind_cols(dplyr::as_tibble(sf::st_coordinates(.))) %>%
  # sf::st_drop_geometry() %>%
  # dplyr::mutate(timestamp = lubridate::with_tz(timestamp, tzone = "UTC")) %>%
  dplyr::left_join(meta_data, by = c("id" = "name")) #%>%
  # dplyr::as_tibble() %>%
  # dplyr::filter(timestamp >= release_date_utc & timestamp <= last_date_utc)

# mov_data_tiete %>%
#   print(width = Inf)

# Load Puma and Jaguar data - Legado das Aguas
mov_data_legado <- readr::read_csv("data/raw_movement_data_Legado_Aguas_all_individuals.csv") %>% 
  dplyr::mutate(timestamp = as.POSIXct(timestamp, format = "%Y/%m/%d %H:%M:%S", tz = "UTC")) %>% 
  dplyr::rename(longitude = X, latitude = Y, tagID = Tag_ID, id_number = ID, id = name) %>% 
  dplyr::select(-name.old) %>% 
  dplyr::filter(species == "Puma concolor") %>% 
  dplyr::select(longitude:timestamp)
  
# mov_data_legado
ncol(mov_data_legado)
str(mov_data_legado)

# take only pumas, filter by release date
mov_data_legado <- mov_data_legado %>% 
  dplyr::left_join(meta_data, by = c("id" = "name"), suffix = c("", ".y")) %>% 
  dplyr::filter(timestamp >= release_date_utc, timestamp <= last_date_utc) %>% 
  dplyr::rename(shortname = name_old)

# --------------- label=combine_data

#' ## Combining data
#' 
#' After loading data, column names were standardized and the `data.frame`s were combined in a single
#' one. 

# --------------- label=combine_data2, eval=TRUE

# Combine information from all individuals into a single data.frame

# Standardize column names
colnames(mov_data_tiete)
colnames(mov_data_legado)

# Columns to be selected
cols_to_keep <- c("tagID", "id", "id_number", "shortname", "longitude", "latitude", "timestamp",
                  "hDOP", "species", "sex", "weight_kg", "estimated_age_months",
                  "study_name", "country", "dispersal_behavior")

# Keep only these columns and merge studies
mov_data <- mov_data_tiete %>% 
  dplyr::select(all_of(cols_to_keep)) %>% 
  dplyr::bind_rows(
    mov_data_legado %>% 
      dplyr::select(all_of(cols_to_keep))
  ) %>% 
  dplyr::arrange(species, id_number, timestamp) %>% 
  dplyr::rename(X = longitude, Y = latitude) %>% 
  dplyr::distinct(id, timestamp, .keep_all = TRUE)

# mov_data

# Check
mov_data$id %>% unique
mov_data$species %>% unique

mov_data %>% 
  dplyr::group_by(id) %>% 
  dplyr::summarise(ID = unique(id_number)) %>% 
  dplyr::arrange(ID) %>% 
  print(n = 25)

# summary
mov_data %>% 
  dplyr::group_by(id) %>% 
  dplyr::summarise(
    start = min(timestamp),
    end = max(timestamp),
    days = difftime(max(timestamp), min(timestamp), units = "days") %>% round(),
    n = n()
  ) #%>% 
  # write.csv("data/summary_data.csv", row.names = F)

#' ## Exporting organized data

#' We export the data in different formats to allow multiple visualizations and analyses.

# --------------- label=export_combined, eval=FALSE

# export as geopackage
mov_data_sf <- mov_data %>% 
  sf::st_as_sf(coords = c("X", "Y"), remove = FALSE, crs = 4326)

sf::st_write(mov_data_sf, dsn = "data/movement_data_pumas_tiete_legado.gpkg",
             delete_dsn = TRUE)

# export as shapefile
# split maps into individual tracks in each year and exporting them as shapefile
mov_data_split <- mov_data_sf %>%
  dplyr::select(id_number, id, X, Y, timestamp, geometry) %>%
  dplyr::mutate_at(c(1:5), as.character) %>%
  split(.$id)

# export shapefiles for each individual
dir_shp <- "data/03_movement_data_shp/"
if(!dir.exists(dir_shp)) 
  dir.create(dir_shp)

lapply(names(mov_data_split), function(x)
  sf::st_write(mov_data_split[[x]],
               dsn = paste0(dir_shp, x, ".shp"),
               delete_dsn = TRUE))

# save as csv file
mov_data %>%
  write.csv(file = "data/movement_data_pumas_tiete_legado.csv", row.names = F)

# export csv for each individual
dir_csv <- "data/04_csv_each_individual/"
if(!dir.exists(dir_csv)) 
  dir.create(dir_csv)

mov_data %>%
  split(.$id) %>%
  purrr::map(., function(x)
    x %>%
      write.csv(file = paste0(dir_csv, "movement_data_", x$id[1], ".csv"),
                row.names = F))

# save only data from Pardas do Tiete project, to upload on MoveBank
mov_data %>% 
  dplyr::filter(study_name == "Pardas_Tiete") %>% 
  write.csv(file = "data/movement_data_pumas_tiete_toMoveBank.csv", row.names = F)

# save as rda
save(mov_data, file = "data/movement_data_pumas_tiete_legado.rda")

#' # Data availabilty in MoveBank and within the `PardasIPC` R package
#' 
#' The organized GPS data is available on MoveBank, under the study name
#' "Projeto Pardas do Tiete/Pumas from Tiete Project", which can be assessed 
#' [here](https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study577226894).
#' Besides, the data are also organized in a R package called `PardasIPC`
#' (after Pardas (pumas) from the studies of the Instituto Pró-Carnívoros, IPC).
#' The R package is private, but if you are given authorization to access it, you
#' can install it and get the data using:
#' 

# --------------- label=install_package, eval=FALSE

library(devtools)
devtools::install_github("bniebuhr/puma_dispersal_homerange_AF", subdir = "PardasIPC")
