#'---
#' title: "Puma space use and dispersal in the tropics: Organize raw data"
#' author: Bernardo Niebuhr
#'---

#-----------------------------
# Set up 

# Remove old stuff
rm(list = ls())

# Load packages
if(!require(install.load)) install.packages('install.load'); library(install.load)
install.load::install_load('tidyverse', 'lubridate', 'sf')

# Directories

# Data folder
datadir <- "data/01_new_data_raw/"
# Output folder - the same as the data folder, since we're just re-organizing data
outdir <- datadir

#-----------------------------
# Load data

# Read metada and transform all non-numeric columns in strings
meta.data <- readr::read_csv('data/meta_data.csv')
meta.data
meta.data %>% print(width = Inf)

# Read movement data

# Female Cougars - Pardas do Tiete
# F1 - Sucuri and F2 - Poran
mov_data_puma_f1e2 <- readr::read_tsv('data/01_new_data_raw/pardas_tiete_f1e2_Sucuri_Pora_2018_02_d20.csv', locale = locale(decimal_mark = ',')) %>% 
  dplyr::mutate(UTC_Date = as.character(UTC_Date),
         UTC_Time = as.character(UTC_Time),
         timestamp = as.POSIXct(paste(UTC_Date, UTC_Time), format = '%d/%m/%Y %H:%M', tz = 'UTC')) %>% 
  dplyr::arrange(Name, timestamp)


mov_data_puma_f1e2
mov_data_puma_f1e2 %>% print(width = Inf)
ncol(mov_data_puma_f1e2)
str(mov_data_puma_f1e2)

# Male cougars - Nick and Aracatuba
mov_data_puma_m1 <- readr::read_csv('data/01_new_data_raw/pardas_tiete_m1_Nick_2018_05_d11.csv', locale = locale(decimal_mark = '.')) %>% 
  dplyr::mutate(UTC_Date = as.character(UTC_Date),
         UTC_Time = as.character(UTC_Time),
         timestamp = as.POSIXct(paste(UTC_Date, UTC_Time), format = '%Y-%m-%d %H:%M:%S', tz = 'UTC')) %>% 
  dplyr::arrange(timestamp)

mov_data_puma_m1 %>% print(width = Inf)

mov_data_puma_m2 <- readr::read_csv('data/01_new_data_raw/pardas_tiete_m2_Aracatuba_2015_09_d20.csv', locale = locale(decimal_mark = '.')) %>% 
  dplyr::mutate(UTC_Date = as.character(UTC_Date),
         UTC_Time = as.character(UTC_Time),
         timestamp = as.POSIXct(paste(UTC_Date, UTC_Time), format = '%Y-%m-%d %H:%M:%S', tz = 'UTC')) %>% 
  dplyr::arrange(Name, timestamp)

mov_data_puma_m2 %>% print(width = Inf)

# Other individuals
files1 <- list.files(datadir, pattern = "csv", full.names = T)[-c(1:3)]
files2 <- list.files(datadir, pattern = "TXT", full.names = T)

# Sirtrack collars
sirtrack <- list()

# for all, export as shp
for(i in files1) {
  name <- str_split(i, "_")[[1]][7]
  
  sirtrack[[name]] <- readr::read_csv(i) %>% 
    dplyr::filter(!is.na(Latitude)) %>% 
    dplyr::select(Tag_ID:MinVolt) %>% 
    dplyr::mutate(timestamp = lubridate::ymd_hms(paste(UTC_Date, UTC_Time)),
                  Name = name) %>% 
    dplyr::rename(hDOP = HDOP)
}

# Lotek collars
lotek <- list()

for(i in files2) {
  name <- str_split(i, "_")[[1]][7]
  
  lotek[[name]] <- read.csv(i, stringsAsFactors = F) %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(!is.na(Latitude)) %>% 
    dplyr::mutate(timestamp = lubridate::mdy_hms(Date...Time..GMT.)) %>% 
    dplyr::rename(Name = 1, Tag_ID = 2, hDOP = DOP, Temperature_C = Temp..C.) %>% 
    dplyr::select(Name:Main..V., timestamp)
}

#-----------------------------
# Organize data

# Select columns

# Columns to be kept
# (cols.to.keep <- ind.data[[1]] %>% 
#    select(c(1:6,11:ncol(.))) %>% 
#    colnames)
cols_to_drop_sirtrack <- c("UTC_Date", "UTC_Time", 'CNR', 'Sats', 'TimeOn(s)', 'MinVolt')
cols_to_drop_lotek <- c("Date...Time..GMT.", "Date...Time..Local.", "Fix.Status", "Main..V.")
cols_to_drop_other <- c("UTC_Date", "UTC_Time", 'TempDeg', 'Activity', 'CNR', 'Sats', 
                        'TimeOn(s)', 'MinVolt', "Sex", "Name_old", "ID")

# Gather datasets in a single data.frame and order in time
mov.datasets <- sirtrack %>% 
  dplyr::bind_rows() %>% 
  dplyr::select(-cols_to_drop_sirtrack) %>%
  dplyr::bind_rows(
    dplyr::bind_rows(mov_data_puma_f1e2, mov_data_puma_m1, mov_data_puma_m2) %>%
      dplyr::select(-cols_to_drop_other)
  ) %>%
  dplyr::bind_rows(
    dplyr::bind_rows(lotek) %>%
      dplyr::select(-cols_to_drop_lotek)
  ) %>%
  dplyr::arrange(Name, timestamp)

mov.datasets$Name[grepl("Tup", mov.datasets$Name)] <- "Tupa"

# Verify
unique(mov.datasets$Name)
mov.datasets

# Add columns with information

# Add columns of information, filter by release date, drop not useful columns
colnames(meta.data)
colnames(mov.datasets)

mov.datasets.all <- mov.datasets %>% 
  dplyr::left_join(meta.data, by = c("Name" = "name")) %>%
  dplyr::filter(timestamp >= release_date_utc, Latitude > -40, timestamp <= last_date_utc) %>%
  dplyr::select(-c(release_date_utc, last_date_utc, observation, obs2, 
                   ends_with('.y'))) %>% 
  dplyr::rename(Tag_ID = Tag_ID.x, name = Name, shortname = name_old)

mov.datasets.all
ncol(mov.datasets.all)
str(mov.datasets.all)

# Verify by plot

# all
ggplot(data = mov.datasets.all) +
  geom_path(aes(x = Longitude, y = Latitude, color = name)) +
  coord_quickmap()

ggplot(data = mov.datasets.all) +
  geom_path(aes(x = Longitude, y = Latitude)) +
  facet_wrap(~ name, scales = 'free')

# Problems with individual Rafiki
# Remove outliers
mov.datasets.all %>% 
  dplyr::filter(name == 'Rafiki') %>% 
  ggplot() +
  geom_path(aes(x = Longitude, y = Latitude))

outlier.lines <- mov.datasets.all %>%
  rownames_to_column() %>%
  filter(name == 'Rafiki', Latitude < -21.950) %>%
  dplyr::pull(rowname) %>%
  as.numeric

mov.datasets.all <- mov.datasets.all %>%
  filter(!(row_number() %in% outlier.lines))

# Check
ggplot(data = mov.datasets.all) +
  geom_path(aes(x = Longitude, y = Latitude)) +
  facet_wrap(~ name, scales = 'free')

#-------------
# Export data

mov.datasets.all <- mov.datasets.all %>% 
  dplyr::select(tagID = Tag_ID, id = name, id_number = ID, shortname, timestamp,
                longitude = Longitude, latitude = Latitude, altitude = Altitude, hDOP, 
                tempC = Temperature_C, species:country) %>% 
  dplyr::arrange(id_number, timestamp)
mov.datasets.all

# Separate studies from Tiete and Legado
studies <- unique(mov.datasets.all$study_name)

# Transform data into spatial
shp <- mov.datasets.all %>% 
  sf::st_as_sf(coords = c('longitude', 'latitude'), crs = 4326, remove = F) # CRS(+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs)

# Include coordinates in SIRGAS2000-UTM 22S
shp.complete <- cbind(
  shp,
  shp %>%
    sf::st_transform(crs = 31982) %>% # CRS(+proj=utm +zone=22 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs)
    sf::st_coordinates()) %>% 
  dplyr::rename(x_GRS80_utm22S = X, y_GRS80_utm22S = Y)

# Include coordinates in SIRGAS2000-Albers
shp.complete <- cbind(
  shp.complete,
  shp.complete %>%
    sf::st_transform(crs = "+proj=aea +lat_1=-2 +lat_2=-22 +lat_0=-12 +lon_0=-54 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs") %>% # CRS(+proj=utm +zone=22 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs)
    sf::st_coordinates()) %>% 
  dplyr::rename(x_GRS80_albers = X, y_GRS80_albers = Y)

# Export data
for(i in studies) {
  to.save <- shp.complete %>% 
    dplyr::filter(study_name == i)
  
  # csv
  to.save %>% 
    sf::st_drop_geometry() %>% 
    write.csv(file = paste0('data/02_new_data_sorted/movement_data_', i, '_all_individuals_',
                               str_replace_all(lubridate::today(), '-', '_'), '.csv'), 
              row.names = F)
  
  # csv for each individual
  to.save %>% 
    sf::st_drop_geometry() %>% 
    split(.$id) %>% 
    map(., function(x) 
      x %>% 
        write.csv(file = paste0('data/02_new_data_sorted/movement_data_', i, "_", x$id[1], '.csv'), 
                        row.names = F))
  
  # shp
  to.save %>% 
    dplyr::select(tagID:sex) %>% 
    sf::st_write(paste0('data/02_new_data_sorted/shape/movement_data_', i, '_all_individuals_',
                               str_replace_all(lubridate::today(), '-', '_'), '.shp'),
               delete_dsn = TRUE)
  # gpkg
  sf::st_write(to.save, paste0('data/02_new_data_sorted/shape/movement_data_', i, '_all_individuals_',
                               str_replace_all(lubridate::today(), '-', '_'), '.gpkg'),
               delete_dsn = TRUE)
  
  
  # shp for each individual
  to.save %>%
    dplyr::select(tagID:sex) %>% 
    dplyr::rename(X = longitude, Y = latitude) %>% 
    dplyr::mutate_if(is.double, as.character) %>% 
    split(.$id) %>% 
    map(., function(x) 
      x %>% 
        sf::st_write(paste0('data/02_new_data_sorted/shape/movement_data_', i, "_", x$id[1], '.shp'), 
                  delete_dsn = TRUE))
}

# csv for all data up to 2019
to.save %>% 
  sf::st_drop_geometry() %>% 
  dplyr::filter(lubridate::year(timestamp) <= 2019) %>% 
  write.csv(file = paste0('data/02_new_data_sorted/movement_data_', i, '_individuals_up_to_2019', '.csv'), 
            row.names = F)
