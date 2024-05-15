
# Load major packages

library(move)
library(ctmm)

library(dplyr)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(sf)
library(data.table)
library(tidyr)
library(mapview)
library(readxl)
library(lubridate)
library(stringr)
library(ggspatial)
library(gridExtra)

mapviewOptions(fgb = FALSE)

# ------------------------------------------------------------------------------
# PART 3 - flag outliers
# ------------------------------------------------------------------------------

# Workflow:
# 1. Subset dataset to animal of interest
# 2. Read the calibration error model
# 3. Remove first 2 and last 2 days from data deployment
# 4. Make telemetry object and assign calibration error model
# 5. Make 4 outlier detection plots 
# 6. Flag outliers
# 7. Export the outliers per individual (if some events are flagged as outliers)



# Read the datasets
# ------------------------------------------------------------------------------

data.y1 <- data.table::fread("Data/MoveBank_both_years/Year_1_processed_data.csv",
                             stringsAsFactors = FALSE,
                             header = TRUE) %>%
  as.data.frame()


data.y2 <- data.table::fread("Data/MoveBank_both_years/Year_2_processed_data.csv",
                             stringsAsFactors = FALSE,
                             header = TRUE) %>%
  as.data.frame()


head(data.y1)


# Subset dataset to animal of interest
# ------------------------------------------------------------------------------

unique(data.y1$`individual-local-identifier`)
unique(data.y2$`individual-local-identifier`)


# Specify the ID 
individual.id <- '6b'


# Choose the individual - pay attention to data.y1 or data.y2
data.subset <- data.y1 %>% 
  dplyr::filter(`individual-local-identifier` == individual.id)

# Then read the calibration model
# ------------------------------------------------------------------------------

calibration_model <- readRDS("Data/undep/calibration.error.model.inventa_no_class.rds")
calibration_model

rownames(calibration_model$UERE) <- 'gps'
rownames(calibration_model$DOF) <- 'gps'


# Remove first 2 and last 2 days from data deployment
# ------------------------------------------------------------------------------

range(data.subset$timestamp)

data.subset %<>% 
  dplyr::arrange(Date) %>%
  dplyr::filter(!(Date %in% seq(min(Date), min(Date + 1), by = "day"))) %>%
  dplyr::filter(!(Date %in% seq(max(as.Date(Date) - 1), max(as.Date(Date)), by = "day")))

range(data.subset$timestamp)


# Make telemetry object and assign calibration error model
# ------------------------------------------------------------------------------

data.subset %<>% dplyr::select(-`gps:fix-type-raw`)

data.tel <- as.telemetry(data.subset, datum = 'WGS84')

# Assign calibration error model
uere(data.tel) <- calibration_model
uere(data.tel)

# Visualization
plot(data.tel)

data.subset.sf <- st_as_sf(data.subset, coords = c("location-long", "location-lat"))
mapviewOptions(fgb = FALSE)

mapview(data.subset.sf)

# ------------------------------------------------------------------------------


# Make 4 outlier detection plots 
# ------------------------------------------------------------------------------

theme_set(theme_bw())

# Source custom functions
source('R_scripts/From_Stefano/S_Functions.R') # all custom functions
# NOTE: change the path to the file S_Functions.R


outlier_plots(telemetry = data.tel, return = TRUE)



# Prepare the data for flagging
# ------------------------------------------------------------------------------

# Add extra columns to be consistent with code
d <- data.subset %>% dplyr::mutate(animal = `individual-local-identifier`, 
                               outlier = FALSE, 
                               species = `individual-taxon-canonical-name`) %>%
  dplyr::rename(location.long = `location-long`,
                location.lat = `location-lat`, 
                dataset_name = `study-name`)


d %<>%
  filter(!is.na(timestamp)) %>%
  mutate(hdop = NA, 
         pdop = NA,
         DOP = NA,
         dop = `gps:dop`,
         original_outliers = outlier,
         outlier = if_else(outlier, 1, 0)) %>% # make numeric
  relocate(dop, .after = pdop) %>%
  nest(tel = ! c(species, dataset_name, animal))



# Check animal
# -----------------------------------

out <- check_animal(individual.id, calibration = TRUE, calibration_model = calibration_model, return = TRUE)

head(out)

out[out$speed > 0.4, ] # filter by speed value
out[out$angle > 170, ] # filter by angle
out[out$distance > 1500, ] # core deviation = distance
out[out$dt > 14400, ] # filter by time (in seconds 3600s = 1hour)


plot_adj(individual.id, max_speed = 0.4, calibration = TRUE, calibration_model = calibration_model)

# Create map 
plot_adj(individual.id, max_angle = 170, max_speed = 0.4, calibration = TRUE, calibration_model = calibration_model, 
         map = TRUE) 

# Add angle boundary (if needed)
out[out$angle > 170 & out$speed > 0.4, ] 
out[out$angle > 170 & out$speed > 0.4 & out$dt > 14400, ] 
out[out$angle > 170 & out$speed > 0.4 & out$distance > 1500, ]

# & and 
# | or

# Then flag that outlier (if desicion is TRUE and it is outlier)
# just based on max_speed, max_angle, max_distance
flag_outlier(id = individual.id, max_speed = 0.4, max_angle = 170, value = 1) 


# Check if there is added flags 
d$tel[[1]] %>% dplyr::filter(outlier == 1) %>% 
  as.data.frame()


# Option 1
# Save the event-ids for flags
# -----------------------------------
out.flags <- d$tel[[1]] %>% dplyr::filter(outlier == 1) %>% as.data.frame()

# Option 2
# Flag based on event-id
# -----------------------------------
events_to_flag <- c(23489271248)
out.flags <- d$tel[[1]] %>% dplyr::filter(as.character(`event-id`) %in% as.character(events_to_flag)) %>% dplyr::mutate(outlier = 1) %>% as.data.frame()


# Change the path 
write.csv(out.flags, paste0("Data/MoveBank_both_years/outliers/Year 1/", individual.id, ".csv"))




