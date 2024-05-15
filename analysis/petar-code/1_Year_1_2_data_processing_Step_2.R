# Author: Petar Bursac
# Date: 13.11.2023.

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
# PART 2 - data processing 
# ------------------------------------------------------------------------------

# Workflow:
# 1. Subset dataset to animal of interest
# 2. Read the calibration error model
# 3. Remove first 2 and last 2 days from data deployment
# 4. Make telemetry object and assign calibration error model
# 5. Make 4 outlier detection plots 



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
head(data.y2)

# ------------------------------------------------------------------------------
# NOTE:
# Example is shown for one animal.
# Later it can be wrapped into the function to process one by one animal - id 
# ------------------------------------------------------------------------------


# Subset dataset to animal of interest
# ------------------------------------------------------------------------------

# 11  (f/SI)
data.11 <- data.y1 %>% dplyr::filter(`individual-local-identifier` == '6b')

# Then read the calibration model
# ------------------------------------------------------------------------------

calibration_model <- readRDS("Data/undep/calibration.error.model.inventa_no_class.rds")
calibration_model

# Remove first 2 and last 2 days from data deployment
# ------------------------------------------------------------------------------

range(data.11$timestamp)

data.11 %<>% 
  dplyr::arrange(Date) %>%
  dplyr::filter(!(Date %in% seq(min(Date), min(Date + 2), by = "day"))) %>%
  dplyr::filter(!(Date %in% seq(max(as.Date(Date) - 2), max(as.Date(Date)), by = "day")))

range(data.11$timestamp)


# Some visualization
# ------------------------------------------------------------------------------

# check if times are correct
mov <- data.11 %>%
  arrange(timestamp) %>%
  transmute(timestamp,
            time = hour(timestamp) + minute(timestamp) / 60,
            dt = as.numeric(timestamp - lag(timestamp)) / (1 %#% 'hours'),
            distance = sqrt((`location-long` - lag(`location-long`))^2 +
                              (`location-lat` - lag(`location-lat`))^2))

# histograms of times between locations
mov %>%
  ggplot() +
  geom_histogram(aes(as.numeric(dt)), bins = 20) +
  scale_x_log10(expression(Time~between~locations~(log[10]~scale))) +
  ylab('Counts') +
  theme_bw()


# histograms of animals' average times weighted by displacement
summarize(mov,
          x = mean(time,
                   weights = distance / mean(distance, na.rm = TRUE),
                   na.rm = TRUE)) %>%
  ggplot() +
  geom_histogram(aes(x))

# hex plot of straight-line displacement by time of day
mov %>%
  ggplot() +
  geom_hex(aes(time, distance), bins = 10) +
  labs(x = 'Time of day', y = 'SLD between consecutive locations') +
  theme_bw() +
  scale_fill_viridis_c('Count', limits = c(1, NA), trans = 'log10')



# Make telemetry object and assign calibration error model
# ------------------------------------------------------------------------------

data.tel.11 <- as.telemetry(data.11, datum = 'WGS84')

# Assign calibration error model
uere(data.tel.11) <- calibration_model
uere(data.tel.11)

plot(data.tel.11)


data.11 %>% dplyr::group_by(`gps:fix-type-raw`) %>%
  summarize(total = n())



data.11.sf <- st_as_sf(data.11, coords = c("location-long", "location-lat"))
mapviewOptions(fgb = FALSE)

mapview(data.11.sf)

# Make 4 outlier detection plots 
# ------------------------------------------------------------------------------

# Link: https://github.com/QuantitativeEcologyLab/bc-mammals-speeds/tree/main/analysis

theme_set(theme_bw())

# source custom functions
source('R_scripts/From_Stefano/S_Functions.R') # all custom functions
# NOTE: change the path to the file S_Functions.R

# Make 4 outlier plots
# -----------------------------------

# Solution 1
# -----------------------------------

png(filename = paste0('Data/MoveBank_both_years/figures/outlier-diagnostics/',
                      "11", '.png'),
    width = 10, height = 10, units = 'in', res = 300)

outlier_plots(telemetry = data.tel.11, return = FALSE)

dev.off()

# Solution 2 - step by step 
# -----------------------------------

# Add extra columns to be consistent with code
d <- data.11 %>% dplyr::mutate(animal = `individual-local-identifier`, 
                               outlier = FALSE, 
                               species = `individual-taxon-canonical-name`) %>%
  dplyr::rename(location.long = `location-long`,
                location.lat = `location-lat`, 
                dataset_name = `study-name`)


d %<>%
  filter(! is.na(timestamp)) %>%
  mutate(hdop = NA, 
         pdop = NA,
         DOP = NA,
         dop = `gps:dop`,
         original_outliers = outlier,
         outlier = if_else(outlier, 1, 0)) %>% # make numeric
  # select(! c(HDOP, DOP, gps.dop)) %>% # remove duplicates
  relocate(dop, .after = pdop) %>%
  nest(tel = ! c(species, dataset_name, animal))


# check which species have DOP values
d %>%
  unnest(tel) %>%
  group_by(species) %>%
  summarise(across(c(dop, hdop, pdop), ~ sum(.x, na.rm = TRUE))) %>%
  mutate(has_dop = dop + hdop + pdop > 0) %>%
  arrange(desc(has_dop))

tel.d <- as.telemetry(d$tel[[1]], mark.rm = TRUE) # create telemetry object
uere(tel.d) <- calibration_model # assign calibration error model

tel.d %>%
  outlier_plots()

# There is no difference ...


# Check animal
# -----------------------------------

# With CEM - Calibration Error Model
out <- check_animal('11', calibration = TRUE, calibration_model = calibration_model, return = TRUE)
out[out$speed > 0.4, ]
plot_adj('11', max_speed = 0.4, calibration = TRUE, calibration_model = calibration_model)

# Without CEM
out <- check_animal('11', calibration = FALSE, calibration_model = calibration_model, return = TRUE)
out[out$speed > 0.4, ]
plot_adj('11', max_speed = 0.4, calibration = FALSE, calibration_model = calibration_model)

# There is no difference ...


plot_adj('11', max_speed = 0.4, map = TRUE)
plot_adj('11', max_angle = 170, max_speed = 0.4, calibration = TRUE, calibration_model = calibration_model) 

out[out$speed > 0.2, ]

plot_adj('11', max_angle = 90, max_speed = 0.2, calibration = TRUE, calibration_model = calibration_model) 
plot_adj('11', max_angle = 90, max_speed = 0.2, map = TRUE, calibration = TRUE, calibration_model = calibration_model) 

out[out$angle > 170 & out$speed > 0.5, ] 



# Process all animals and export outlier plots
# ------------------------------------------------------------------------------

ids.y1 <- unique(data.y1$`individual-local-identifier`)

for(i in 1:length(ids.y1)) {
  
  data.id <- data.y1 %>% dplyr::filter(`individual-local-identifier` == ids.y1[i])
  
  data.id %<>% 
    dplyr::arrange(Date) %>%
    dplyr::filter(!(Date %in% seq(min(Date), min(Date + 1), by = "day"))) %>%
    dplyr::filter(!(Date %in% seq(max(as.Date(Date) - 1), max(as.Date(Date)), by = "day")))
  
  data.tel.id <- as.telemetry(data.id, datum = 'WGS84')
  
  png(filename = paste0('Data/MoveBank_both_years/figures/outlier-diagnostics_new/Year 1/',
                        ids.y1[i], '.png'),
      width = 10, height = 10, units = 'in', res = 300)
  
  outlier_plots(telemetry = data.tel.id, return = FALSE)
  
  dev.off()
  
  print(paste("Exported: ", i, " / ", length(ids.y1)))
  
}




ids.y2 <- unique(data.y2$`individual-local-identifier`)

for(i in 1:length(ids.y2)) {
  
  data.id <- data.y2 %>% dplyr::filter(`individual-local-identifier` == ids.y2[i])
  
  data.id %<>% 
    dplyr::arrange(Date) %>%
    dplyr::filter(!(Date %in% seq(min(Date), min(Date + 1), by = "day"))) %>%
    dplyr::filter(!(Date %in% seq(max(as.Date(Date) - 1), max(as.Date(Date)), by = "day")))
  
  data.tel.id <- as.telemetry(data.id, datum = 'WGS84')
  
  png(filename = paste0('Data/MoveBank_both_years/figures/outlier-diagnostics_new/Year 2/',
                        ids.y2[i], '.png'),
      width = 10, height = 10, units = 'in', res = 300)
  
  outlier_plots(telemetry = data.tel.id, return = FALSE)
  
  dev.off()
  
  print(paste("Exported: ", i, " / ", length(ids.y2)))
  
}



