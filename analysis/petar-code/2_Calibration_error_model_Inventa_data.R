# Author: Petar Bursac (bursac.petar3@gmail.com)
# Date: 14.10.2023.

# Load major packages
library(dplyr)
library(magrittr)
library(ggplot2)
library(sf)
library(data.table)
library(tidyr)
library(mapview)
library(readxl)
library(lubridate)
library(ctmm)

mapviewOptions(fgb = FALSE)

# Workflow from paper "A comprehensive framework for handling location error in animal tracking data"
# Fleming et. al, 2021 - Figure 1


# Calibration data import
# ------------------------------------------------------------------------------

calibration.data <- readr::read_csv(".../gps_inventa.csv") %>%
  as.data.frame()

# Check range "2021-08-01 21:01:09 UTC" TO "2021-08-09 03:02:47 UTC"
range(calibration.data$timestamp) 

# Create summary table with n= number of events per GPS fix type. 
calibration.data %>% 
  dplyr::group_by(`Fix-type`) %>%
  dplyr::summarise(n_per_type = n()) %>% 
  as.data.frame()


# Visualization
# ------------------------------------------------------------------------------

# Make sf - spatial class object from data
cal.data.sf <- st_as_sf(calibration.data, # specify dataset
                        coords = c("location-long", # specify columns with coordinates
                                   "location-lat"), # specify coordinate reference system
                        crs = 4326)

# view data interactive on web map
mapview(cal.data.sf, # specify data - sf spatial object
        zcol = "individual.id", # specify column with which you want to color the points by value
        layer.name = "Tag ID") # optional: specify layer name which will be displayed in map legend


# Error calibration - model
# ------------------------------------------------------------------------------

# This approach uses data pulled directly from Inventa with 75 items labelled "No Fix" removed.

head(calibration.data)
calibration.data %<>% dplyr::select(-`Fix-type`)


# ------------------------------------------------------------------------------

# Make telemetry object from data
calibration.tel <- as.telemetry(calibration.data, 
                                datum = 'WGS84') 

# Warnings, using DOP, minimum sampling interval warning occurs due to multiple serial numbers in one file.
plot(calibration.tel)
names(calibration.tel)

# Plot the data
plot(calibration.tel, col=rainbow(length(unique(calibration.data$`individual.id`))))

# Use the uere command to estimate the RMS UERE parameter(s) from calibration data.
# Fit the model
UERE <- uere.fit(calibration.tel)
#Warning occurs due to low number of fix types in certain categories
summary(UERE)
UERE$DOF
# 2D and 3D are low sample size and low degrees of freedom
# The UERE parameters are assigned to a dataset with the uere()<- command.

UERES <- lapply(calibration.tel, uere.fit)
#warnings small sample size for individual
summary(list(joint = UERE, individual = UERES)) #  individualized and joint error models
# results in joint error model preferred for horizontal

# ------------------------------------------------------------------------------

# Choose joint error model.
UERE

# Save the calibration error model
# ------------------------------------------------------------------------------
saveRDS(UERE, ".../calibration.error.model.inventa_no_class.rds")




