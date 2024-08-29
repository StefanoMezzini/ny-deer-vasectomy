# Author: Petar Bursac (bursac.petar3@gmail.com)
# Date: 11.11.2023.

# Load major packages

library(move)
library(ctmm)

library(tidyverse)
library(dplyr)
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

# ------------------------------------------------------------------------------
# PART 1 - data analysis and pre-processing
# ------------------------------------------------------------------------------

# Workflow:
# 1. Read the dataset - from MoveBank
# 2. Remove events without coordinates (zero) which is in the same time "NO-FIX"
# 3. Add date column 
# 4. Read table with additional information per individual, filter individual ids and Join - reference data
# 5. Remove by transmission protocol
# 6. Final check for duplicated events - timestamps
# 7. Check for devices that only have 1 fix/every 4 hours
# 8. Export the final data per Year (1 and 2)


# Read the dataset - from MoveBank
# ------------------------------------------------------------------------------

# Dataset for both years
data.raw <- data.table::fread("Data/MoveBank_both_years/Odocoileus virginianus DeNicola Staten Island, NY and Rockefeller Park, NY.csv",
                              stringsAsFactors = FALSE,
                              header = TRUE) %>%
  as.data.frame()

# Check the names
names(data.raw)

# Check the length of unique ids and tags
length(unique(data.raw$`tag-local-identifier`)) 
length(unique(data.raw$`individual-local-identifier`)) 


# Remove events without coordinates
# ------------------------------------------------------------------------------

# First check number of events per type and number of locations without coordinates
data.raw %>%
  dplyr::group_by(`gps:fix-type-raw`) %>%
  dplyr::summarise(total = n(),
                   at_least_one_na_coordinate = sum(is.na(`location-lat`) | is.na(`location-lat`)))


data.raw %<>% dplyr::filter(`gps:fix-type-raw` != "NO FIX") # remove events without fixed type of GPS observation
# and also in this way, zero coordinates are removed

# Check again
data.raw %>%
  dplyr::group_by(`gps:fix-type-raw`) %>%
  dplyr::summarise(total = n(),
                   at_least_one_na_coordinate = sum(is.na(`location-lat`) | is.na(`location-lat`)))


# Add date column 
# ------------------------------------------------------------------------------
data.raw %<>% 
  dplyr::mutate(Date = as.Date(timestamp), # add Date column 
                Time = format(as.POSIXct(timestamp), format = "%H:%M:%S")) %>% # add Time column
  dplyr::mutate(Year = substr(Date, 1, 4), # add Year, Month and Day columns separately from Date column with substring function
                Month = substr(Date, 6, 7),
                Day = substr(Date, 9, 10))

range(data.raw$timestamp) # check

# check start and end date per individual
data.raw %>%
  dplyr::group_by(`individual-local-identifier`) %>%
  dplyr::summarise(start_date = min(timestamp), 
                   end_date = max(timestamp)) %>%
  as.data.frame()


# Read table with additional information per individual
# ------------------------------------------------------------------------------

reference_data <- readxl::read_xlsx(path = "Data/MoveBank_both_years/reference_data.xlsx", 
                                    sheet = "Data") %>%
  as.data.frame()

# Remove rows with tag - do not consider
reference_data %<>% dplyr::filter(is.na(`X = do not consider`))

head(reference_data)

# Check refrence data - ids in both years
# ------------------------------------------------------------------------------
reference_data %>% dplyr::filter(duplicated(`animal-id`)) %>% dplyr::select(`animal-id`)
ids.dup <- c("51", "1123", "25", "4b", "16", "27", "29", "30", "43", "1")
reference_data %>% dplyr::filter(`animal-id` %in% ids.dup)

# 10 ids available in both years

# Split to Year 1 and 2 and filter by unique id and date on-off
# ------------------------------------------------------------------------------

reference_data %<>% dplyr::rename(`animal-age` = Age)

reference_data_y1 <- reference_data %>% dplyr::filter(`study-year` == 1)
reference_data_y2 <- reference_data %>% dplyr::filter(`study-year` == 2)

# Year 1
# -------------------------------------
list.y1 <- list()
for(i in 1:dim(reference_data_y1)[1]) {
  
  ref.data <- reference_data_y1[i, ]
  
  data.id <- data.raw %>% dplyr::filter(`individual-local-identifier` == ref.data$`animal-id`)
  data.id %<>% filter(between(Date, as.Date(ref.data$`deploy-on-date`), as.Date(ref.data$`deploy-off-date`)))
  
  # join attributes (study area, sex, year)
  data.id %<>% left_join(., ref.data %>% dplyr::select(`animal-id`, `study-site`, `animal-sex`, `study-year`, `animal-age`), by = c("individual-local-identifier" = "animal-id"))
  data.id %<>% dplyr::rename(study.area = `study-site`, study.year = `study-year`)
  
  list.y1[[i]] <- data.id
  
  print(paste0("Added: ", i, " / ", dim(reference_data_y1)[1], " --- Year 1"))
  
}

data.y1 <- data.table::rbindlist(list.y1) %>%
  as.data.frame()


# Year 2
# -------------------------------------
list.y2 <- list()
for(j in 1:dim(reference_data_y2)[1]) {
  
  ref.data <- reference_data_y2[j, ]
  
  data.id <- data.raw %>% dplyr::filter(`individual-local-identifier` == ref.data$`animal-id`)
  data.id %<>% filter(between(Date, as.Date(ref.data$`deploy-on-date`), as.Date(ref.data$`deploy-off-date`)))
  
  # join attributes (study area, sex, year)
  data.id %<>% left_join(., ref.data %>% dplyr::select(`animal-id`, `study-site`, `animal-sex`, `study-year`, `animal-age`), by = c("individual-local-identifier" = "animal-id"))
  data.id %<>% dplyr::rename(study.area = `study-site`, study.year = `study-year`)
  
  list.y2[[j]] <- data.id
  
  print(paste0("Added: ", j, " / ", dim(reference_data_y2)[1], " --- Year 2"))
  
}

data.y2 <- data.table::rbindlist(list.y2) %>%
  as.data.frame()



# Analyse by fix type - keep all classes
# ------------------------------------------------------------------------------

data.y1 %>%
  dplyr::group_by(`gps:fix-type-raw`) %>%
  dplyr::summarise(total = n(),
                   min_dop = min(`gps:dop`, na.rm = TRUE),
                   max_dop = max(`gps:dop`, na.rm = TRUE),
                   min_sat = min(`gps:satellite-count`, na.rm = TRUE),
                   max_sat = max(`gps:satellite-count`, na.rm = TRUE))


data.y2 %>%
  dplyr::group_by(`gps:fix-type-raw`) %>%
  dplyr::summarise(total = n(),
                   min_dop = min(`gps:dop`, na.rm = TRUE),
                   max_dop = max(`gps:dop`, na.rm = TRUE),
                   min_sat = min(`gps:satellite-count`, na.rm = TRUE),
                   max_sat = max(`gps:satellite-count`, na.rm = TRUE))


# Analyse by altitude - height (above or below 500 meters) - gross errors
# ------------------------------------------------------------------------------

data.y1 %>%
  group_by(`gps:fix-type-raw`) %>%
  filter(abs(`height-above-ellipsoid`) > 500) %>%
  summarise(n = n(),
            median_long = median(`location-long`),
            median_lat = median(`location-lat`))

data.y2 %>%
  group_by(`gps:fix-type-raw`) %>%
  filter(abs(`height-above-ellipsoid`) > 500) %>%
  summarise(n = n(),
            median_long = median(`location-long`),
            median_lat = median(`location-lat`))


# Only analyse by gps:sat and gps:dop values
# ------------------------------------------------------------------------------
range(data.y1$`gps:satellite-count`, na.rm = T)
range(data.y1$`gps:dop`)

data.y1 %>% dplyr::group_by(`gps:satellite-count`) %>% dplyr::summarize(n_count = n())
data.y1 %>% dplyr::group_by(`gps:dop`) %>% dplyr::summarize(n_count = n())

data.y1 %>%
  group_by(`gps:satellite-count`) %>%
  summarize(n_count = n()) %>%
  mutate(percent = round(n_count / sum(n_count) * 100, 4)) %>%
  arrange(desc(percent))


range(data.y2$`gps:satellite-count`, na.rm = T)
range(data.y2$`gps:dop`)

data.y2 %>% dplyr::group_by(`gps:satellite-count`) %>% dplyr::summarize(n_count = n())
data.y2 %>% dplyr::group_by(`gps:dop`) %>% dplyr::summarize(n_count = n())

data.y2 %>%
  group_by(`gps:satellite-count`) %>%
  summarize(n_count = n()) %>%
  mutate(percent = round(n_count / sum(n_count) * 100, 4)) %>%
  arrange(desc(percent))


# Format timestamp attribute 
# ------------------------------------------------------------------------------
data.y1 %<>%
  dplyr::mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS"))

data.y2 %<>%
  dplyr::mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%OS"))


# Analyse and remove events which is duplicated
# eg. remove by transmission protocol
# ------------------------------------------------------------------------------

# Number of NA satellite count per fix type
data.y1 %>%
  dplyr::filter(is.na(`gps:satellite-count`)) %>%
  group_by(`gps:fix-type-raw`) %>%
  dplyr::summarise(total = n())

data.y2 %>%
  dplyr::filter(is.na(`gps:satellite-count`)) %>%
  group_by(`gps:fix-type-raw`) %>%
  dplyr::summarise(total = n())

# Remove by transmission protocol
data.y1 %>%  
  dplyr::group_by(`transmission-protocol`, `gps:fix-type-raw`) %>%
  dplyr::summarise(total = n())

data.y2 %>%  
  dplyr::group_by(`transmission-protocol`, `gps:fix-type-raw`) %>%
  dplyr::summarise(total = n()) 

data.y1 %>%
  dplyr::filter(`transmission-protocol` == "")

data.y2 %>%
  dplyr::filter(`transmission-protocol` == "") # nothing to remove


# Remove only blanks in Year 1
data.y1 %<>%
  dplyr::filter(`transmission-protocol` != "")


# Final check for duplicated events - timestamps
# ------------------------------------------------------------------------------

data.y1 %>%
  group_by(`individual-local-identifier`) %>%
  summarize(dupl_count = sum(duplicated(timestamp))) %>%
  as.data.frame()

data.y2 %>%
  group_by(`individual-local-identifier`) %>%
  summarize(dupl_count = sum(duplicated(timestamp))) %>%
  as.data.frame()

# Average time per individual and 
# Check for devices that only have 1 fix/every 4 hours
# ------------------------------------------------------------------------------

# Year 1
average_time1 <-
  data.y1 %>% 
  group_by(`individual-local-identifier`) %>% 
  mutate(timestamp_Next = lead(timestamp)) %>% 
  # ungroup() %>% 
  mutate(
    diff_days = difftime(timestamp_Next, timestamp, units = 'days'),
    diff_hours = difftime(timestamp_Next, timestamp, units = 'hours'),
    diff_mins = difftime(timestamp_Next, timestamp, units = 'mins'),
    diff_secs = difftime(timestamp_Next, timestamp, units = 'secs')
  ) %>% 
  summarise(avg_hour = mean(diff_hours, na.rm = TRUE), study.year = first(study.year)) %>%
  as.data.frame()

average_time1

average_time1 %>% 
  ggplot() +
  geom_density(aes(x = avg_hour), lwd = 1)

gy1 <- ggplot(data = average_time1) +
  geom_point(data = average_time1, aes(y = avg_hour, x = as.factor(`individual-local-identifier`))) +
  geom_text(data = average_time1 %>% filter(avg_hour > 2), 
            aes(y = avg_hour, x = as.factor(`individual-local-identifier`), label = `individual-local-identifier`), 
            color = "red", 
            nudge_y = 0.2, 
            nudge_x = 0.2) +
  #facet_wrap(~study.year) +
  labs(x = "Individual ID", title = "Average time - Year 1") +
  #coord_flip()
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 0.5))
gy1

# Year 2 
average_time2 <-
  data.y2 %>% 
  group_by(`individual-local-identifier`) %>% 
  mutate(timestamp_Next = lead(timestamp)) %>% 
  # ungroup() %>% 
  mutate(
    diff_days = difftime(timestamp_Next, timestamp, units = 'days'),
    diff_hours = difftime(timestamp_Next, timestamp, units = 'hours'),
    diff_mins = difftime(timestamp_Next, timestamp, units = 'mins'),
    diff_secs = difftime(timestamp_Next, timestamp, units = 'secs')
  ) %>% 
  summarise(avg_hour = mean(diff_hours, na.rm = TRUE), study.year = first(study.year)) %>%
  as.data.frame()

average_time2

average_time2 %>% 
  ggplot() +
  geom_density(aes(x = avg_hour), lwd = 1)

gy2 <- ggplot(data = average_time2) +
  geom_point(data = average_time2, aes(y = avg_hour, x = as.factor(`individual-local-identifier`))) +
  geom_text(data = average_time2 %>% filter(avg_hour > 2), 
            aes(y = avg_hour, x = as.factor(`individual-local-identifier`), label = `individual-local-identifier`), 
            color = "red", 
            nudge_y = 0.2, 
            nudge_x = 0.2) +
  #facet_wrap(~study.year) +
  labs(x = "Individual ID", title = "Average time - Year 2") +
  #coord_flip()
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 0.5))


gr1 <- grid.arrange(gy1, gy2, ncol = 1)

# export the plot
ggsave(plot = gr1,
       filename = "Data/MoveBank_both_years/average_time_Year_1_and_2.jpg",
       width = 32,
       height = 24,
       units = "cm",
       device = "jpeg",
       dpi = 700)



# Check for those far away points
# ------------------------------------------------------------------------------

data.y1.sf <- st_as_sf(data.y1, coords = c("location-long", "location-lat"), crs = 4326)
data.y2.sf <- st_as_sf(data.y2, coords = c("location-long", "location-lat"), crs = 4326)

sf::st_write(data.y1.sf, "Data/MoveBank_both_years/sf_year1.gpkg")
sf::st_write(data.y2.sf, "Data/MoveBank_both_years/sf_year2.gpkg")

# Disjoint - selected in QGIS (first drawn boundaries)
dj.year1 <- sf::st_read("Data/MoveBank_both_years/sf_year1_disjoint.gpkg")
dj.year2 <- sf::st_read("Data/MoveBank_both_years/sf_year2_disjoint.gpkg")

data.y1 %<>% dplyr::filter(!(as.numeric(`event-id`) %in% as.numeric(dj.year1$event.id)))
data.y2 %<>% dplyr::filter(!(as.numeric(`event-id`) %in% as.numeric(dj.year2$event.id)))


# Export the results
# ------------------------------------------------------------------------------
write.csv(data.y1, file = "Data/MoveBank_both_years/Year_1_processed_data.csv")
write.csv(data.y2, file = "Data/MoveBank_both_years/Year_2_processed_data.csv")


