library('dplyr')     # for data wrangling (mutate(), %>%, etc.)
library('tidyr')     # for data wrangling (nest(), unnest(), pivot_*, etc.)
library('purrr')     # for functional programming (map_***(), etc.)
library('ctmm')      # for movement models
library('lubridate') # for working with dates
library('mapview')   # for interactive maps
library('ggplot2')   # for fancy plots
theme_set(theme_bw())

# source custom functions
setwd('~/GitHub/bc-mammals-speeds')
source('functions/outlier_plots.R') # to plot outlier diagnostic plots
source('functions/check_animal.R') # to run diagnostic plots
source('functions/plot_adj.R') # to plot 20 adjacent locations
source('functions/flag_outlier.R') # to mark outliers
source('functions/remove_outlier_flags.R') # to start over with an animal

# change wd back to current folder
setwd('~/Work/vickie-de-nicola')

# import the full dataset ----
d <- bind_rows(readr::read_csv('data/Year_1_processed_data.csv',
                               col_types = '?') %>%
                 mutate(dataset = 1),
               readr::read_csv('data/Year_2_processed_data.csv',
                               col_types = '?') %>%
                 mutate(dataset = 2)) %>%
  select(! c(...1, `import-marked-outlier`)) %>% # drop rowname column
  rename_with(\(name) {
    gsub(':', '_', name) %>%
      gsub(pattern = '-', replacement = '_', x = .) %>%
      tolower()
  },
  everything()) %>%
  mutate(animal = individual_local_identifier,
         outlier = 0,
         original_outliers = FALSE) %>%
  rename(location.long = location_long, location.lat = location_lat) %>%
  nest(tel = ! c(individual_taxon_canonical_name, animal, study_name,
                 study.area, animal_sex, study.year, animal_age, dataset))

# check duplicated animal IDs
d$animal[which(duplicated(d$animal))]
length(d$animal[which(duplicated(d$animal))])

ID <- '1782'

out <- check_animal(ID)
plot_adj(ID, max_speed = 0.3, max_angle = 170)
plot_adj(ID, max_speed = 0.3, max_angle = 170, map = TRUE)
flag_outlier(ID, max_speed = 0.3, max_angle = 170, value = 1)
filter(filter(d, animal == ID)$tel[[1]], outlier > 0)

out <- check_animal(ID)
plot_adj(ID, max_speed = 0.5, max_angle = 170)
flag_outlier(ID, max_speed = 0.5, max_angle = 170, value = 1)

plot_adj(id = ID, max_angle = 160, max_speed = 0.175, n_adj = 5)
plot_adj(id = ID, max_angle = 160, max_speed = 0.175, map = TRUE)
plot_adj(id = ID, max_angle = 160, max_speed = 0.3, map = TRUE, n_adj = 50)
plot_adj(id = ID, max_dt = 4 %#% 'hours')

plot_adj(id = ID, max_angle = 160, max_speed = 0.5, min_angle = 175)
plot_adj(id = ID, max_angle = 160, max_speed = 0.5, min_angle = 175,
         map = TRUE)
plot_adj(id = ID, max_angle = 175, max_speed = 0.5, map = TRUE)
plot_adj(id = ID, max_angle = 160, max_speed = 0.3, map = TRUE)
plot_adj(id = ID, max_angle = 160, max_speed = 0.3,
         map = TRUE)

plot_adj(id = ID, map = TRUE, many = TRUE)
