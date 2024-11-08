library('dplyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for working with data
library('ctmm')      # for movement modeling
library('ggplot2')   # for fancy plots
source('analysis/figures/default-theme.R')

# find names of moving window files
files <- list.files(path = 'data',
                    pattern = '-window-7-days-dt-3-days.rds',
                    full.names = TRUE, recursive = TRUE)
length(files) # should be 124

# import moving window data
d <- map_dfr(
  files,
  function(.filename) {
    year <- substr(.filename,
                   nchar('data/Year-x'),
                   nchar('data/Year-x'))
    
    .d <- readRDS(.filename) %>%
      select(! distance_est:distance_upr) %>%
      select(! ind.id) %>% # ID is incorrect
      mutate(study_year = as.numeric(year),
             animal =
               substr(.filename, nchar('data/Year X/X'), nchar(files)) %>%
               gsub(pattern = '-.*', perl = TRUE,
                    replacement = '')) %>%
      relocate(study_year, .before = 1) %>%
      relocate(animal, .before = 1)
    
    # some datasets extra unnecessary columns
    if('d.length' %in% colnames(.d)) {
      .d <- select(.d, ! d.length)
    }
    if('d.length.int' %in% colnames(.d)) {
      .d <- select(.d, ! d.length.int)
    }
    
    return(.d)
  })

colnames(d)
n_distinct(d$animal) # should be 114
n_distinct(paste(d$animal, d$study_year)) # should be 124
unique(d$animal)

# drop windows with insufficient data
filter(d, map_chr(model, \(x) class(x)) != 'ctmm')$model
d <- filter(d, map_chr(model, \(x) class(x)) == 'ctmm')

# add speed and diffusion estimates, and DOF for each parameter ----
d <-
  mutate(d,
         # estimates
         diff_lwr = map_dbl(model, \(.m) ctmm:::diffusion(.m)[1]),
         diff_est = map_dbl(model, \(.m) ctmm:::diffusion(.m)[2]),
         diff_upr = map_dbl(model, \(.m) ctmm:::diffusion(.m)[3]),
         speed_obj = map(model, \(.m) speed(.m, units = FALSE) %>%
                           suppressWarnings()), # for non-OUF models
         speed_est = map_dbl(speed_obj, \(.s) .s$CI[, 'est']),
         #' change missing speeds from `Inf` to `NA`
         speed_est = if_else(is.finite(speed_est), speed_est, NA_real_),
         # add units and convert accordingly
         area_units = 'km^2', # already converted
         diff_units = 'km^2/day',
         diffusion_est = diff_units %#% diffusion_est,
         speed_units = 'km/day',
         speed_est = speed_units %#% speed_est,
         # find degrees of freedom (returns 0 if value is NA)
         dof_area = map_dbl(model, \(.m) ctmm:::DOF.area(.m)),
         dof_diff = map_dbl(model, \(.m) summary(.m)$DOF['diffusion']),
         dof_speed = map_dbl(model, \(.m) ctmm:::DOF.speed(.m))) %>%
  # only keep necessary columns
  rename(tel = dataset) %>%
  select(study_year, animal, posixct, date,
         tel, model, akde,
         dof_area, hr_units, hr_est_95, hr_lwr_95, hr_upr_95,
         dof_speed, speed_units, speed_est,
         dof_diff, diffusion_units, diffusion_est) %>%
  # add metadata on each individual
  left_join(read.csv('data/reference_data.csv'),
            by = c('animal' = 'id', 'study_year')) %>%
  # add factors for models
  mutate(animal = factor(animal),
         study_year = factor(study_year),
         animal_year = factor(paste(animal, study_year)),
         sex_treatment = factor(paste(sex, study_site)),
         s_t_y = factor(paste(sex, study_site, study_year))) %>%
  # drop first ~10 days for each animal to remove odd behaviors
  group_by(animal, study_year) %>%
  filter(date >= min(date) + 10) %>% # equivalent to using the first day
  ungroup() %>%
  # create a column of days since August 1st
  mutate(days_since_aug_1 =
           if_else(
             # if in Aug, Sept, Oct, Nov, or Dec
             month(date) >= 8,
             # then calculate days since Aug 1
             date - as.Date(paste0(year(date), '-08-01')),
             # otherwise calculate days since Aug 1 of previous year 
             date - as.Date(paste0(year(date) - 1, '-08-01'))) %>%
           as.numeric(), # convert difftime to numeric
         has_fawn =
           (study_year == 1 & animal %in% c('12', '14b', '15b', '16')) |
           (study_year == 2 & animal %in% c('1', '16', '123')))

# count number of deer per site and sex
d %>%
  group_by(sex_treatment, study_year) %>%
  summarise(n_tel = n_distinct(animal_year))

# check if fawns have been assigned correctly
filter(d, has_fawn) %>%
  group_by(study_year) %>%
  summarize(n = n_distinct(animal),
            which = paste(unique(animal), collapse = ', '))

# check which treatments the does with fawns are in
d %>%
  filter(sex == 'f') %>%
  group_by(study_site, study_year) %>%
  summarize(prop = mean(has_fawn))

# save the final dataset ----
# version without akdes to push to keep < 100 MB and push to GitHub
d %>%
  select(! akde) %>%
  saveRDS(file = 'data/years-1-and-2-data-no-akde.rds')

# version with AKDEs
saveRDS(object = d, file = 'data/years-1-and-2-data-akde.rds')

# create a csv of which animal is included in which analysis
d_summary <-
  d %>%
  group_by(animal, study_year) %>%
  summarise(
    #' check if at least 1 value is not `NA`
    has_hr = any(! is.na(hr_est_95)),
    has_speed = any(! (is.na(speed_est) | is.infinite(speed_est))),
    has_diffusion = any(! is.na(diffusion_est)),
    #' check proportion of values that are not `NA`
    prop_finite_hr = mean(! is.na(hr_est_95)) %>%
      round(2),
    prop_finite_speed = mean(! (is.na(speed_est) | is.infinite(speed_est))) %>%
      round(2),
    prop_finite_diffusion = mean(! is.na(diffusion_est)) %>%
      round(2),
    .groups = 'drop')

write.csv(d_summary, 'data/ids-for-each-analysis.csv', row.names = FALSE)
