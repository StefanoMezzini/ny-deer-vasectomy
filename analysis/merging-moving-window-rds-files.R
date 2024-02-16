library('dplyr') # for data wrangling
library('purrr') # for map_dfr()
library('ctmm')  # for movement modeling

files <- list.files(path = 'data',
                    pattern = '-window-7-days-dt-3-days.rds',
                    full.names = TRUE, recursive = TRUE)
length(files) # should be 124

d <- map_dfr(files,
             function(.filename) {
               year <- substr(.filename,
                              nchar('data/Year-x'),
                              nchar('data/Year-x'))
               
               .d <- readRDS(.filename) %>%
                 select(! distance_est:distance_upr) %>%
                 mutate(study_year = as.numeric(year)) %>%
                 relocate(study_year, .before = 1)
               # some datasets have a column of window size called d.length
               if('d.length' %in% colnames(.d)) {
                 .d <- select(.d, ! d.length)
               }
               if('d.length.int' %in% colnames(.d)) {
                 .d <- select(.d, ! d.length.int)
               }
               return(.d)
             }) %>%
  rename(animal = ind.id)

colnames(d)

filter(d, is.na(animal))$model # two windows with insufficient data
d <- filter(d, ! is.na(animal)) # drops rows with no movement model

# all units are consistent
unique(d$diffusion_units)
unique(d$tau_p_units)
unique(d$tau_v_units)
unique(d$hr_units)
unique(d$speed_model_units) # from the speed() function
unique(d$speed_gauss_units) # from the movement model

# only keep necessary columns
d <- d %>%
  select(
    study_year,
    animal,
    tau_p_units, tau_p_est, tau_p_lwr, tau_p_upr,
    tau_v_units, tau_v_est, tau_v_lwr, tau_v_upr,
    speed_gauss_units, speed_gauss_est, speed_gauss_lwr, speed_gauss_upr,
    speed_model_units, speed_model_est, speed_model_lwr, speed_model_upr,
    hr_units, hr_est_50, hr_lwr_50, hr_upr_50,
    hr_est_95, hr_lwr_95, hr_upr_95,
    diffusion_units, diffusion_est, diffusion_lwr, diffusion_upr,
    posixct, date)

# add metatada on each individual
metadata <- read.csv('data/reference_data.csv') %>%
  rename(animal = id)

d_final <- left_join(d, metadata, by = c('animal', 'study_year'))

# save the final dataset
saveRDS(object = d_final, file = 'data/years-1-and-2-data.rds')
