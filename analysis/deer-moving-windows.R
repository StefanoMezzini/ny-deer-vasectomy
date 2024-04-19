###########################################################################
#' **not run; the moving widows were fit by Petar.**
#' the code should produce something equivalent
###########################################################################
library('tidyr') #' for `nest()`
library('furrr') # for fast parallel computations
source('functions/window_ctmm.R')

# model for the estimated root-mean-square User Equivalent Range Error
deer_uere <- readRDS('models/calibration.error.model.inventa_no_class.rds')

# import both years of data
d <- bind_rows(read.csv('data/Year_1_processed_data.csv'),
               read.csv('data/Year_2_processed_data.csv')) %>%
  mutate(animal = individual.local.identifier,
         collar = paste0(animal, '_y', study.year)) %>%
  group_by(study.year, animal) %>% # animals were collared each year
  nest(tel = - c(study.year, animal, collar)) %>%
  mutate(tel = map(tel, \(.tel) {
    .tel <- as.telemetry(.tel)
    uere(.tel) <- deer_uere # add collar error calibration
    return(.tel)
  })) %>%
  ungroup()

# ensure years 1 and 2 are labelled separately
unique(d$study.year)

# for a quick test
if(FALSE) {
  test <- window_ctmm(d$tel[[1]][c(1:6, 120:145), ],
                      window = 1 %#% 'day',
                      dt = 1 %#% 'day',
                      fig_path = NULL,
                      rds_path = NULL,
                      cores = 1, # cannot parallelize on Windows
                      progress = 1) # can't have progress if parallelized
}

# set number of cores to run on -- will be parallelizing by collar
#' **do not set CORES at or above the number of cores you have available**
CORES <- availableCores(logical = FALSE) - 1
plan(multisession, workers = CORES)

future_map2(d$tel, d$study.year, function(.telemetry, .sy) {
  window_ctmm(.tel = .telemetry,
              study_year = .sy,
              window = 14 %#% 'day',
              dt = 14 %#% 'day',
              fig_path = paste0('figures/moving-window/', d$study.year),
              rds_path = paste0('models/moving-window/', d$study.year),
              cores = 1, # cannot parallelize on Windows
              progress = 0) # can't have progress if parallelized
},
.progress = TRUE,
.options = furrr_options(seed = TRUE))
