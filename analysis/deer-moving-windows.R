library('tidyr') #' for `nest()`
library('furrr') # for fast parallel computations
source('functions/window_ctmm.R')

# model for the estimated root-mean-square User Equivalent Range Error
deer_uere <- readRDS('models/calibration.error.model.inventa_no_class.rds')

# import both years of data
d <- bind_rows(read.csv('data/Year_1_processed_data.csv'),
               read.csv('data/Year_2_processed_data.csv')) %>%
  mutate(animal = individual.local.identifier) %>%
  group_by(study.year, animal) %>% # animals were collared each year
  nest(tel = - c(study.year, animal)) %>%
  mutate(tel = map(tel, \(.tel) {
    .tel <- as.telemetry(.tel)
    uere(.tel) <- deer_uere # add collar error calibration
    return(.tel)
  }))

# ensure years 1 and 2 are labelled separately
unique(d$study.year)

# set number of cores to run on -- will be parallelizing by collar
#' **do not set CORES at the number of cores you have available or higher**
CORES <- availableCores(logical = FALSE) - 2

plan(multisession, workers = CORES)

future_map(d$tel, function(.d) {
  .d <- window_ctmm(.d,
                    window = 14 %#% 'day',
                    dt = 14 %#% 'day',
                    fig_path = 'figures/moving-window',
                    rds_path = 'models/moving-window',
                    cores = 1) # cannot parallelize on Windows
},
.progress = TRUE,
.options = furrr_options(seed = TRUE))
