library('ctmm')  # for movement models
library('dplyr') # for data wrangling
library('tidyr') # for data wrangling
library('purrr') # for functional programming
library('furrr') # for parallel computing

extract_hr <- function(.ud, par, l.ud) {
  summary(.ud, units = FALSE, level = 0.95,
          level.UD = l.ud)$CI['area (square meters)', par]
}

# set up parallel computation plan (separate R sessions)
CORES <- availableCores(logical = FALSE) - 1
plan(multisession, workers = CORES)
plan() # check plan

# fit movement models ----
# calibrated tels took too long (one was not done fitting after a month)
d <-
  bind_rows(read.csv('data/Year_1_processed_data.csv'),
            read.csv('data/Year_2_processed_data.csv')) %>%
  as_tibble() %>%
  mutate(animal = individual.local.identifier) %>% # duplicate column
  select(! X) %>% # drop row number
  nest(tel = ! c(animal, study.year, study.area, mortality.status,
                 animal.sex, animal.age)) %>% # create col of tel
  mutate(
    tel = map(tel, \(.tel) {
      # outliers already removed manually earlier
      as.telemetry(.tel, mark.rm = TRUE)
    }),
    models = future_map(tel, \(.tel) {
      tibble(
        # find initial guesses for models (assuming calibrated error)
        guess = ctmm.guess(data = .tel, interactive = FALSE) %>%
          list(),
        # select best movement model based on subset of .tel
        model = ctmm.select(data = .tel, CTMM = guess[[1]]) %>%
          list(),
        # estimate autocorrelated kernel density estimate
        akde = akde(data = .tel, CTMM = model[[1]],
                    weights = TRUE) %>% # weights using sampling frequency
          list(),
        # calculate estimated average density at each location 
        tel = map2(tel, akde, \(.tel, .akde) {
          .r <- raster(.akde, DF = 'PDF')
          .tel$density <- raster::extract(.r, SpatialPoints.telemetry(.tel))
        }),
        # find home range estimate
        hr_est_95 = extract_hr(.ud = akde[[1]], par='est', l.ud=0.95),
        # # find degrees of freedom
        dof_area = ctmm:::DOF.area(model[[1]]),
        dof_diff = summary(model[[1]])$DOF['diffusion'],
        dof_speed = ctmm:::DOF.speed(model[[1]]))
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))) %>%
  tidyr::unnest(models) %>%
  mutate(mean_diffusion = map_dbl(model, \(.m) { # find diffusion estimates
    if(any(grepl('diffusion',
                 rownames(summary(.m)$CI)))) {
      return(summary(.m, units = FALSE)$
               CI['diffusion (square meters/second)', 'est'])
    } else {
      return(NA_real_)
    }}),
    # find speed
    mean_speed = map_dbl(model, \(.m) {
      if(any(grepl('speed',
                   rownames(summary(.m)$CI)))) {
        return(summary(.m, units = FALSE)$
                 CI['speed (meters/second)', 'est'])
      } else {
        return(NA_real_)
      } #' close `else`
    } #' close `function`
    ),
    #' #' `ctmm::%#%` syntax: `"new units" %#% value in SI units`
    units_area = 'km^2',
    hr_est_95 = units_area %#% hr_est_95,
    units_diff = 'km^2/day',
    mean_diffusion = units_diff %#% mean_diffusion,
    units_speed = 'km/day',
    mean_speed = units_speed %#% mean_speed)

# save movement models ----
saveRDS(d, 'models/full-telemetry-movement-models.rds')
