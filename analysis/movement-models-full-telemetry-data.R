library('ctmm')  # for movement models
library('dplyr') # for data wrangling
library('tidyr') # for data wrangling
library('purrr') # for functional programming
library('furrr') # for parallel computing

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
          list())
    },
    .progress = TRUE)) %>%
  unnest(models)

# extract movement parameters ----
#' # find home range estimate
#' summary(full_ud, units = FALSE)$CI['area (square meters)', 'est']
#' hr_est_95 = extract_hr(a = akde[[1]], par='est', l.ud=0.95),
#' hr_lwr_95 = extract_hr(a = akde[[1]], par='low', l.ud=0.95),
#' hr_upr_95 = extract_hr(a = akde[[1]], par='high', l.ud=0.95),
#' # find diffusion estimates
#' diffusion = map_dbl(model, \(.m) {
#'   if(any(grepl('diffusion',
#'                rownames(summary(.m)$CI)))) {
#'     return(summary(.m, units = FALSE)$
#'              CI['diffusion (square meters/second)', 'est'])
#'   } else {
#'     return(NA_real_)
#'   }}),
#' # find speed
#' speed = map_dbl(model, \(.m) {
#'   if(any(grepl('speed',
#'                rownames(summary(.m)$CI)))) {
#'     return(summary(.m, units = FALSE)$
#'              CI['speed (meters/second)', 'est'])
#'   } else {
#'     return(NA_real_)
#'   }}),
#' # find degrees of freedom
#' dof_area = summary(model[[1]])$DOF['area'],
#' dof_diff = if_else(
#'   condition = grepl('OU', summary(model[[1]])$name),
#'   true = summary(model[[1]])$DOF['speed'],
#'   false = NA_real_),
#' dof_speed = if_else(
#'   condition = grepl('OUF', toupper(summary(model[[1]])$name)),
#'   true = summary(model[[1]])$DOF['speed'],
#'   false = NA_real_))
#' })) %>% # close function for imap()
#'   tidyr::unnest(models) %>%
#'   #' `ctmm::%#%` syntax: `"new units" %#% value in SI units`
#'   mutate(units_area = 'km^2',
#'          area = units_diff %#% area,
#'          units_diff = 'km^2/day',
#'          diffusion = units_diff %#% diffusion,
#'          units_speed = 'km/day',
#'          speed = units_speed %#% speed,
#'          t_center = (t_start + t_end) / 2,
#'          posixct = as.POSIXct(t_center, origin = '1970-01-01',
#'                               tz = .tel@info$timezone),
#'          date = as.Date(posixct))
