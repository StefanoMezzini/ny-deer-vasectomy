# function that loops through a moving window of tracking data
# computes movement parameters and CIs for each window
# plots ML point estimate of HR size and 95% CIs
# written by Stefano Mezzini based on code by Mark Bidwell and Chris Fleming

library('ctmm')  # for continuous movement modeling
library('dplyr') # for data wrangling
library('purrr') # for functional programming
library('ggplot2') # for fancy plots
theme_set(theme_bw() + theme(legend.position = 'none'))

window_hr <- function(tel, window, dt, projection, full_ud = NULL,
                      fig_path = NULL, rds_path = NULL, cores = 1,
                      akde_weights = FALSE, plot_50_q = FALSE) {
  
  cat('Assuming the telemetry error is calibrated.\n')
  
  if(.Platform$OS.type == 'Windows' & cores > 1) {
    warning('cores > 1 only works on Unix machines.')
  }
  
  # extract 95% AKDE for the full telemetry, if available
  if(! is.null(full_ud) & class(full_ud) == 'UD') {
    HR_0 <- summary(full_ud, units = FALSE)$CI['area (square meters)', 'est'] / 1e6
  } else {
    HR_0 <- NA_real_
  }
  
  #' moving window of size `window` slides by `dt`
  #' `---`.... --> .`---`... --> ..`---`.. --> ...`---`. --> ....`---`
  times <- seq(min(tel$t), max(tel$t) - window, by = dt)
  N <- length(times)
  
  # to extract home ranges later
  extract_hr <- function(a, par, l.ud) {
    # convert HR to km^2
    summary(a, units = FALSE, level.UD = l.ud)$CI['area (square meters)', ][par] / 1e6
  }
  
  out <-
    tibble(
      # add start and end times
      t_start = times, # left bound
      t_end = t_start + window, # right bound
      #' subset times within window (can't use `filter()` on telemetry obj)
      dataset = map2(t_start, t_end,
                     function(t_1, t_2) tel[tel$t >= t_1 & tel$t <= t_2, ]),
      models = imap(
        dataset,
        \(d, i) {
          
          if(nrow(d) > 1) {
            cat('Analyzing dataset ', i, ' of ', N, '.\n', sep = '')
            tibble(
              # find initial guesses for models (assuming calibrated error)
              guess = ctmm.guess(data = d, interactive = FALSE,
                                 CTMM = ctmm(error = TRUE)) %>%
                list(),
              # select best movement model based on subset of tel
              model = ctmm.select(data = d, CTMM = guess[[1]],
                                  cores = cores) %>%
                list(),
              # estimate autocorrelated kernel density estimate
              akde = akde(data = d, CTMM = model[[1]],
                          weights = akde_weights) %>%
                list(),
              # find home range estimate
              hr_est_95 = extract_hr(a = akde[[1]], par='est', l.ud=0.95),
              hr_lwr_95 = extract_hr(a = akde[[1]], par='low', l.ud=0.95),
              hr_upr_95 = extract_hr(a = akde[[1]], par='high', l.ud=0.95))
          } else {
            tibble(
              guess = list('Insufficient data.'),
              model = list('Insufficient data.'),
              akde = list('Insufficient data.'),
              hr_est_50 = NA_real_,
              hr_lwr_50 = NA_real_,
              hr_upr_50 = NA_real_,
              hr_est_95 = NA_real_,
              hr_lwr_95 = NA_real_,
              hr_upr_95 = NA_real_)
          } # close else
        })) %>% # close function for imap()
    tidyr::unnest(models) %>%
    mutate(diffusion = map_dbl(model, \(.m) {
      summary(mw$model[[1]], units = FALSE)$
        CI['diffusion (square meters/second)', 'est']
    }),
    speed = map_dbl(model, \(.m) {
      summary(mw$model[[1]], units = FALSE)$
        CI['speed (meters/second)', 'est']
    }),
    #' `ctmm::%#%` syntax: `"new units" %#% value in SI units`
    units_area = 'km^2', # already converted
    units_diff = 'km^2/day',
    diffusion = units_diff %#% diffusion,
    units_speed = 'km/day',
    speed = units_speed %#% speed,
    dof_area = map_dbl(model, \(.m) summary(.m)$DOF['area']),
    dof_diff = map_dbl(model, \(.m) summary(.m)$DOF['diffusion']),
    dof_speed = map_dbl(model, \(.m) summary(.m)$DOF['speed'])) %>%
    mutate(t_center = (t_start + t_end) / 2,
           posixct = as.POSIXct(t_center, origin = '1970-01-01',
                                tz = tel@info$timezone),
           date = as.Date(posixct))
  
  if(! is.null(rds_path)) {
    saveRDS(out,
            file.path(rds_path,
                      paste0(tel@info['identity'],
                             '-window-', window / (1 %#% 'day'), '-days',
                             '-dt-', dt / (1 %#% 'day'), '-days.rds')))
  }
  
  # plot results ----
  # tracking data
  plt_a <-
    ggplot(tel) +
    coord_equal() +
    geom_point(aes(longitude, latitude, color = timestamp)) +
    geom_path(aes(longitude, latitude), alpha = 0.1) +
    scale_color_viridis_c('Time') +
    labs(x = '', y = NULL)
  
  plt_b <-
    ggplot(out) +
    
    # 95% CIs for 95% home range estimates
    geom_ribbon(aes(date, ymin = hr_lwr_95, ymax = hr_upr_95), alpha = 0.3) +
    
    # 95% home range estimates
    geom_line(aes(date, hr_est_95), linewidth = 1.25) +
    geom_line(aes(date, hr_est_95, color = posixct)) +
    
    # home range from model fit to full dataset
    geom_hline(yintercept = HR_0, color = 'darkorange', na.rm = TRUE) +
    
    scale_x_date(NULL, date_labels = '%b %Y') +
    scale_color_viridis_c() +
    labs(y = expression(Home~range~(km^2)))
  
  if(plot_50_q) {
    plt_b <- plt_b +
      # core home range 95% CIs
      geom_ribbon(aes(date, ymin = hr_lwr_50, ymax = hr_upr_50), alpha = 0.3) +
      # core home range estimates
      geom_line(aes(date, hr_est_50), linewidth = 1.25) +
      geom_line(aes(date, hr_est_50, color = posixct))
  }
  
  plt <- cowplot::plot_grid(plt_a, plt_b, labels = c('A', 'B'), nrow = 1,
                            align = 'hv')
  
  if (! is.null(fig_path)) {
    # Save figure as a png using the animal's name
    file.path(fig_path, paste0(tel@info['identity'],
                               '-window-', window / (1 %#% 'day'),
                               '-days-dt-', dt / (1 %#% 'day'),
                               '-days.png')) %>%
      ggsave(plot = plt, units = 'in', width = 7, height = 3,
             dpi = 600, bg = 'white')
  } else {
    print(plt)
  }
}

if(FALSE) {
  # Test the function on the buffalo data
  data('buffalo')
  
  # Apply to a single individual
  test <- window_hr(tel = buffalo$Queen, # buffalo with smallest dataset
                    window = 7 %#% 'day', # size of the moving window
                    dt = 1 %#% 'day', # step size for the moving window
                    fig_path = NULL) # can specify where to save the figure
  
  # apply to all buffaloes in the dataset
  results <- lapply(buffalo,
                    window_hr,
                    # arguments to pass to window_hr
                    window = 30 %#% 'day',
                    dt = 2 %#% 'day')
}
