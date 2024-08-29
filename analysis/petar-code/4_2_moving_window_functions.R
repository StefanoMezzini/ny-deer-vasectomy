# ------------------------------------------------------------------------------
# Some functions are adjusted - 
# based on work from Stefano Mezzini: 
# https://github.com/QuantitativeEcologyLab/hr-resource-stoch/blob/main/functions/window_hr.R
# ------------------------------------------------------------------------------


# Function: moving_window
# Author: Petar Bursac (bursac.petar3@gmail.com)
# Date: 05.01.2024.
# Description: Function to call window_hr_do_parallel and go trough ids and years

# Parameters: 
# individual.id - individual id of interest
# year.id - study year of interest
# data.y1 - dataset of Year 1
# data.y2 - dataset of Year 2
# calibration_model - calibration error model
# outliers.files.path - path to the files which contain outliers - flags
# window - size of the moving window
# dt - step size for the moving window - slide
# fig_path - path to save the figures
# cores - number of the cores used for the parallelization of processes
# rds_path - path to save the rds files


moving_window <- function(individual.id, year.id, data.y1 = data.y1, data.y2 = data.y2, calibration_model, outliers.files.path, window, dt, fig_path = NA, cores, rds_path = NA) {
  
  # Choose the individual - pay attention to data.y1 or data.y2
  
  print(paste0("Subset the data based on ind id: ", individual.id))
  st_time <- Sys.time()
  print(paste0("Start time: ", st_time))
  
  if(year.id == 1) {
    data.subset <- data.y1 %>% 
      dplyr::filter(`individual-local-identifier` == individual.id)
  } else if (year.id == 2) {
    data.subset <- data.y2 %>% 
      dplyr::filter(`individual-local-identifier` == individual.id)
  } else {
    print("Not available year...")
    return()
  }
  
  # Remove first 2 and last 2 days from data deployment
  # ----------------------------------------------------------------------------
  print("Removing first 2 and last 2 days from data deployment.")
  data.subset %<>% 
    dplyr::arrange(Date) %>%
    dplyr::filter(!(Date %in% seq(min(Date), min(Date + 1), by = "day"))) %>%
    dplyr::filter(!(Date %in% seq(max(as.Date(Date) - 1), max(as.Date(Date)), by = "day")))
  
  # Remove the outliers flagged manually
  # ----------------------------------------------------------------------------
  print("Removing the outliers flagged manually.")
  
  out1 <- list.files(paste0(outliers.files.path, "Year 1/")) %>%
    as.data.frame()
  
  out1 %<>% dplyr::rename(name = 1) %>%
    dplyr::mutate(id_name = str_remove(name, "v1.csv"),
                  full = list.files(paste0(outliers.files.path, "Year 1/"), full.names = TRUE))
  
  out1 %<>% dplyr::mutate(Year = 1)
  
  
  out2 <- list.files(paste0(outliers.files.path, "Year 2/")) %>%
    as.data.frame()
  
  out2 %<>% dplyr::rename(name = 1) %>%
    dplyr::mutate(id_name = str_remove(name, "v1.csv"),
                  full = list.files(paste0(outliers.files.path, "Year 2/"), full.names = TRUE))
  
  out2 %<>% dplyr::mutate(Year = 2)
  
  outliers <- rbind(out1, out2)
  
  if(individual.id %in% outliers$id_name) {
    
    outliers.id <- outliers %>% dplyr::filter(id_name == individual.id & Year == year.id)
    
    if(dim(outliers.id)[1] == 0) {
      print(paste0("Individual: ", individual.id, " does not have the maunnaly flagged outliers!"))
    } else {
      rm_event <- read.csv(outliers.id$full)
      
      data.subset %<>%
        dplyr::filter(!(`event-id`) %in% rm_event$`event-id`)
    }
    
  } else {
    
    print(paste0("Individual: ", individual.id, " does not have the maunnaly flagged outliers!"))
    
  }
  
  # Make telemetry object and assign calibration error model
  # ----------------------------------------------------------------------------
  print("Make telemetry object and assign calibration error model.")
  
  data.subset %<>% dplyr::select(-`gps:fix-type-raw`)
  data.tel <- as.telemetry(data.subset, datum = 'WGS84')
  # Assign calibration error model
  uere(data.tel) <- calibration_model

  # Run the moving window
  # ----------------------------------------------------------------------------
  print("Running moving window.")
  
  window_hr_do_parallel(tel = data.tel,
                        window = window, # size of the moving window
                        dt = dt, # step size for the moving window
                        fig_path = fig_path, # can specify where to save the figure 
                        cores = cores, 
                        rds_path = rds_path) 
                     
            
  print(paste0("Done processing for individual: ", individual.id, " from year: ", year.id))
  
  end_time <- Sys.time()
  print(paste0("End time: ", end_time))
  print(end_time - st_time)
  
  print("---------------------------------------------------------------------------------")
}


# ------------------------------------------------------------------------------


# Function: window_hr_do_parallel
# Author: Petar Bursac (bursac.petar3@gmail.com)
# Date: 05.01.2024.
# Description: Function to run the moving window analysis in parallel

# Parameters: 
# tel - individual id of interest
# window - size of the moving window
# dt - step size for the moving window - slide
# full_ud - if full ud is needed
# fig_path - path to save the figures
# rds_path - path to save the rds files
# cores - number of the cores used for the parallelization of processes
# weights - TRUE/FALSE - if weights are needed to estimate autocorrelated kernel density
# plot_50_q - TRUE/FALSE - if 50q plot is needed


window_hr_do_parallel <- function(tel, window, dt, full_ud = NULL,
                                  fig_path = NULL, rds_path = NULL, cores = 1,
                                  weights = FALSE, plot_50_q = FALSE) {
                               
  
  print(paste0("running in parallel on # of cores: ", cores))
  
  if(! is.null(full_ud) & class(full_ud) == 'UD') {
    HR_0 <- summary(full_ud, units = FALSE)$CI['area (square meters)',
                                               'est'] / 1e6
  } else {
    HR_0 <- NA_real_
  }
  
  # moving window created at beginning of t; slides forward as far as possible:
  # |---|.... --> .|---|... --> ..|---|.. --> ...|---|. --> ....|---|
  
  times <- seq(min(tel$t), max(tel$t) - window, by = dt)
  N <- length(times)
  
  # to extract home ranges later
  extract_hr <- function(a, par, l.ud) {
    # convert HR to km^2
    if("area (square meters)" %in% rownames(summary(a, units = FALSE, level.UD = l.ud)$CI)) {
      summary(a, units = FALSE, level.UD = l.ud)$CI['area (square meters)', ][par] / 1e6
    } else if ("area (hectares)" %in% rownames(summary(a, units = FALSE, level.UD = l.ud)$CI)) {
      summary(a, units = FALSE, level.UD = l.ud)$CI['area (hectares)', ][par] / 100
    } else if("area (square kilometers)" %in% rownames(summary(a, units = FALSE, level.UD = l.ud)$CI)){
      summary(a, units = FALSE, level.UD = l.ud)$CI['area (square kilometers)', ][par]
    } else {
      NA_real_
    }
  }
  
  # to extract speeds - Gaussian later
  extract_speed_gauss <- function(a, par) {
    # convert speed to kilometers/day
    if("speed (meters/second)" %in% rownames(a$CI)) {
      a$CI['speed (meters/second)', ][par] * (3600/1000*24)
    } else if ("speed (kilometers/day)" %in% rownames(a$CI)){
      a$CI['speed (kilometers/day)', ][par]
    } else{
      NA_real_
    }
  }
  
  # to extract diffusions later
  extract_diffusion <- function(a, par) {
    
    if("CI" %in% names(summary(a))){
      a = summary(a)
    } else {
      a = summary(a[[1]])
    }
    
    # convert diffusion to square kilometers/day
    if("diffusion (hectares/day)" %in% rownames(a$CI)) {
      a$CI['diffusion (hectares/day)', ][par] / 100
    } else if("diffusion (square meters/day)" %in% rownames(a$CI)) {
      a$CI['diffusion (square meters/day)', ][par] / 1000000
    } else if("diffusion (square kilometers/day)" %in% rownames(a$CI)) {
      a$CI['diffusion (square kilometers/day)', ][par]
    } else {
      NA_real_
    }
  }
  
  # to extract speeds - model later
  extract_speed_model <- function(a, par) {
    
    if("CI" %in% names(summary(a))){
      a = summary(a)
    } else {
      a = summary(a[[1]])
    }
    
    # convert speed to kilometers/day
    if("speed (meters/second)" %in% rownames(a$CI)) {
      a$CI['speed (meters/second)', ][par] * (3600/1000*24)
    } else if("speed (kilometers/day)" %in% rownames(a$CI)){
      a$CI['speed (kilometers/day)', ][par]
    } else {
      NA_real_
    }
  }
  
  # to extract tau p later
  extract_tau_p <- function(a, par) {
    
    if("CI" %in% names(summary(a))){
      a = summary(a, units = FALSE)
    } else {
      a = summary(a[[1]], units = FALSE)
    }
    
    
    if(any(grepl("[position]", rownames(a$CI), fixed = TRUE))) {
      a$CI[grepl("[position]", rownames(a$CI), fixed = TRUE), ][par]
    } else {
      NA_real_
    }
  }
  
  # to extract tau v later
  extract_tau_v <- function(a, par) {
    
    if("CI" %in% names(summary(a))){
      a = summary(a, units = FALSE)
    } else {
      a = summary(a[[1]], units = FALSE)
    }
    
    
    if(any(grepl("[velocity]", rownames(a$CI), fixed = TRUE))) {
      a$CI[grepl("[velocity]", rownames(a$CI), fixed = TRUE), ][par]
    } else {
      NA_real_
    }
  }
  
  
  times.df <- data.frame(t_start = times) %>%
    dplyr::mutate(t_end = t_start + window)
  
  
  registerDoParallel(cores=cores)
  
  # go trough windows
  
  models.list <- foreach(j = 1:dim(times.df)[1], .packages = c("dplyr", "magrittr", "ctmm")) %dopar% { 
                           
                           t_1 <- times.df[j, 1]
                           t_2 <- times.df[j, 2]
                           
                           d <- tel[tel$t >= t_1 & tel$t <= t_2, ]
                           
                           if(nrow(d) > 1) {
                             res.tb <- tibble(
                               
                               ind.id = tel@info['identity'][[1]],
                               
                               # create variogram
                               variogram = ctmm::variogram(d, error = TRUE) %>%
                                 list(),
                               
                               # find initial guesses for models
                               guess = ctmm.guess(data = d, variogram = variogram[[1]], interactive = FALSE) %>%
                                 list(),
                               
                               # select best model based on subset of tel
                               model = ctmm.select(data = d, CTMM = guess[[1]]) %>%
                                 list(),
                               
                               # estimate autocorrelated kernel density estimate
                               akde = ctmm::akde(data = d, CTMM = model[[1]], weights = weights) %>% 
                                 list(),
                               
                               # estimate Gaussian speed 
                               speed = ctmm::speed(object = d, CTMM = model[[1]], units = FALSE) %>%
                                 list(),
                               
                               # extract diffusion estimates
                               diffusion_units = "square kilometers/day",
                               
                               diffusion_est = extract_diffusion(a = model[[1]], par = 'est'),
                               diffusion_lwr = extract_diffusion(a = model[[1]], par = 'low'),
                               diffusion_upr = extract_diffusion(a = model[[1]], par = 'high'),
                               
                               # extract tau p estimates
                               tau_p_units = "seconds",
                               
                               tau_p_est = extract_tau_p(a = model[[1]], par = 'est'),
                               tau_p_lwr = extract_tau_p(a = model[[1]], par = 'low'),
                               tau_p_upr = extract_tau_p(a = model[[1]], par = 'high'),
                               
                               # extract tau p estimates
                               tau_v_units = "seconds",
                               
                               tau_v_est = extract_tau_v(a = model[[1]], par = 'est'),
                               tau_v_lwr = extract_tau_v(a = model[[1]], par = 'low'),
                               tau_v_upr = extract_tau_v(a = model[[1]], par = 'high'),
                               
                               # extract speed estimates - Gaussian
                               speed_gauss_units = "kilometers/day",
                               speed_gauss_est = extract_speed_gauss(a = speed[[1]], par = 'est'),
                               speed_gauss_lwr = extract_speed_gauss(a = speed[[1]], par = 'low'),
                               speed_gauss_upr = extract_speed_gauss(a = speed[[1]], par = 'high'),
                               
                               # extract speed from model
                               speed_model_units = "kilometers/day",
                               
                               speed_model_est = extract_speed_model(a = model[[1]], par = 'est'),
                               speed_model_lwr = extract_speed_model(a = model[[1]], par = 'low'),
                               speed_model_upr = extract_speed_model(a = model[[1]], par = 'high'),
                               
                               # find home range estimate
                               hr_units = "square kilometers",
                               hr_est_50 = extract_hr(a = akde[[1]], par='est', l.ud=0.50),
                               hr_lwr_50 = extract_hr(a = akde[[1]], par='low', l.ud=0.50),
                               hr_upr_50 = extract_hr(a = akde[[1]], par='high', l.ud=0.50),
                               hr_est_95 = extract_hr(a = akde[[1]], par='est', l.ud=0.95),
                               hr_lwr_95 = extract_hr(a = akde[[1]], par='low', l.ud=0.95),
                               hr_upr_95 = extract_hr(a = akde[[1]], par='high', l.ud=0.95),
                               
                               d.length = as.numeric(max(d$timestamp) - min(d$timestamp)),
                               
                               distance_est = speed_model_est * d.length,
                               distance_lwr = speed_model_lwr * d.length,
                               distance_upr = speed_model_upr * d.length
                               
                             )
                             
                          } else {
                            res.tb <- tibble(
                              ind.id = NA_real_,
                              # create variogram
                              variogram = list('Insufficient data.'),
                              # find initial guesses for models
                              guess = list('Insufficient data.'),
                              # select best model based on subset of tel
                              model = list('Insufficient data.'),
                              # estimate autocorrelated kernel density estimate
                              akde = list('Insufficient data.'),
                              # estimate Gaussian speed
                              speed = list('Insufficient data.'),
                              # extract diffusion estimates
                              diffusion_units = NA_real_,
                              diffusion_est = NA_real_,
                              diffusion_lwr = NA_real_,
                              diffusion_upr = NA_real_,
                              # extract tau p estimates
                              tau_p_units = NA_real_,
                              tau_p_est = NA_real_,
                              tau_p_lwr = NA_real_,
                              tau_p_upr = NA_real_,
                              # extract tau v estimates
                              tau_v_units = NA_real_,
                              tau_v_est = NA_real_,
                              tau_v_lwr = NA_real_,
                              tau_v_upr = NA_real_,
                              # extract speed estimates - Gaussian
                              speed_gauss_units = NA_real_,
                              speed_gauss_est = NA_real_,
                              speed_gauss_lwr = NA_real_,
                              speed_gauss_upr = NA_real_,
                              # extract speed from model
                              speed_model_units = NA_real_,
                              speed_model_est = NA_real_,
                              speed_model_lwr = NA_real_,
                              speed_model_upr = NA_real_,
                              # find home range estimate
                              hr_units = NA_real_,
                              hr_est_50 = NA_real_,
                              hr_lwr_50 = NA_real_,
                              hr_upr_50 = NA_real_,
                              hr_est_95 = NA_real_,
                              hr_lwr_95 = NA_real_,
                              hr_upr_95 = NA_real_,
                              d.length = NA_real_,
                              #d.length.int = NA_real_,
                              distance_est = NA_real_,
                              distance_lwr = NA_real_,
                              distance_upr = NA_real_)
                          } # close else
                           
                           
                           
  } # end of do parallel - foreach
  
  stopImplicitCluster()
  
  models <- data.table::rbindlist(models.list) %>%
    as_tibble()
  
  out <- tibble(
    t_start = times, 
    t_end = t_start + window, 
    dataset = purrr::map2(t_start, t_end, function(t_1, t_2) tel[tel$t >= t_1 & tel$t <= t_2, ]),
    models = models) %>%
    tidyr::unnest(models) %>%
    mutate(t_center = (t_start + t_end) / 2,
           posixct = as.POSIXct(t_center, origin = '1970-01-01',
                                tz = tel@info$timezone),
           date = as.Date(posixct), 
           t_start_date = as.POSIXct(t_start, origin = '1970-01-01',
                                     tz = tel@info$timezone),
           t_end_date = as.POSIXct(t_end, origin = '1970-01-01',
                                   tz = tel@info$timezone))
  
  if(!is.null(rds_path)) {
    saveRDS(out, file.path(rds_path,
                           paste0(tel@info['identity'],
                                  '-window-', window / (1 %#% 'day'), '-days',
                                  '-dt-', dt / (1 %#% 'day'), '-days.rds')))
    writexl::write_xlsx(out %>% dplyr::select(-dataset, -guess, -variogram, -model, -akde, -speed) %>%
                          dplyr::select(ind.id, everything()),  
                        file.path(rds_path,
                                  paste0(tel@info['identity'],
                                         '-window-', window / (1 %#% 'day'), '-days',
                                         '-dt-', dt / (1 %#% 'day'), '-days.xlsx')))
    
  }
  
  # plot results
  library('ggplot2') # for fancy plots
  theme_set(theme_bw() + theme(legend.position = 'none'))
  Sys.setlocale("LC_ALL", "English")
  
  # tracking data
  plt_a <-
    ggplot(tel) +
    coord_equal() +
    geom_point(aes(longitude, latitude, color = timestamp)) +
    geom_path(aes(longitude, latitude), alpha = 0.1) +
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
    scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 week",
                 date_labels = "%b %Y") +
    labs(y = expression(Home~range~(km^2))) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
  
  if(plot_50_q) {
    plt_b <- plt_b +
      # core home range 95% CIs
      geom_ribbon(aes(date, ymin = hr_lwr_50, ymax = hr_upr_50), alpha = 0.3) +
      # core home range estimates
      geom_line(aes(date, hr_est_50), linewidth = 1.25) +
      geom_line(aes(date, hr_est_50, color = posixct))
  }
  
  plt <- cowplot::plot_grid(plt_a, plt_b, labels = c('a.', 'b.'), nrow = 1,
                            align = 'hv')
  
  if (! is.null(fig_path)) {
    # Save figure as a png using the animal's name
    file.path(fig_path, paste0(tel@info['identity'],
                               '-window-', window / (1 %#% 'day'),
                               '-days-dt-', dt / (1 %#% 'day'),
                               '-days.png')) %>%
      ggsave(plot = plt, units = 'in', width = 18, height = 8,
             dpi = 600, bg = 'white')
  } else {
    print(plt)
  }
  
}

# ------------------------------------------------------------------------------




