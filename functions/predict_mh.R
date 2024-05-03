#' to make predictions using the Metropolis Hastings sampler via `gam.mh()` 
library('mgcv')    #' for modeling and `gam.mh()`
library('gratia')  # for convenient functions to work with GAMs
library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('ggplot2') # for fancy plots
theme_set(theme_bw())

RUN_TEST <- FALSE
DIAGNOSTICS <- TRUE

# plot the estimated common trends using the MH sampler ----
predict_mh <- function(model, newdata = NULL, terms = NULL, exclude = NULL,
                       quantiles = c(0.025, 0.975), n_samples = 1e4,
                       burn_in = 1e3, t_df = 40, rw_scale = 0.25,
                       thin_by = 0, chain_diagnostics = TRUE,
                       plot_smooths = FALSE, ...) {
  
  if(! is.null(terms) & ! is.null(exclude)) {
    stop('Cannot have both non-NULL `terms` and `exclude.')
  }
  
  if(is.null(newdata)) newdata <- model$model[, -1] # drop response column
  
  if('linkinv' %in% names(model$family)) {
    
    ls_family <- FALSE # not location-scale
    ilink_mu <- inv_link(model, parameter = 'location')
    i_mu <- rep(TRUE, length(coef(model))) # take all coefficients
    
    # find linear predictor matrix
    Xp_mu <- predict(object = model, newdata = newdata, type = 'lpmatrix',
                     terms = terms, exclude = exclude, ...)
    
  } else if (model$family$family %in% c('gaulss', 'gammals')) {
    
    #' can change condition to `inherits(b$family, "general.family")` later
    ls_family <- TRUE
    
    ilink_mu <- inv_link(model, parameter = 'location')
    ilink_scale <- inv_link(model, parameter = 'scale')
    
    # find full linear predictor matrix
    Xp <- predict(object = model, newdata = newdata, type = 'lpmatrix',
                  terms = terms, exclude = exclude, ...)
    
    # indices of coefficients for each linear predictor
    i_scale <- grepl('^s\\.1', colnames(Xp)) | # starts with s.1
      colnames(Xp) %in% c('(Intercept).1') # or scale parameter intercept
    i_mu <- ! i_scale # all others (currently only supports location-scale)
    
    # find linear predictor matrix for mean and scale parts
    Xp_mu <- Xp[, i_mu, drop = FALSE]
    Xp_scale <- Xp[, i_scale, drop = FALSE]
    
  } else {
    stop('Family ', model$family$family, ' not yet supported.')
  }
  
  # estimate distribution of coefficients using the M-H sampler
  mh_samples <- gam.mh(model, ns = n_samples, burn = burn_in,
                       t.df = t_df, thin = thin_by, rw.scale = rw_scale)
  
  if(chain_diagnostics) {
    # estimated effective sample size
    eess <- tibble(n = coda::effectiveSize(coda::as.mcmc(mh_samples$bs)),
                   parameter = if_else(i_mu, 'Location', 'Scale'))
    
    p_hist <-
      ggplot(eess) +
      facet_grid(. ~ parameter) +
      geom_histogram(aes(x = n), bins = 10, color = 'black',
                     fill = 'grey') +
      labs(x ='Estimated sample size for each coefficient', y = 'Count')
    
    # check effective sample sizes for each coefficient
    p_chains <- 
      as.data.frame(mh_samples$bs) %>%
      mutate(iter = 1:n()) %>%
      pivot_longer(cols = ! iter) %>%
      mutate(parameter = if_else(! grepl('^s\\.1', name) |
                                   name %in% c('(Intercept).1'),
                                 'Location', 'Scale')) %>%
      # z-transform the coefficients
      group_by(name) %>%
      mutate(value = (value - mean(value)) / sd(value)) %>%
      ungroup() %>%
      ggplot(aes(iter, value, group = name)) +
      facet_grid(. ~ parameter) +
      geom_line(alpha = 0.05) +
      labs(x = 'Iteration (after warmup and thinning)',
           y = 'Z-transformed coefficient estimate')
    
    if(plot_smooths) {
      p_chains <- p_chains +
        geom_smooth(color = '#ff8c0040', method = 'gam',
                    formula = y ~ s(x), se = FALSE)
    }
    
    plot(cowplot::plot_grid(p_hist, p_chains, ncol = 1))
  }
  
  if(! ls_family) {
    preds_mu <-
      Xp_mu %*% t(mh_samples$bs[, i_mu]) %>% # calculate predictions
      as.data.frame() %>% # convert matrix to data frame
      # pivot to long format
      mutate(row = 1:n()) %>%
      pivot_longer(! row, names_to = 'sim', values_to = 'mean_samples') %>%
      mutate(mean_samples = ilink_mu(mean_samples)) %>%
      # calculate quantiles and mean
      group_by(row) %>%
      summarize(
        mu_lwr = quantile(mean_samples, min(quantiles), na.rm = TRUE),
        mu_est = quantile(mean_samples, 0.5, na.rm = TRUE),
        mu_upr = quantile(mean_samples, max(quantiles), na.rm = TRUE)) %>%
      bind_cols(newdata, .) %>%
      as_tibble()
    
    return(preds_mu)
  } else {
    samples_mu <-
      Xp_mu %*% t(mh_samples$bs[, i_mu]) %>% # calculate predictions
      as.data.frame() %>% # convert matrix to data frame
      # pivot to long format
      mutate(row = 1:n()) %>%
      pivot_longer(! row, names_to = 'sim', values_to = 'mean_samples')
    
    samples_scale <-
      Xp_scale %*% t(mh_samples$bs[, i_scale]) %>% # calculate predictions
      as.data.frame() %>% # convert matrix to data frame
      # pivot to long format
      mutate(row = 1:n()) %>%
      pivot_longer(! row, names_to = 'sim', values_to = 'scale_samples')
    
    preds <- full_join(samples_mu, samples_scale, by = c('row', 'sim'))
    
    if(model$family$family == 'gaulss') {
      # convert SD to variance (independent of mean)
      get_var <- function(mu, scale) scale^2
    } else if(model$family$family == 'gammals') {
      get_var <- function(mu, scale) mu * (mu * scale) 
    } else {
      stop('Family', model$family$family, ' not supported.')
    }
    
    preds <-
      preds %>%
      group_by(row) %>%
      summarize(
        # estimated mean and corresponding CIs
        mu_lwr = quantile(mean_samples, min(quantiles), na.rm = TRUE),
        mu_est = quantile(mean_samples, 0.5, na.rm = TRUE),
        mu_upr = quantile(mean_samples, max(quantiles), na.rm = TRUE),
        # estimated variance and corresponding CIs
        scale_lwr = quantile(scale_samples, min(quantiles), na.rm = TRUE),
        scale_est = quantile(scale_samples, 0.5, na.rm = TRUE),
        scale_upr = quantile(scale_samples, max(quantiles), na.rm = TRUE)) %>%
      mutate(s2_lwr = get_var(mu_lwr, scale_lwr),
             s2_est = get_var(mu_est, scale_est),
             s2_upr = get_var(mu_upr, scale_upr)) %>%
      bind_cols(newdata, .)
    
    return(preds)
  }
}

# test the function ----
if(RUN_TEST) {
  #' coefs are less Gaussian more with low `N`; try `N in {10, 100, 1000}`
  N <- 1000
  Q <- c(0.25, 0.75)
  d <- tibble(x1 = runif(N),
              x2 = rgamma(N, 3)) %>%
    mutate(mu = exp(5 - x1 * 4 - cospi(x2)^2 * 2),
           y = rgamma(n = N, shape = (mu) * x2, rate = 1 / x2))
  
  if(DIAGNOSTICS) {
    layout(t(1:2))
    plot(y ~ x1, d)
    plot(y ~ x2, d)
  }
  
  m <- gam(y ~ s(x1, k = 10),# + s(x2, k = 10),
           family = Gamma(link = 'log'),
           data = d,
           method = 'REML')
  
  if(DIAGNOSTICS) {
    plot(m, pages = 1, scheme = 1, scale = 0)
    layout(matrix(1:4, ncol = 2))
    gam.check(m, pch = 1)
    abline(0, 1, col = 'red', lwd = 2)
    layout(1)
  }
  
  newd <- tibble(x1 = seq(0, 1, length.out = 400),
                 x2 = seq(min(d$x2), max(d$x2), length.out = 400),
                 x3 = 'new level')
  
  preds <- predict(m, newdata = newd, type = 'link', se.fit = TRUE,
                   terms = c('s(x1)', '(Intercept)')) %>%
    as.data.frame() %>%
    bind_cols(newd, .) %>%
    mutate(mu_lwr = m$family$linkinv(fit + qnorm(Q[1]) * se.fit),
           mu_est = m$family$linkinv(fit),
           mu_upr = m$family$linkinv(fit + qnorm(Q[2]) * se.fit))
  
  preds_mh <-
    predict_mh(model = m,
               newdata = newd,  #' new data to predict from
               terms = NULL,    #' can predict for specific terms
               exclude = NULL,  #' can exclude specific terms
               quantiles = Q,   #' quantiles for credible intervals
               n_samples = 1e4, #' number of samples excluding `burn_in`
               burn_in = 100,   #' n samples to discard (can be low)
               t_df = 40,       #' t df for initial proposal
               thin_by = 1,     #' take every `thin_by`th sample
               rw_scale = 0.6,  #' adjust so that RW acceptance ~= 0.25
               chain_diagnostics = DIAGNOSTICS)
  
  p <-
    ggplot() +
    geom_ribbon(aes(x1, ymin = mu_lwr, ymax = mu_upr, fill = 'mh'),
                preds_mh, alpha = 0.3) +
    geom_line(aes(x1, mu_est, color = 'mh'), preds_mh, lty = 'dashed') +
    geom_ribbon(aes(x1, ymin = mu_lwr, ymax = mu_upr, fill = 'gaus'),
                preds, alpha = 0.3) +
    geom_line(aes(x1, mu_est, color = 'gaus'), preds, lty = 'dotted') +
    scale_color_brewer('Method', type = 'qual', palette = 6,
                       aesthetics = c('color', 'fill')) +
    ylab('Response')
  
  if(DIAGNOSTICS) p <- p + geom_point(aes(x1, y), d, alpha = 0.5)
  
  cowplot::plot_grid(gratia::qq_plot(m, point_alpha = sqrt(1/N) * 3), p)
}
 