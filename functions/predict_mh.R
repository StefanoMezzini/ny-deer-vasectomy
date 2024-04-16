#' to make predictions using the Metropolis Hastings sampler via `gam.mh()` 
library('mgcv')    #' for modeling and `gam.mh()`
library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('ggplot2') # for fancy plots
theme_set(theme_bw())

RUN_TEST <- FALSE

# plot the estimated common trends using the MH sampler ----
lp_mat_mh <- function(model, newdata = NULL, terms = NULL, exclude = NULL,
                      quantiles = c(0.025, 0.975),
                      median = TRUE, n_samples = 1e4, burn_in = 1e3,
                      t_df = 40, thin_by = 0, rw_scale = 0.25, ...) {
  
  if(is.null(newdata)) newdata <- model$model[, -1] # drop response column
  
  mh_samples <- gam.mh(model, ns = n_samples, burn = burn_in,
                       t.df = t_df, thin = thin_by, rw.scale = rw_scale)
  
  # estimated effective sample size
  eess <- coda::effectiveSize(coda::as.mcmc(mh_samples$bs))
  
  p_hist <-
    ggplot() +
    geom_histogram(aes(x = eess), bins = 10, color = 'black',
                   fill = 'grey') +
    labs(x ='Estimated sample size', y = 'Count')
  
  # check effective sample sizes for each coefficient
  p_chains <- 
    as.data.frame(mh_samples$bs) %>%
    mutate(iter = 1:n()) %>%
    pivot_longer(cols = ! iter) %>%
    # z-transform the coefficients
    group_by(name) %>%
    mutate(value = (value - mean(value)) / sd(value)) %>%
    ungroup() %>%
    filter(name == first(name)) %>%
    ggplot(aes(iter, value, group = name)) +
    geom_line(alpha = 0.2) +
    geom_point(alpha = 0.2) +
    geom_smooth(color = 'darkorange', method = 'gam', formula = y ~ s(x)) +
    labs(x = 'Iteration (after any thinning)',
         y = 'Z-transformed coefficient estimate')
  plot(cowplot::plot_grid(p_hist, p_chains))
  
  X <- predict(object = model, newdata = newdata, type = 'lpmatrix',
               terms = terms, exclude = exclude, ...)
  
  preds <- X %*% t(mh_samples$bs) %>%
    as.data.frame() %>%
    mutate(row = 1:n()) %>%
    pivot_longer(! row, names_to = 'sim', values_to = 'samples') %>%
    group_by(row) %>%
    summarize(mu_lwr = quantile(samples, min(quantiles)),
              mu_est = quantile(samples, 0.5),
              mu_upr = quantile(samples, max(quantiles))) %>%
    mutate(across(c(mu_lwr, mu_est, mu_upr), model$family$linkinv)) %>%
    bind_cols(newdata, .)
  
  return(preds)
}

# test the function
if(RUN_TEST) {
  set.seed(1)
  N <- 500
  d <- tibble(x1 = runif(N),
              x2 = rgamma(N, 3),
              x3 = sample(letters, N, replace = TRUE) %>%
                factor(),
              y = rpois(n = N, exp(5 - x1 * 4 - cospi(x2) * 2)))
  
  layout(t(1:2))
  plot(y ~ x1, d)
  plot(y ~ x2, d)
  
  m <- gam(y ~ s(x1) + s(x2) + s(x3, bs = 're'),
           family = poisson(link = 'log'),
           data = d,
           method = 'REML')
  
  plot(m, pages = 1, scheme = 1)
  layout(matrix(1:4, ncol = 2))
  gam.check(m)
  layout(1)
  
  preds <- predict(m, type = 'link', se.fit = TRUE,
                   terms = c('s(x1)', '(Intercept)')) %>%
    as.data.frame() %>%
    bind_cols(d, .) %>%
    mutate(mu_lwr = exp(fit - 1.96 * se.fit),
           mu_est = exp(fit),
           mu_upr = exp(fit + 1.96 * se.fit))
  
  set.seed(1)
  preds_mh <-
    lp_mat_mh(m,
              terms = c('s(x1)', '(Intercept)'),
              n_samples = 1e3, #' number of samples excluding `burn_in`
              burn_in = 100,   #' n samples to discard (can be low)
              t_df = 40,       #' t df for initial proposal
              thin_by = 1,     #' take every `thin_by`th sample
              rw_scale = 0.25)  #' random walk scale
  
  ggplot() +
    geom_point(aes(x1, y), d, alpha = 0.5) +
    geom_ribbon(aes(x1, ymin = mu_lwr, ymax = mu_upr, fill = 'mh'),
                preds_mh, alpha = 0.3) +
    geom_line(aes(x1, mu_est, color = 'mh'), preds_mh) +
    geom_ribbon(aes(x1, ymin = mu_lwr, ymax = mu_upr, fill = 'gaus'),
                preds, alpha = 0.3) +
    geom_line(aes(x1, mu_est, color = 'gaus'), preds) +
    scale_color_brewer('Method', type = 'qual', palette = 6,
                       aesthetics = c('color', 'fill')) +
    ylab('Response')
}

