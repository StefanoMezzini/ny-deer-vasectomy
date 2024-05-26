#' this code was originally written by Gavin Simpson for simulations from
#' Gamma location-scale (H)GAMs `gammals`. I (Stefano Mezzini) have
#' adapted the code to work with the `betals` family of distributions that
#' Simon Wood coded for our lab. However, you may want to look at the 
#' `{gratia}` package that Gavin maintains, since he intended to include
#' all of these functions in the package when I wrote this. The functions
#' for derivatives for generic location (H)GAMs aleady exist, but he was
#' working on writing generic functions for location-scale families, too.
#' *original code*: https://github.com/simpson-lab/wpg-mb-lakes/blob/main/analysis/variance-simulation-and-derivatives.R

#' *NOTE:* while the `betals()` function says the link functions are the
#'         identity functions, they are actually logit functions for both
#'         mu and phi

EXAMPLES <- FALSE # should examples be run?

##' Simulate from the posterior distribution of the mean of a Beta LS GAM
##'  using Gaussian approximation to the posterior
##'
##' @param model         the fitted GAM
##' @param data          the new data locations you want to get variance for
##' @param nsims         the number of posterior draws wanted - low by default
##'                      to avoid excessive computation, but needs to be
##'                      10,000+ for quantile-based intervals
##' @param unconditional logical; use the smoothness selection corrected version
##'                      of the Bayesian covariance matrix of the model?
`betals_mean` <- function(model, data, nsims = 100, unconditional = FALSE,
                          exclude = NULL) {
  ## Simulate variance from posterior
  sim <- sim_betals_mean(model = model, data = data, nsims = nsims,
                         unconditional = unconditional, exclude = exclude)
  ## process results into a tibble
  colnames(sim) <- paste0("sim", seq_len(nsims))
  tbl <- as_tibble(sim) %>%
    bind_cols(data) %>%
    tibble::rowid_to_column(var = "row")
  tbl <- pivot_longer(tbl,
                      cols = matches("^sim"),
                      names_to = "simulation",
                      values_to = "mean")
  tbl
}

##' Simulate from the posterior distribution of the variance of a Beta LS GAM
##'  using Gaussian approximation to the posterior
##'
##' @param model         the fitted GAM
##' @param data          the new data locations you want to get variance for
##' @param nsims         the number of posterior draws wanted - low by default
##'                      to avoid excessive computation, but needs to be
##'                      10,000+ for quantile-based intervals
##' @param unconditional logical; use the smoothness selection corrected version
##'                      of the Bayesian covariance matrix of the model?
`betals_var` <- function(model, data, nsims = 100, unconditional = FALSE,
                         exclude = NULL) {
  ## Simulate variance from posterior
  sim <- sim_betals_var(model = model, data = data,
                        nsims = nsims, unconditional = unconditional,
                        exclude = exclude)
  ## process results into a tibble
  colnames(sim) <- paste0("sim", seq_len(nsims))
  tbl <- as_tibble(sim) %>%
    bind_cols(data) %>%
    tibble::rowid_to_column(var = "row")
  tbl <- pivot_longer(tbl,
                      cols = matches("^sim"),
                      names_to = "simulation",
                      values_to = "variance")
  tbl
}

##' The internal workhorse does all the cool stuff 
`sim_betals_mean` <- function(model, data, nsims = 100,
                              unconditional = FALSE, exclude = NULL) {
  ## prediction matrix
  Xp <- predict(model, newdata = data, type = 'lpmatrix', exclude = exclude)
  ## model parameters
  coefs <- coef(model)
  ## Bayesian covariance matrix
  Vb <- vcov(model, unconditional = unconditional)
  ## which coefs go with the phi linear predictor
  phi_take <- grepl('^s\\.1', colnames(Xp)) |
    colnames(Xp) %in% c('(Intercept).1')
  
  ## Simulate from posterior using Gaussian approximation
  betas <- mvnfast::rmvn(n = nsims,
                         mu = coefs,
                         sigma = Vb)
  
  ## simplify later code so form the compliment to select mean
  ## linear predictor
  mu_take <- !phi_take
  ## subset Xp matrix into mean and phi parts
  Xp_mu <- Xp[, mu_take, drop = FALSE]
  
  ## Predict for mean
  fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
  ilink_mu <- brms:::inv_logit
  fit_mu <- ilink_mu(fit_mu) # apply g-1() = inv_logit
  fit_mu
}

##' The internal workhorse does all the cool stuff 
`sim_betals_var` <- function(model, data, nsims = 100,
                             unconditional = FALSE, exclude = NULL) {
  ## prediction matrix
  Xp <- predict(model, newdata = data, type = 'lpmatrix', exclude = exclude)
  ## model parameters
  coefs <- coef(model)
  ## Bayesian covariance matrix
  Vb <- vcov(model, unconditional = unconditional)
  ## which coefs go with the phi linear predictor
  phi_take <- grepl('^s\\.1', colnames(Xp)) |
    colnames(Xp) %in% c('(Intercept).1')
  
  ## Simulate from posterior using Gaussian approximation
  betas <- mvnfast::rmvn(n = nsims,
                         mu = coefs,
                         sigma = Vb)
  
  ## simplify later code so form the compliment to select mean
  ## linear predictor
  mu_take <- !phi_take
  ## subset Xp matrix into mean and phi parts
  Xp_mu <- Xp[, mu_take, drop = FALSE]
  Xp_phi <- Xp[, phi_take, drop = FALSE]
  
  ## Predict for mean
  fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
  ilink_mu <- brms:::inv_logit
  fit_mu <- ilink_mu(fit_mu) # apply g-1() = inv_logit
  
  ## Predict for phi
  fit_phi <- Xp_phi %*%
    t(betas[, phi_take, drop = FALSE]) # predict on internal scale
  ilink_phi <- brms:::inv_logit
  fit_phi <- ilink_phi(fit_phi) # apply g-1() = inv_logit
  
  ##' variance is `sigma2 = mu * (1 - mu) * phi` where:
  ##' `mu` is the mean, defined as `mu = a / (a + b)`, and
  ##' `phi` is the scale parameter, defined as `phi = 1 / (a + b + 1)`,
  ##' such that
  ##' `sigma2 = mu * (1 - mu) * phi = a/(a+b) * b/(a+b) * 1/(1+a+b) =`
  ##' `       = (a * b) / ((a+ b)^2 * (1 + a + b))`,
  ##' where `a` and `b` are the two shape parameters of a Beta distribution
  ##' see: *https://en.wikipedia.org/wiki/Beta_distribution*
  fit_var_draws <- fit_mu * (1 - fit_mu) * fit_phi
  ## return
  fit_var_draws
}

## Examples -----------------------------------------------------------

if(EXAMPLES) {
  
  library('dplyr')   # for data wrangling
  library('tidyr')   # for data wrangling
  library('ggplot2') # for fancy plots
  library('mgcv')    # for GAMs
  library('gratia')  # for visualizing GAMs
  theme_set(theme_bw())
  source('functions/betals.r')
  
  ## number of simulations
  K <- 25
  
  # simulate some data with known mean and variance
  set.seed(1)
  d <- tibble(t = seq(0, 25, by = 0.001),
              a = sin(t / 2) * 20 + 30,
              b = sin(t / 3) * 20 + 25,
              mu = a / (a + b),
              phi = 1 / (1 + a + b),
              sigma2 = a * b / ((a + b)^2 * (1 + a + b)),
              y = rbeta(length(t), shape1 = a, shape2 = b))
  
  # ensure parameterization is correct
  with(d, all(round(mu * (1 - mu) * phi, 10) == round(sigma2, 10)))
  
  # preview the data
  pivot_longer(d, ! t) %>%
    ggplot(aes(t, value)) +
    facet_wrap(~ name, scales = 'free_y') +
    geom_line()
  
  # fit a Beta model
  m_beta <- gam(y ~ s(t),
                family = betar(),
                data = d,
                method = 'REML')
  
  # new data for predictions
  newd <- select(d, ! y)
  
  # predictions are reasonable but off because they assume a constant phi
  mutate(newd,
         mu_hat = predict(m_beta, newdata = newd, type = 'response')) %>%
    ggplot() +
    geom_point(aes(t, y), d, alpha = 0.2) +
    geom_line(aes(t, mu_hat), color = 'cyan', linewidth = 1) +
    geom_line(aes(t, mu), d, color = 'darkorange', linewidth = 1)
  
  # fit a betals model
  m_betals <- gam(list(y ~ s(t, k = 20),
                       ~ s(t, k = 20)),
                  family = betals(),
                  data = d,
                  method = 'REML')
  
  preds_betals <-
    bind_cols(newd,
              predict(m_betals, newdata = newd, type = 'response') %>%
                as.data.frame() %>%
                transmute(mu_hat = V1,
                          phi_hat = V2,
                          sigma2_hat = mu_hat * (1 - mu_hat) * phi_hat))
  
  # simulate for mean
  set.seed(1)
  mu_sim <- betals_mean(m_betals, data = newd, nsims = K)
  
  # plot the simulations for the mean
  ggplot(mu_sim) +
    geom_line(aes(t, mean, group = simulation), mu_sim, alpha = 0.2,
              linewidth = 1) +
    geom_line(aes(t, mu), d, color = 'darkorange', linewidth = 1) +
    geom_line(aes(t, mu_hat), preds_betals, color = 'cyan', linewidth = 1)
  
  # simulations from posterior of variance of a Beta LS model
  set.seed(1)
  var_sim <- betals_var(m_betals, data = newd, nsims = K)
  
  # plot the simulations for the variance
  ggplot() +
    geom_line(aes(t, variance, group = simulation), var_sim, alpha = 0.2,
              linewidth = 1) +
    geom_line(aes(t, sigma2), d, color = 'darkorange', linewidth = 1) +
    geom_line(aes(t, sigma2_hat), preds_betals, color = 'cyan', linewidth = 1)
  
  ## derivatives of Beta LS mean using Gaussian approximation
  ## note these are all being done on the mean scale itself and only via
  ## left finite differences
  set.seed(1)
  mu_d <- betals_mean_deriv(m_betals, data = newd, nsims = K,
                            var = 't', eps = 1e-5)
  
  ## plot derivatives
  cowplot::plot_grid(
    ggplot(mu_sim, aes(x = t, y = mean, group = simulation)) +
      geom_line(alpha = 0.2),
    ggplot(mu_d, aes(x = t, y = derivative, group = simulation)) +
      geom_area(alpha = 0.2) +
      geom_hline(yintercept = 0),
    ncol = 1)
  
  ## derivatives of Beta LS variance using Gaussian approximation
  ## note these are all being done on the variance scale itself ad only via
  ## left finite differences
  var_d <- betals_var_deriv(m_betals, data = newd, nsims = K,
                            var = 't', eps = 1e-5)
  
  ## plot derivatives
  cowplot::plot_grid(
    ggplot(var_sim, aes(x = t, y = variance, group = simulation)) +
      geom_line(alpha = 0.2),
    ggplot(var_d, aes(x = t, y = derivative, group = simulation)) +
      geom_area(alpha = 0.2) +
      geom_hline(yintercept = 0),
    ncol = 1)
}
