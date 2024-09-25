#' this code was originally written by Gavin Simpson for simulations.
#' You may want to look at the `{gratia}` package that Gavin maintains,
#' since he intended to include all of these functions in the package when
#' I edited this. The functions for derivatives for generic location
#' (H)GAMs aleady exist, but he was working on writing generic functions
#' for location-scale families, too.
#' *original code*: https://github.com/simpson-lab/wpg-mb-lakes/blob/main/analysis/variance-simulation-and-derivatives.R

##' Simulate from the posterior distribution of the mean of a Gamma LS GAM
##'  using Gaussian approximation to the posterior
##'
##' @param model         the fitted GAM
##' @param data          the new data locations you want to get variance for
##' @param nsims         the number of posterior draws wanted - low by default
##'                      to avoid excessive computation, but needs to be
##'                      10,000+ for quantile-based intervals
##' @param unconditional logical; use the smoothness selection corrected version
##'                      of the Bayesian covariance matrix of the model?
`gammals_mean` <- function(model, data, nsims = 100,
                         unconditional  = FALSE, exclude = NULL, ...) {
    ## Simulate variance from posterior
    sim <- sim_gammals_mean(model = model, data = data,
                            nsims = nsims, unconditional = unconditional,
                            exclude = exclude, ...)
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

##' Simulate from the posterior distribution of the variance of a Gamma LS GAM
##'  using Gaussian approximation to the posterior
##'
##' @param model         the fitted GAM
##' @param data          the new data locations you want to get variance for
##' @param nsims         the number of posterior draws wanted - low by default
##'                      to avoid excessive computation, but needs to be
##'                      10,000+ for quantile-based intervals
##' @param unconditional logical; use the smoothness selection corrected version
##'                      of the Bayesian covariance matrix of the model?
`gammals_var` <- function(model, data, nsims = 100,
                          unconditional  = FALSE, exclude = NULL, ...) {
    ## Simulate variance from posterior
    sim <- sim_gammals_var(model = model, data = data,
                           nsims = nsims, unconditional = unconditional,
                           exclude = exclude, ...)
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
`sim_gammals_mean` <- function(model, data, nsims = 100, check_data = TRUE,
                               unconditional = FALSE, exclude = NULL) {
    ## prediction matrix
    Xp <- predict(model, newdata = data, type = 'lpmatrix', exclude = exclude,
                  newdata.guaranteed = check_data)
    ## model parameters
    coefs <- coef(model)
    ## Bayesian covariance matrix
    Vb <- vcov(model, unconditional = unconditional)
    ## which coefs go with the theta linear predictor
    theta_take <- grepl('^s\\.1', colnames(Xp)) |
        colnames(Xp) %in% c('(Intercept).1',
                            'sex_treatmentf staten_island.1',
                            'sex_treatmentm rockefeller.1',
                            'sex_treatmentm staten_island.1')

    ## Simulate from posterior using Gaussian approximation
    betas <- mvnfast::rmvn(n = nsims,
                           mu = coefs,
                           sigma = Vb)
    
    ## simplify later code so form the compliment to select mean
    ## linear predictor
    mu_take <- !theta_take
    ## subset Xp matrix into mean and theta parts
    Xp_mu <- Xp[, mu_take, drop = FALSE]

    ## Predict for mean
    fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
    ilink_mu <- inv_link(model, parameter = "location") # link function
    fit_mu <- ilink_mu(fit_mu) # apply g-1() this is just identity so redundant
    ## This model is parameterised in terms of the log-mean, so we still need to
    ## transform to the actual data scale using exp()
    fit_mu <- exp(fit_mu)
    fit_mu
}

##' The internal workhorse does all the cool stuff 
`sim_gammals_var` <- function(model, data, nsims = 100, check_data = TRUE,
                              unconditional = FALSE, exclude = NULL) {
    ## prediction matrix
    Xp <- predict(model, newdata = data, type = 'lpmatrix', exclude = exclude,
                  newdata.guaranteed = check_data)
    ## model parameters
    coefs <- coef(model)
    ## Bayesian covariance matrix
    Vb <- vcov(model, unconditional = unconditional)
    ## which coefs go with the theta linear predictor
    theta_take <- grepl('^s\\.1', colnames(Xp)) |
      colnames(Xp) %in% c('(Intercept).1',
                          'sex_treatmentf staten_island.1',
                          'sex_treatmentm rockefeller.1',
                          'sex_treatmentm staten_island.1')
    
    ## Simulate from posterior using Gaussian approximation
    betas <- mvnfast::rmvn(n = nsims,
                           mu = coefs,
                           sigma = Vb)
    
    ## simplify later code so form the compliment to select mean
    ## linear predictor
    mu_take <- !theta_take
    ## subset Xp matrix into mean and theta parts
    Xp_mu <- Xp[, mu_take, drop = FALSE]
    Xp_theta <- Xp[, theta_take, drop = FALSE]

    ## Predict for mean
    fit_mu <- Xp_mu %*% t(betas[, mu_take, drop = FALSE]) # predict on internal scale
    ilink_mu <- inv_link(model, parameter = "location") # link function
    fit_mu <- ilink_mu(fit_mu) # apply g-1() this is just identity so redundant
    ## This model is parameterised in terms of the log-mean, so we still need to
    ## transform to the actual data scale using exp()
    fit_mu <- exp(fit_mu)
    
    ## Predict for theta
    fit_theta <- Xp_theta %*%
        t(betas[, theta_take, drop = FALSE]) # predict on internal scale
    ilink_theta <- inv_link(model, parameter = "scale") # theta inverse link
    fit_theta <- ilink_theta(fit_theta) # apply g-1()
    ## fit scale even after using inverse link is log for theta parameter, so
    ## more back transforming
    fit_theta <- exp(fit_theta)

    ## variance is mu * s where s is the scale in sense of rgamma, not theta
    ## From ?rgamma Var(y) = shape * scale^2 = (1 / theta) * (mu * theta)^2
    ## From ?gammals Var(y) = mu * scale = mu * s
    ##   where scale = s = mu * theta. Hence from ?gammals we arrive finally
    ##   at: Var(y) = mu * s = mu * (mu * theta)
    fit_var_draws <- fit_mu * (fit_mu * fit_theta)
    ## return
    fit_var_draws
}

# test the function
if(FALSE) {
  library(dplyr)
  library(tidyr)
  library(mgcv)
  library(gratia)
  
  d <- tibble(x = seq(0, 1, by = 0.001),
              mu = sinpi(2 * x) + 3,
              scale = cospi(2 * x) + 2,
              var = mu * (mu * scale),
              y = rgamma(n = length(x),
                         shape = mu^2 / var,
                         scale = var / mu))
  
  plot(mu ~ x, d)
  plot(scale ~ x, d)
  plot(var ~ x, d)
  plot(y ~ x, d)
  
  m <- gam(list(y ~ s(x), ~ s(x)), family = gammals(), data = d,
           method = 'REML')
  plot(m, pages = 1, scheme = 1)
  
  preds_mu <-
    gammals_mean(m, data = d, nsims = 1e4, unconditional = TRUE) %>%
    group_by(x) %>%
    summarize(mean = median(mean))
  
  preds_s2 <-
    gammals_var(m, data = d, nsims = 1e4, unconditional = TRUE) %>%
    group_by(x) %>%
    summarize(s2 = median(variance))
  
  plot(mean ~ x, preds_mu)
  plot(s2 ~ x, preds_s2)
  
  plot(preds_mu$mean ~ d$mu)
  plot(preds_s2$s2 ~ d$var)
}
