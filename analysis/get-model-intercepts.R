library('dplyr') # for data wrangling
library('mgcv')  # for GAMs

get_cis <- function(m) {
  if(grepl('density', m)) {
    l_inv <- brms:::inv_logit
  } else {
    l_inv <- exp
  }
  
  m <- readRDS(m)
  
  tibble(sex_treatment = unique(m$model$sex_treatment),
         days_since_aug_1 = 0,
         study_year = 1,
         animal_year = m$model$animal_year[1]) %>%
    bind_cols(.,
              predict(m, newdata = ., se.fit = TRUE, discrete = FALSE,
                      terms = c('(Intercept)', 'sex_treatment')) %>%
                as.data.frame()) %>%
    mutate(lwr_95 = l_inv(fit.1 - 1.96 * se.fit.1),
           mu = l_inv(fit.1),
           upr_95 = l_inv(fit.1 + 1.96 * se.fit.1)) %>%
    select(sex_treatment, lwr_95, mu, upr_95)
}

# extract model summaries
get_cis('models/m_hr-hgamls.rds')
get_cis('models/m_diffusion-hgamls.rds')
get_cis('models/m_speed-hgamls.rds')
get_cis('models/m_density-hgamls.rds')
