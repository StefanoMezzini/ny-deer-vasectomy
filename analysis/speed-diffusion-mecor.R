library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('ggplot2')   # for fancy plots
library('lme4')      # for models with random effects but no error calibration
library('mecor')     # for models with measurement error correction
library('ctmm')      # for movement models
library('khroma')    # for colorblind-friendly color palettes
source('analysis/figures/default-theme.R')

if(file.exists('data/log-speed-diffusion.rds')) {
  d <- readRDS('data/log-speed-diffusion.rds')
} else {
  d <-
    readRDS('data/years-1-and-2-data-no-akde.rds') %>%
    filter(is.finite(speed_est)) %>% # loosing ~40% of the windows
    #' `Log` corrects for Jensen's Inequality and returns a list with the
    #' logged values (`log`) & the variance of the logged values (`VAR.log`)
    #' see `ctmm::Exp()` for back-transforming
    bind_cols(
      Log(.$model, variable = 'speed', debias = TRUE, units = FALSE) %>%
        rename(log_speed = log, var_log_speed = VAR.log),
      Log(.$model, variable = 'diffusion', debias = TRUE, units = FALSE) %>%
        rename(log_diffusion = log, var_log_diffusion = VAR.log)) %>%
    select(sex_treatment, study_year, animal, Age, date,
           log_speed, var_log_speed, dof_speed,
           log_diffusion, var_log_diffusion, dof_diff) %>%
    mutate(dt_h = if_else(animal %in% c('22', '33', '38', '41', '57', '91'),
                          4, 1))
  saveRDS(d, 'data/log-speed-diffusion.rds')
}
d

ggplot(mutate(d, dt_h = factor(dt_h))) +
  facet_grid(dt_h ~ .) +
  # vertical bars for 1 SD
  geom_errorbar(aes(x = log_diffusion,
                    ymin = log_speed - sqrt(var_log_speed),
                    ymax = log_speed + sqrt(var_log_speed),
                    color = dt_h),
                width = 0, alpha = 0.05) +
  # horizontal bars for 1 SD
  geom_errorbarh(aes(xmin = log_diffusion - sqrt(var_log_diffusion),
                     xmax = log_diffusion + sqrt(var_log_diffusion),
                     y = log_speed, color = dt_h),
                 height = 0, alpha = 0.05) +
  # data points
  geom_point(aes(log_diffusion, log_speed, color = dt_h), alpha = 0.3) +
  # models
  geom_smooth(aes(log_diffusion, log_speed), color = 'black',
              method = 'lm', formula = y ~ x) +
  labs(x = 'log(diffusion)', y = 'log(speed)') +
  scale_color_brewer(name = 'Sampling interval (hours)', type = 'qual',
                     palette = 6)

# fit the naive model with lme4 ----
m_lme4 <- lmer(log_speed ~ 1 + log_diffusion + dt_h +
                 (1 | animal) + (log_diffusion | animal),
               data = d, REML = TRUE)

# fit the error-corrected model ----
m_mecor <- mecor(log_speed ~ dt_h +
                   MeasErrorRandom(substitute = log_diffusion,
                                   variance = mean(d$var_log_diffusion)),
                 data = d,
                 method = "standard",
                 B = 1e4)

# histograms and qqplots of bootstrap outputs for each coefficient
layout(1:3)
hist(m_mecor$corfit$boot$coef[, '(Intercept)'], breaks = 100)
hist(m_mecor$corfit$boot$coef[, 'cor_log_diffusion'], breaks = 100)
hist(m_mecor$corfit$boot$coef[, 'dt_h'], breaks = 100)

qqnorm(m_mecor$corfit$boot$coef[, '(Intercept)'])
qqline(m_mecor$corfit$boot$coef[, '(Intercept)'], col = 'red')
qqnorm(m_mecor$corfit$boot$coef[, 'cor_log_diffusion'])
qqline(m_mecor$corfit$boot$coef[, 'cor_log_diffusion'], col = 'red')
qqnorm(m_mecor$corfit$boot$coef[, 'dt_h'])
qqline(m_mecor$corfit$boot$coef[, 'dt_h'], col = 'red')
layout(1)

# plot the model (can't predict for new RE with the lme4 model)
coefs <- tibble(
  model = c('LM', 'lme4 HLM', 'mecor LM'),
  b_0 = c(summary(m_mecor)$uc$ci['(Intercept)', 'Estimate'],
          summary(m_lme4)$coefficients['(Intercept)', 'Estimate'],
          summary(m_mecor)$c$coefficients['(Intercept)', 'Estimate']),
  b_dt_h = c(summary(m_mecor)$uc$ci['dt_h', 'Estimate'],
             summary(m_lme4)$coefficients['dt_h', 'Estimate'],
             summary(m_mecor)$c$coefficients['dt_h', 'Estimate']),
  b_log_diffusion =
    c(summary(m_mecor)$uc$ci['log_diffusion', 'Estimate'],
      summary(m_lme4)$coefficients['log_diffusion', 'Estimate'],
      summary(m_mecor)$c$coefficients['cor_log_diffusion', 'Estimate']),
  se_log_diffusion = c(
    summary(m_mecor)$uc$coefficients['log_diffusion', 'Std. Error'],
    summary(m_lme4)$coefficients['log_diffusion', 'Std. Error'],
    summary(m_mecor)$c$coefficients['cor_log_diffusion', 'SE (btstr)']))

newd <- tibble(dt_h = 1, log_diffusion = seq(-1.3, 4.1, length.out = 400))

preds <- mutate(coefs, data = list(newd)) %>%
  unnest(data) %>%
  mutate(mu = b_0 + b_dt_h + b_log_diffusion * log_diffusion,
         lwr = mu - 2 * se_log_diffusion,
         upr = mu + 2 * se_log_diffusion)

# plot the estimated model effects
ggplot(mutate(d, dt_h = factor(dt_h))) +
  # vertical bars for 2 SD
  geom_errorbar(aes(x = log_diffusion,
                    ymin = log_speed - 2 * sqrt(var_log_speed),
                    ymax = log_speed + 2 * sqrt(var_log_speed)),
                width = 0, alpha = 0.05) +
  # horizontal bars for 2 SD
  geom_errorbarh(aes(xmin = log_diffusion - 2 * sqrt(var_log_diffusion),
                     xmax = log_diffusion + 2 * sqrt(var_log_diffusion),
                     y = log_speed),
                 height = 0, alpha = 0.05) +
  # data points
  geom_point(aes(log_diffusion, log_speed), alpha = 0.3) +
  # model estimates
  geom_ribbon(aes(log_diffusion, ymin = lwr, ymax = upr, fill = model),
              preds, alpha = 0.3) +
  geom_line(aes(log_diffusion, mu, color = model), preds) +
  scale_color_bright(name = 'Model') +
  scale_fill_bright(name = 'Model') +
  labs(x = 'log(diffusion)', y = 'log(speed)') +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.1, 0.9))

ggsave('figures/speed-diffusion-regression.png', height = 6, width = 6,
       dpi = 600, bg = 'white')

# looking for trends in the residuals
d <- mutate(
  d,
  mu_hat =
    summary(m_mecor)$c$coefficients['(Intercept)', 'Estimate'] +
    summary(m_mecor)$c$coefficients['dt_h', 'Estimate'] * dt_h +
    log_diffusion * b1_mecor,
  e = log_speed - mu_hat)
plot(e ~ sqrt(dof_speed), d)
plot(e ~ sqrt(dof_diff), d)
plot(e ~ Age, d)
plot(e ~ var_log_speed, d)
plot(e ~ var_log_diffusion, d)
plot(e ~ sex_treatment, d)
plot(e ~ animal, d)
