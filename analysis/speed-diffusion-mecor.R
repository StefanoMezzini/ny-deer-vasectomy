library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('ggplot2')   # for fancy plots
library('lme4')      # for models with random effects but no error calibration
library('mecor')     # for models with measurement error correction
library('ctmm')      # for movement models
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
           log_diffusion, var_log_diffusion, dof_diff)
  saveRDS(d, 'data/log-speed-diffusion.rds')
}
d

ggplot(d) +
  # vertical bars for 1 SD
  geom_errorbar(aes(x = log_speed,
                    ymin = log_diffusion - sqrt(var_log_diffusion),
                    ymax = log_diffusion + sqrt(var_log_diffusion)),
                width = 0, alpha = 0.05) +
  # horizontal bars for 1 SD
  geom_errorbarh(aes(xmin = log_speed - sqrt(var_log_speed),
                     xmax = log_speed + sqrt(var_log_speed),
                     y = log_diffusion),
                 height = 0, alpha = 0.05) +
  # data points
  geom_point(aes(log_speed, log_diffusion), alpha = 0.3) +
  # animal-level models
  geom_smooth(aes(log_speed, log_diffusion, group = animal),
              method = 'lm', se = FALSE, color = '#6495EDAA',
              formula = y ~ x) +
  labs(x = 'log(speed)', y = 'log(diffusion)')

# fit the naive model with lme4 ----
m_lm <- lm(log_diffusion ~ log_speed, data = d)

b0_lm <- coef(m_lm)['(Intercept)']
b1_lm <- coef(m_lm)['log_speed']

plot(log_diffusion ~ log_speed, d)
abline(a = b0_lm, b = b1_lm, col = 'forestgreen', lwd = 2)

# fit the naive model with lme4 ----
m_lme4 <- lmer(log_diffusion ~ 1 + log_speed +
                  (1 | animal) + (log_speed | animal),
                data = d, REML = TRUE)

b0_lme4 <- mean(coef(m_lme4)$animal[ , '(Intercept)'])
b1_lme4 <- mean(coef(m_lme4)$animal[ , 'log_speed'])

plot(log_diffusion ~ log_speed, d)
abline(a = b0_lm, b = b1_lm, col = 'forestgreen', lwd = 2)
abline(a = b0_lme4, b = b1_lme4, col = 'cornflowerblue', lwd = 2)

# fit the error-corrected model ----
m_mecor <- mecor(log_diffusion ~
                   MeasErrorRandom(substitute = log_speed,
                                   variance = mean(d$var_log_speed)),
                 data = d,
                 method = "standard",
                 B = 1e4)

# histograms and qqplots of bootstrap outputs for each coefficient
layout(1:2)
hist(m_mecor$corfit$boot$coef[, '(Intercept)'], breaks = 100)
hist(m_mecor$corfit$boot$coef[, 'cor_log_speed'], breaks = 100)

qqnorm(m_mecor$corfit$boot$coef[, '(Intercept)'])
qqline(m_mecor$corfit$boot$coef[, '(Intercept)'], col = 'red')
qqnorm(m_mecor$corfit$boot$coef[, 'cor_log_speed'])
qqline(m_mecor$corfit$boot$coef[, 'cor_log_speed'], col = 'red')
layout(1)

# plot the model
coefs <- summary(m_mecor)$c$coefficients
b0_mecor <- coefs['(Intercept)', 'Estimate']
b1_mecor <- coefs['cor_log_speed', 'Estimate']

# should include random intercepts, but can't
plot(log_diffusion ~ log_speed, d)
abline(a = b0_lm, b = b1_lm, col = 'forestgreen', lwd = 2)
abline(a = b0_lme4, b = b1_lme4, col = 'cornflowerblue', lwd = 2)
abline(a = b0_mecor, b = b1_mecor, col = 'darkorange', lwd = 2)

# looking for trends in the residuals
d <- mutate(d,
            mu_hat = b0_mecor + log_speed * b1_mecor,
            e = log_diffusion - mu_hat)
plot(e ~ sqrt(dof_speed), d)
plot(e ~ sqrt(dof_diff), d)
plot(e ~ Age, d)
plot(e ~ var_log_diffusion, d)
plot(e ~ var_log_speed, d)
plot(e ~ sex_treatment, d)
plot(e ~ animal, d)
abline(h = 1, lwd = 2, col = 'red')

d %>%
  group_by(animal) %>%
  summarize(sex_treatment = unique(sex_treatment),
            study_year = mean(as.numeric(study_year)),
            animal = unique(animal),
            mean_resid = mean(e),
            n_days = diff(range(date)),
            Age = mean(Age)) %>%
  filter(mean_resid > 1)
