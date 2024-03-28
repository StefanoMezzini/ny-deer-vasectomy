library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('mgcv')      # for modeling
library('lubridate') # for working with dates
library('ggplot2')    # for fancy plots
library('gratia')    # for predicting from models
theme_set(theme_bw() + theme(text = element_text(face = 'bold')))
source('functions/gammals-variance-simulation-cis.R')
source('analysis/ref_dates.R')

d <- readRDS('data/years-1-and-2-data.rds') %>%
  mutate(animal = factor(animal),
         study_year = factor(study_year),
         sex_treatment = factor(paste(sex, study_site)),
         s_t_y = factor(paste(sex, study_site, study_year))) %>%
  # drop first 10 days for each animal to remove odd behaviors
  group_by(animal, study_year) %>%
  filter(date >= min(date) + 10) %>%
  ungroup() %>%
  mutate(dof_area = map_dbl(model, \(.m) summary(.m)$DOF['area']))

quantile(d$dof_area, probs = c(0, 0.1, 0.2, 0.25))

# data is not continuous throughout the year
hist(yday(d$date))

# create a column of days since August 1st
d <- mutate(d,
            days_since_aug_1 = if_else(
              # if in Aug, Sept, Oct, Nov, or Dec
              month(date) >= 8,
              # then calculate days since Aug 1
              date - as.Date(paste0(year(date), '-08-01')),
              # otherwise calculate days since Aug 1 of previous year 
              date - as.Date(paste0(year(date) - 1, '-08-01'))) %>%
              as.numeric()) # convert difftime to numeric

# plot the data ----
# plot overall raw data
ggplot(d, aes(days_since_aug_1, hr_est_95, group = animal)) +
  facet_wrap(~ sex_treatment + study_year, ncol = 2) +
  # coord_cartesian(ylim = c(0, 15)) +
  geom_line()

# using free scales to view cycles more clearly
p_hr <-
  mutate(d,
         sex = if_else(sex == 'f', 'females', 'males'),
         treatment = if_else(study_site == 'rockefeller', 'Rockefeller',
                             'Staten Island'),
         t_s = paste(treatment, sex),
         study_year = paste('Year', study_year)) %>%
  ggplot(aes(date, hr_est_95)) +
  facet_grid(t_s ~ study_year, scales = 'free') +
  geom_vline(xintercept = REF_DATES, col = 'red') +
  geom_line(aes(group = animal)) +
  labs(x = NULL,
       y = expression(bold('Estimated 7-day space-use requirements'~(km^2))))
p_hr

ggsave('figures/hr-estimates.png', p_hr, width = 8, height = 8, dpi = 600,
       bg = 'white')

# find males with large HRs on S Island
filter(d, sex_treatment == 'm staten_island', hr_est_95 > 15) %>%
  select(study_year, animal, Age)

# find males with large HRs on Rockefeller
filter(d, sex_treatment == 'm rockefeller', hr_est_95 > 7.5) %>%
  select(study_year, animal, Age)

# fit a Hierarchical Generalized Additive Model for Location and Scale ----
# not using cyclic cubic splines because the gap is too big
# for year 1 the gap is even bigger
range(d$days_since_aug_1) # not close to 0 to 365
365 - diff(range(d$days_since_aug_1))

# some very high HRs
quantile(d$hr_est_95, c(0.95, 0.97, 0.98, 0.99, 1))
# dropping HRs > 7.5 barely drops < 2% of the data
round(mean(d$hr_est_95 > 7.5), 3)

# location-scale model ----
if(FALSE) {
  m_hr <- gam(formula = list(
    # linear predictor for the mean
    hr_est_95 ~
      # temporal sex- and treatment-level trends with different smoothness
      #' using different smoothness for each `sex_treatment` and high `k`
      #' because females have cyclical estrous periods, while males do not
      s(days_since_aug_1, by = sex_treatment, k = 30, bs = 'ad') +
      # accounts for deviation from average between years
      #' keeping `by = sex_treatment` and high `k` to account for full
      #' differences between years
      s(days_since_aug_1, by = sex_treatment, study_year, k = 30, bs = 'sz') +
      # accounts for differences between individuals
      s(animal, bs = 're'),
    
    # linear predictor for the scale (sigma2 = mu^2 * scale)
    # allows mean-variance relationship to be different between sexes
    # sex- and treatment-level trends over season
    ~ s(days_since_aug_1, sex_treatment, k = 6, bs = 'fs')),
    
    family = gammals(),
    data = d,
    subset = hr_est_95 < 7.5, # dropping extreme HR estimates
    method = 'REML',
    control = gam.control(trace = TRUE))
  
  beepr::beep()
  appraise(m_hr, point_alpha = 0.1) # overdispersed residuals
  
  saveRDS(m_hr, paste0('models/m_hr-hgamls-', Sys.Date(), '.rds'))
} else {
  m_hr <- readRDS('models/m_hr-hgamls-2024-03-27.rds')
}

plot(m_hr, pages = 1, scheme = c(rep(1, 4), rep(0, 6)))
summary(m_hr)

ggplot() +
  geom_point(aes(m_hr$fitted.values[, 1], m_hr$model$hr_est_95),
             alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = 'red')

# check what groups cause oddly large outliers
d %>%
  filter(hr_est_95 < 7.5) %>%
  filter(hr_est_95 > 4, m_hr$fitted.values[, 1] < 4) %>%
  group_by(sex_treatment) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(prop = round(n / sum(n), 2))

d %>%
  filter(hr_est_95 < 7.5) %>%
  mutate(sex = if_else(sex == 'f', 'Females', 'Males'),
       treatment = if_else(study_site == 'rockefeller', 'Rockefeller',
                           'Staten Island'),
       fitted = m_hr$fitted.values[, 1],
       large = resid(m_hr) > 3) %>%
  ggplot() +
  facet_grid(treatment ~ sex) +
  geom_point(aes(fitted, hr_est_95, color = large, alpha = large)) +
  labs(x = 'Fitted values', y = 'Observed values') +
  scale_color_manual('Deviance residuals > 3', values = 1:2) +
  scale_alpha_manual('Deviance residuals > 3', values = c(0.3, 1)) +
  theme(legend.position = 'top')

ggsave('figures/hr-model-obs-fitted.png', width = 8, height = 8, dpi = 600,
       bg = 'white')

# check periods of oddly large outliers
d %>%
  filter(hr_est_95 < 7.5) %>%
  filter(hr_est_95 > 4, m_hr$fitted.values[, 1] < 4) %>%
  ggplot() +
  facet_grid(study_site ~ sex) +
  geom_histogram(aes(days_since_aug_1), color = 'black', fill = 'grey',
                 bins = 6)

# plot the estimated trends common between the two years ----
newd <-
  expand.grid(date = seq(as.Date('2021-09-01'),
                         as.Date('2022-05-30'),
                         length.out = 400),
              sex_treatment = unique(d$sex_treatment),
              study_year = 'new year',
              animal = 'new animal',
              s_t_y = 'new s_t_y') %>%
  mutate(days_since_aug_1 = if_else(
    # if in Aug, Sept, Oct, Nov, or Dec
    month(date) >= 8,
    # then calculate days since Aug 1
    date - as.Date(paste0(year(date), '-08-01')),
    # otherwise calculate days since Aug 1 of previous year 
    date - as.Date(paste0(year(date) - 1, '-08-01'))) %>%
      as.numeric(),
    sex = substr(sex_treatment, 1, 1),
    site = substr(sex_treatment, 3,
                  nchar(as.character(sex_treatment))))

preds_mu <- gammals_mean(model = m_hr, data = newd, nsims = 1e4,
                         unconditional = FALSE,
                         # not excluding the term results in NA
                         exclude = c('s(days_since_aug_1,study_year):sex_treatmentf rockefeller',
                                     's(days_since_aug_1,study_year):sex_treatmentf staten_island',
                                     's(days_since_aug_1,study_year):sex_treatmentm rockefeller',
                                     's(days_since_aug_1,study_year):sex_treatmentm staten_island')) %>%
  group_by(date, sex_treatment, days_since_aug_1, sex, site) %>%
  summarize(lwr_95 = quantile(mean, 0.025),
            mu = quantile(mean, 0.500),
            upr_95 = quantile(mean, 0.975),
            .groups = 'drop') %>%
  mutate(sex = case_when(sex == 'f' ~ 'Female',
                         sex == 'm' ~ 'Male'),
         site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                          site == 'staten_island' ~ 'Staten Island'))

preds_s <- gammals_var(model = m_hr, data = newd, nsims = 1e4,
                        unconditional = FALSE,
                        # not excluding the term results in NA
                        exclude = c('s(days_since_aug_1,study_year):sex_treatmentf rockefeller',
                                    's(days_since_aug_1,study_year):sex_treatmentf staten_island',
                                    's(days_since_aug_1,study_year):sex_treatmentm rockefeller',
                                    's(days_since_aug_1,study_year):sex_treatmentm staten_island')) %>%
  group_by(date, sex_treatment, days_since_aug_1, sex, site) %>%
  summarize(lwr_95 = quantile(variance, 0.025),
            s2 = quantile(variance, 0.500),
            upr_95 = quantile(variance, 0.975),
            .groups = 'drop') %>%
  mutate(sex = case_when(sex == 'f' ~ 'Female',
                         sex == 'm' ~ 'Male'),
         site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                          site == 'staten_island' ~ 'Staten Island'))

# mean HR
DATES <- as.Date(c('2021-09-15', '2021-11-15', '2022-01-15', '2022-03-15',
                   '2022-05-15'))
LABS <- format(DATES, '%B 15')
p_mu <-
  ggplot(preds_mu, aes(group = sex_treatment)) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab(expression(bold('Mean 7-day space-use requirements'~(km^2)))) +
  theme(legend.position = 'top'); p_mu

ggsave('figures/hr-mean.png', p_mu, width = 8, height = 8, dpi = 600,
       bg = 'white')

# standard deviation in HR ---
p_s <-
  ggplot(preds_s) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = sqrt(lwr_95), ymax = sqrt(upr_95),
                  fill = site), alpha = 0.3) +
  geom_line(aes(date, sqrt(s2), color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab(expression(bold('SD in 7-day space-use requirements'~(km^2)))) +
  theme(legend.position = 'top'); p_s

ggsave('figures/hr-sd.png', p_s, width = 8, height = 8, dpi = 600,
       bg = 'white')

# plot the estimated trends for each year ----
newd_years <-
  expand.grid(date = seq(as.Date('2021-09-01'),
                         as.Date('2022-05-30'),
                         length.out = 400),
              sex_treatment = unique(d$sex_treatment),
              study_year = 1:2,
              animal = 'new animal') %>%
  mutate(s_t_y = paste(sex_treatment, study_year),
         days_since_aug_1 = if_else(
           # if in Aug, Sept, Oct, Nov, or Dec
           month(date) >= 8,
           # then calculate days since Aug 1
           date - as.Date(paste0(year(date), '-08-01')),
           # otherwise calculate days since Aug 1 of previous year 
           date - as.Date(paste0(year(date) - 1, '-08-01'))) %>%
           as.numeric(),
         sex = substr(sex_treatment, 1, 1),
         site = substr(sex_treatment, 3,
                       nchar(as.character(sex_treatment))))

preds_mu_years <-
  gammals_mean(model = m_hr, data = newd_years, nsims = 1e4,
               unconditional = FALSE) %>%
  group_by(date, sex_treatment, days_since_aug_1, sex, site, study_year)%>%
  summarize(lwr_95 = quantile(mean, 0.025),
            mu = quantile(mean, 0.500),
            upr_95 = quantile(mean, 0.975),
            .groups = 'drop') %>%
  mutate(sex = case_when(sex == 'f' ~ 'Female',
                         sex == 'm' ~ 'Male'),
         site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                          site == 'staten_island' ~ 'Staten Island'))

preds_s_years <-
  gammals_var(model = m_hr, data = newd_years, nsims = 1e4,
              unconditional = FALSE) %>%
  group_by(date, sex_treatment, days_since_aug_1, sex, site, study_year)%>%
  summarize(lwr_95 = quantile(variance, 0.025),
            s2 = quantile(variance, 0.500),
            upr_95 = quantile(variance, 0.975),
            .groups = 'drop') %>%
  mutate(sex = case_when(sex == 'f' ~ 'Female',
                         sex == 'm' ~ 'Male'),
         site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                          site == 'staten_island' ~ 'Staten Island'))

# mean HR
p_mu_y <-
  ggplot(preds_mu_years) +
  facet_grid(sex ~ paste('Year', study_year)) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab(expression(bold('Mean 7-day space-use requirements'~(km^2)))) +
  theme(legend.position = 'top'); p_mu_y

ggsave('figures/hr-mean-years.png',
       p_mu_y, width = 16, height = 8, dpi = 600, bg = 'white')

# variance in HR ---
p_s <-
  ggplot(preds_s_years, aes(group = sex_treatment)) +
  facet_grid(sex ~ paste('Year', study_year)) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = sqrt(lwr_95), ymax = sqrt(upr_95),
                  fill = site), alpha = 0.3) +
  geom_line(aes(date, sqrt(s2), color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab(expression(bold('SD in 7-day space-use requirements'~(km^2)))) +
  theme(legend.position = 'top'); p_s

ggsave('figures/hr-sd-years.png', p_s, width = 16, height = 8, dpi = 600,
       bg = 'white')
