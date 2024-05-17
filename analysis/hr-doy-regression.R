library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('mgcv')      # for modeling
library('lubridate') # for working with dates
library('ggplot2')    # for fancy plots
library('gratia')    # for predicting from models
source('functions/gammals-variance-simulation-cis.R')
source('analysis/ref_dates.R')
source('analysis/figures/default-theme.R')

d <- readRDS('data/years-1-and-2-data-no-akde.rds')

# figure of HRs and ESS before filtering
d %>%
  mutate(sex = if_else(sex == 'f', 'Female', 'Male'),
         study_site = if_else(study_site == 'rockefeller',
                              'Rockefeller', 'Staten Island')) %>%
  ggplot(aes(dof_area, hr_est_95)) +
  facet_grid(study_site ~ sex) +
  geom_hline(yintercept = 10, color = 'grey') +
  geom_point(alpha = 0.2) +
  xlab('Estimated number of range crossings') +
  scale_y_log10(expression(bold(paste('Estimated 7-day space use (k',
                                      m^2,', ', log[10], ' axis)'))),
                breaks = c(1e-4, 1e-2, 1, 100),
                labels = c('0.0001', '0.01', '1', '100'))

ggsave('figures/hr-ess-relationship.png', width = 12, height = 6,
       units = 'in', dpi = 600, bg = 'white')

# HR estimates with an ESS < 3 are likely biased
hist(d$dof_area, breaks = 100)
quantile(d$dof_area, probs = c(0, 0.01, 0.02, 0.05))
mean(d$dof_area < 3)

# HR estimates above 10 km^2 are unreliable for the study areas
hist(d$hr_est_95, breaks = 100)
quantile(d$hr_est_95, probs = c(0.95, 0.98, 0.99, 1))
mean(d$hr_est_95 > 10)

# check relationship between HR and ESS
plot(map_dbl(d$model, ctmm:::DOF.area), d$hr_est_95)
abline(v = 3, col = 'red') # ESS limit
abline(h = 10, col = 'red') # arbitrary HR limit

# plot the raw data ----
p_hr <-
  mutate(d,
         sex = if_else(sex == 'f', 'females', 'males'),
         treatment = if_else(study_site == 'rockefeller', 'Rockefeller',
                             'Staten Island'),
         t_s = paste(treatment, sex),
         study_year = paste('Year', study_year)) %>%
  ggplot(aes(date, hr_est_95)) +
  facet_grid(t_s ~ study_year, scales = 'free') +
  geom_hline(yintercept = 10, col = 'grey') +
  geom_vline(xintercept = REF_DATES, col = 'red') +
  geom_line(aes(group = animal, color = has_fawn)) +
  scale_color_manual('Confirmed fawn', values = c('black', 'darkorange2'),
                     labels = c('No', 'Yes')) +
  labs(x = NULL,
       y = expression(bold('Estimated 7-day space-use requirements'~
                             (km^2)))) +
  theme(legend.position = 'top')
p_hr

ggsave('figures/hr-estimates.png',
       p_hr, width = 8, height = 8, dpi = 600, bg = 'white')

# dropping excessively high HRs ----
# percent of data dropped by group when filtering to HR < 10
d %>%
  group_by(sex_treatment) %>%
  summarise(perc_outliers = round(mean(hr_est_95 > 10) * 100, 2))

# estimates should be unbiased for DOF ~>= 3, if movement is homogeneous
plot(hr_est_95 ~ dof_area, d)
abline(v = 3, col = 'red')

d <- filter(d, hr_est_95 < 10)

# estimates look more reasonable
plot(hr_est_95 ~ dof_area, d)

# fit a Hierarchical Generalized Additive Model for Location and Scale ----
# not using cyclic cubic splines because the gap is too big
# for year 1 the gap is even bigger
range(d$days_since_aug_1) # not close to 0 to 365
365 - diff(range(d$days_since_aug_1))

if(file.exists('models/m_hr-hgam-2024-05-12-poor-fit.rds')) {
  m_hr <- readRDS('models/m_hr-hgam-2024-05-12-poor-fit.rds')
} else {
  # fits in ~ 1 minute
  m_hr <- bam(
    hr_est_95 ~
      # by smooths require a separate explicit intercept for each group
      0 + sex_treatment +
      # average trends for each group
      s(days_since_aug_1, by = sex_treatment, k = 10, bs = 'tp') +
      # average year-level deviations from the average for each group
      s(days_since_aug_1, study_year, by = sex_treatment, k = 10, bs = 'sz') +
      # invidual- and year-level deviations from the average
      s(days_since_aug_1, animal_year, k = 10, bs = 'fs',
        xt = list(bs = 'cr')),
    family = Gamma(link = 'log'),
    data = d,
    method = 'fREML',
    discrete = TRUE,
    control = gam.control(trace = TRUE)); beepr::beep()
  
  saveRDS(m_hr, paste0('models/m_hr-hgam-', Sys.Date(), '-poor-fit.rds'))
}

#' *need a location-scale model*
appraise(m_hr, point_alpha = 0.3)
plot(m_hr, pages = 1, scheme = 0, scale = 0)
summary(m_hr, re.test = FALSE)

d <- mutate(d, e_hgam = resid(m_hr))

# differences in variance between groups
d %>%
  group_by(sex_treatment) %>%
  summarise(mean = mean(e_hgam),
            sd = sd(e_hgam))

ggplot(d) +
  facet_wrap(~ sex_treatment) +
  geom_histogram(aes(e_hgam), fill = 'grey80', color = 'black')

# over-dispersion in all groups
ggplot(d) +
  facet_wrap(~ sex_treatment) +
  geom_qq(aes(sample = e_hgam), alpha = 0.2) +
  geom_qq_line(aes(sample = e_hgam), color = 'red')

# the overdispersion occurs due to these large, sudden excursions
d %>%
  group_by(sex_treatment) %>%
  filter(animal_year == first(animal_year)) %>%
  ggplot() +
  facet_wrap(~ sex_treatment) +
  geom_line(aes(days_since_aug_1, e_hgam, group = animal_year))

# SI female trends would be different if we could account for excursions
mutate(d,
       excursion = e_hgam > 0.5) %>%
  ggplot(aes(days_since_aug_1, e_hgam)) +
  facet_wrap(~ sex_treatment + excursion) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 10))

if(FALSE) {
  # bs = 'tp', k = 15: dev.expl = 60.1%
  # bs = 'ad', k = 30: shows clear estrous cylcles (SIF), dev.expl = 68.2%
  # bs = 'ad', k = 30, s(age) and ti(age, doy) by group: dev.expl = 69.3%
  # bs = 'tp', k = 10, with REs for scale: dev.expl = 76.1%
  # bs = 'tp', k = 15, with REs for scale: dev.expl = 76.4%
  m_hr <- gam(formula = list(
    # linear predictor for the mean
    hr_est_95 ~
      # by smooths require a separate explicit intercept for each group
      0 + sex_treatment +
      # temporal sex- and treatment-level trends with different smoothness
      #' `k = 30` gives excessive wiggliness after estrous period
      s(days_since_aug_1, by = sex_treatment, k = 10, bs = 'tp') +
      # accounts for deviation from average between years
      s(days_since_aug_1, study_year, by = sex_treatment, k = 10, bs = 'sz') +
      s(animal_year, bs = 're'),
      # invidual- and year-level deviations from the average
      # s(days_since_aug_1, animal_year, k = 10, bs = 'fs',
        # xt = list(bs = 'cr')),
    
    # linear predictor for the scale (sigma2 = mu^2 * scale)
    # allows mean-variance relationship to be different between sexes
    # sex- and treatment-level trends over season
    ~ 0 + sex_treatment +
      s(days_since_aug_1, by = sex_treatment, k = 5, bs = 'tp') +
      s(animal_year, bs = 're')),
    
    family = gammals(),
    data = d,
    method = 'REML',
    control = gam.control(trace = TRUE)); beepr::beep()
  
  saveRDS(m_hr, paste0('models/m_hr-hgamls-', Sys.Date(), '.rds'))
} else {
  m_hr <- readRDS('models/m_hr-hgamls-2024-05-12.rds')
}

# residuals are still overdispersed
appraise(m_hr, point_alpha = 0.1, n_bins = 30)
plot(m_hr, pages = 1, scheme = 0)
summary(m_hr, re.test = FALSE)

# there are trends in the SI females' residuals
d <- mutate(d,
            sex_treatment = as.character(sex_treatment),
            site = if_else(substr(sex_treatment, 3, 100) == 'staten_island',
                           'Staten Island',
                           'Rockefeller State Park Preserve'),
            e = resid(m_hr)) %>%
  as_tibble()

cowplot::plot_grid(
  ggplot(d, aes(days_since_aug_1, e)) +
    facet_grid(sex ~ site) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = 'gam', formula = y ~ s(x)) +
    labs(x = expression(bold(Days~since~August~1^'st')),
         y = 'Deviance residuals'),
  ggplot(d, aes(days_since_aug_1, e^2)) +
    facet_grid(sex ~ site) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = 'gam', formula = y ~ s(x)) +
    labs(x = expression(bold(Days~since~August~1^'st')),
         y = 'Squared deviance residuals'))

# check what groups cause oddly large outliers
d %>%
  mutate(sex = if_else(sex == 'f', 'Females', 'Males'),
         treatment = if_else(study_site == 'rockefeller', 'Rockefeller',
                             'Staten Island'),
         fitted = m_hr$fitted.values[, 1],
         large = resid(m_hr) > 3,
         sex_fawn = paste0(
           sex,
           case_when(sex == 'Males' ~ '',
                     sex == 'Females' & has_fawn ~ ' with known fawn',
                     sex == 'Females' & ! has_fawn ~ ' with unknown fawn'))) %>%
  filter(dof_area > 10) %>%
  ggplot() +
  facet_grid(treatment ~ sex_fawn) +
  geom_point(aes(fitted, hr_est_95, color = large, alpha = large)) +
  geom_abline(slope = 1, intercept = 0, color = 'grey') +
  labs(x = 'Fitted values', y = 'Observed values') +
  xlim(c(0, max(m_hr$model$hr_est_95))) +
  ylim(c(0, max(m_hr$model$hr_est_95))) +
  scale_color_manual('Deviance residuals > 3', values = 1:2) +
  scale_alpha_manual('Deviance residuals > 3', values = c(0.2, 1)) +
  theme(legend.position = c(0.1, 0.93))

ggsave('figures/hr-model-obs-fitted.png', width = 12, height = 8,
       dpi = 600, bg = 'white')

# check periods of oddly large outliers
d %>%
  filter(hr_est_95 > 4, m_hr$fitted.values[, 1] < 4) %>%
  ggplot() +
  facet_grid(study_site ~ sex) +
  geom_histogram(aes(days_since_aug_1), color = 'black', fill = 'grey',
                 bins = 6)

# plot the estimated common trends between the two years ----
newd <-
  expand.grid(date = seq(as.Date('2021-09-01'),
                         as.Date('2022-05-30'),
                         length.out = 400),
              sex_treatment = unique(d$sex_treatment),
              study_year = 'new year',
              animal_year = 'new animal',
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
              alpha = 0.5) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab(expression(bold('Mean space-use requirements'~(km^2)))) +
  theme(legend.position = 'top'); p_mu

ggsave('figures/hr-mean.png',
       p_mu, width = 8, height = 8, dpi = 600, bg = 'white')

# standard deviation in HR ---
p_s <-
  ggplot(preds_s) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = sqrt(lwr_95), ymax = sqrt(upr_95),
                  fill = site), alpha = 0.5) +
  geom_line(aes(date, sqrt(s2), color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab(expression(bold('SD in space-use requirements'~(km^2)))) +
  theme(legend.position = 'top'); p_s

ggsave('figures/hr-sd.png',
       p_s, width = 8, height = 8, dpi = 600, bg = 'white')

# plot the estimated trends for each year ----
newd_years <-
  expand.grid(date = seq(as.Date('2021-09-01'),
                         as.Date('2022-05-30'),
                         length.out = 400),
              sex_treatment = unique(d$sex_treatment),
              study_year = 1:2,
              animal_year = 'new animal') %>%
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
              alpha = 0.5) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab(expression(bold('Mean space-use requirements'~(km^2)))) +
  theme(legend.position = 'top'); p_mu_y

ggsave('figures/hr-mean-years.png',
       p_mu_y, width = 16, height = 8, dpi = 600, bg = 'white')

# variance in HR ---
p_s_y <-
  ggplot(preds_s_years, aes(group = sex_treatment)) +
  facet_grid(sex ~ paste('Year', study_year)) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = sqrt(lwr_95), ymax = sqrt(upr_95),
                  fill = site), alpha = 0.5) +
  geom_line(aes(date, sqrt(s2), color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab(expression(bold('SD in space-use requirements'~(km^2)))) +
  theme(legend.position = 'top'); p_s_y

ggsave('figures/hr-sd-years.png',
       p_s_y, width = 16, height = 8, dpi = 600, bg = 'white')
