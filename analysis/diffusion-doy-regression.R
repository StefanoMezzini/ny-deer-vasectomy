library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('mgcv')      # for modeling
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
library('gratia')    # for predicting from models
library('ctmm')      # for movement models
source('functions/gammals-variance-simulation-cis.R')
source('analysis/ref_dates.R')
source('analysis/figures/default-theme.R')

d <- readRDS('data/years-1-and-2-data-no-akde.rds')

# most windows with a non-resolved diffusion have 1 sample per hour 
d_int <- d %>%
  filter(is.na(diffusion_est)) %>%
  transmute(interval = map(tel, \(.tel) {
    median(diff(.tel$timestamp), na.rm = TRUE)
  }),
  interval = map_dbl(interval, \(int) {
    units <- case_when(attr(int, 'units') == 'mins' ~ 'minutes',
                       attr(int, 'units') == 'secs' ~ 'seconds',
                       TRUE ~ attr(int, 'units'))
    
    'hours' %#% as.numeric(int) %#% units
  }))

ggplot(d_int, aes(interval)) + 
  geom_histogram(center = 0, bins = 50) +
  scale_x_continuous(trans = 'sqrt', breaks = c(0, 1, 10, 20, 30, 40, 50)) +
  scale_y_continuous(trans = 'sqrt') +
  labs(x = 'Sampling interval (hours; sqrt axis)', y = 'Count (sqrt axis)')

mean(d_int$interval > 2)
mean(d_int$interval > 1.05)

# drop rows with NA diffusion
sum(is.na(d$diffusion_est))
d <- filter(d, ! is.na(diffusion_est)) # 68 NA values

# some low values, but common between both males and females in Rockefeller
quantile(d$dof_diff)
ggplot(d, aes(dof_diff)) +
  facet_grid(sex ~ study_site) +
  geom_histogram()

# plot the raw data ----
p_diffusion <-
  mutate(d,
         sex = if_else(sex == 'f', 'females', 'males'),
         treatment = if_else(study_site == 'rockefeller', 'Rockefeller',
                             'Staten Island'),
         t_s = paste(treatment, sex),
         study_year = paste('Year', study_year)) %>%
  ggplot(aes(date, diffusion_est, group = animal)) +
  facet_grid(t_s ~ study_year, scales = 'free') +
  geom_vline(xintercept = REF_DATES, col = 'red') +
  geom_line() +
  labs(x = NULL, y = expression(bold('Estimated diffusion'~(km^2/day))))
p_diffusion

ggsave('figures/diffusion-estimates.png', p_diffusion, width = 8,
       height = 8, dpi = 600, bg = 'white')

# fit a Hierarchical Generalized Additive Model for Location and Scale ----
# not using cyclic cubic splines because the gap is too big
# for year 1 the gap is even bigger
range(d$days_since_aug_1) # not close to 0 to 365
365 - diff(range(d$days_since_aug_1))

if(file.exists('models/m_diffusion-hgamls.rds')) {
  m_diffusion <- readRDS('models/m_diffusion-hgamls.rds')
} else {
  # no clear temporal trends
  ggplot(d, aes(days_since_aug_1, diffusion_est, group = animal_year)) +
    facet_wrap(~ sex_treatment, scales = 'free') +
    geom_line(alpha = 0.3)
  
  # setup takes ~15-20 minutes
  m_diffusion <- gam(list(
    # linear predictor for the mean
    diffusion_est ~
      # by smooths require a separate explicit intercept for each group
      sex_treatment +
      # temporal sex- and treatment-level trends with different smoothness
      #' `k = 30` gives excessive wiggliness after estrous period
      s(days_since_aug_1, by = sex_treatment, k = 15, bs = 'tp') +
      # accounts for deviation from average between years
      s(days_since_aug_1, study_year, by = sex_treatment, k = 15, bs = 'sz') +
      # invidual- and year-level deviations from the average
      s(days_since_aug_1, animal_year, k = 15, bs = 'fs',
        xt = list(bs = 'cr')),
    
    # linear predictor for the scale (sigma2 = mu^2 * scale)
    # allows mean-variance relationship to be different between sexes
    # sex- and treatment-level trends over season
    ~ sex_treatment +
      s(days_since_aug_1, by = sex_treatment, k = 15, bs = 'tp') +
      s(days_since_aug_1, study_year, by = sex_treatment, k = 15, bs = 'sz') +
      s(days_since_aug_1, animal_year, k = 15, bs = 'fs',
        xt = list(bs = 'cr'))),
    
    family = gammals(),
    data = d,
    method = 'REML',
    control = gam.control(trace = TRUE))
  
  saveRDS(m_diffusion,
          paste0('models/m_diffusion-hgamls.rds'))
}

appraise(m_diffusion, point_alpha = 0.05, n_bins = 30)
plot(m_diffusion, pages = 1, scale = 0)
summary(m_diffusion, re.test = FALSE)

# check what groups cause oddly large outliers
mutate(d,
       sex = if_else(sex == 'f', 'Females', 'Males'),
       treatment = if_else(study_site == 'rockefeller', 'Rockefeller',
                           'Staten Island'),
       fitted = m_diffusion$fitted.values[, 1],
       large = resid(m_diffusion) > 3) %>%
  ggplot() +
  facet_grid(treatment ~ sex) +
  geom_point(aes(fitted, diffusion_est, color = large, alpha = large)) +
  geom_abline(slope = 1, intercept = 0, color = 'grey') +
  labs(x = 'Fitted values', y = 'Observed values') +
  scale_color_manual('Deviance residuals > 3', values = 1:2) +
  scale_alpha_manual('Deviance residuals > 3', values = c(0.3, 1)) +
  theme(legend.position = 'top')

ggsave('figures/diffusion-model-obs-fitted.png', width = 8, height = 8,
       dpi = 600, bg = 'white')

p_diffusion_2 <-
  mutate(d,
         sex = if_else(sex == 'f', 'females', 'males'),
         treatment = if_else(study_site == 'rockefeller', 'Rockefeller',
                             'Staten Island'),
         t_s = paste(treatment, sex),
         study_year = paste('Year', study_year),
         large = resid(m_diffusion) > 3) %>%
  ggplot(aes(date, diffusion_est, group = animal)) +
  facet_grid(t_s ~ study_year, scales = 'free') +
  geom_vline(xintercept = REF_DATES, col = 'red') +
  geom_line(aes(color = large)) +
  labs(x = NULL, y = expression(bold('Estimated diffusion'~(km^2/day)))) +
  scale_color_manual('Deviance residuals > 3',
                     values = c('black', 'darkorange2')) +
  theme(legend.position = 'top')
p_diffusion_2

ggsave('figures/diffusion-estimates-outliers.png', p_diffusion_2,
       width = 8, height = 8, dpi = 600, bg = 'white')

# check periods of oddly large outliers
d %>%
  filter(resid(m_diffusion) > 3) %>%
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
              animal_year = 'new animal') %>%
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

if(file.exists('models/predictions/diffusion-preds_mu.rds')) {
  preds_mu <- readRDS('models/predictions/diffusion-preds_mu.rds')
} else {
  preds_mu <- gammals_mean(model = m_diffusion, data = newd, nsims = 1e4,
                           unconditional = FALSE,
                           # not excluding the term results in NA
                           exclude =
                             c(paste0('s(days_since_aug_1,study_year):sex_treatment',
                                      unique(d$sex_treatment)),
                               's(days_since_aug_1,animal_year)',
                               paste0('s.1(days_since_aug_1,study_year):sex_treatment',
                                      unique(d$sex_treatment)),
                               's.1(days_since_aug_1,animal_year)')) %>%
    group_by(date, sex_treatment, days_since_aug_1, sex, site) %>%
    summarize(lwr_95 = quantile(mean, 0.025),
              mu = quantile(mean, 0.500),
              upr_95 = quantile(mean, 0.975),
              .groups = 'drop') %>%
    mutate(sex = case_when(sex == 'f' ~ 'Female',
                           sex == 'm' ~ 'Male'),
           site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                            site == 'staten_island' ~ 'Staten Island'))
  saveRDS(preds_mu, 'models/predictions/diffusion-preds_mu.rds')
}

if(file.exists('models/predictions/diffusion-preds_s2.rds')) {
  preds_s2 <- readRDS('models/predictions/diffusion-preds_s2.rds')
} else {
  preds_s2 <- gammals_var(model = m_diffusion, data = newd, nsims = 1e4,
                          unconditional = FALSE,
                          # not excluding the term results in NA
                          exclude =
                            c(paste0('s(days_since_aug_1,study_year):sex_treatment',
                                     unique(d$sex_treatment)),
                              's(days_since_aug_1,animal_year)',
                              paste0('s.1(days_since_aug_1,study_year):sex_treatment',
                                     unique(d$sex_treatment)),
                              's.1(days_since_aug_1,animal_year)')) %>%
    group_by(date, sex_treatment, days_since_aug_1, sex, site) %>%
    summarize(lwr_95 = quantile(variance, 0.025),
              s2 = quantile(variance, 0.500),
              upr_95 = quantile(variance, 0.975),
              .groups = 'drop') %>%
    mutate(sex = case_when(sex == 'f' ~ 'Female',
                           sex == 'm' ~ 'Male'),
           site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                            site == 'staten_island' ~ 'Staten Island'))
  saveRDS(preds_s2, 'models/predictions/diffusion-preds_s2.rds')
}

# mean diffusion
DATES <- as.Date(c('2021-10-01', '2021-12-01', '2022-02-01', '2022-04-01'))
LABS <- format(DATES, '%B 1')
p_mu <-
  ggplot(preds_mu) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.5) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab(expression(bold('Mean diffusion'~(km^2/day)))) +
  theme(legend.position = 'top'); p_mu

ggsave('figures/diffusion-mean.png',
       p_mu, width = 8, height = 8, dpi = 600, bg = 'white')

# standard deviation in diffusion
p_s <-
  ggplot(preds_s2) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = sqrt(lwr_95), ymax = sqrt(upr_95),
                  fill = site), alpha = 0.5) +
  geom_line(aes(date, sqrt(s2), color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab(expression(bold('SD in diffusion'~(km^2/day)))) +
  theme(legend.position = 'top'); p_s

ggsave('figures/diffusion-sd.png', p_s, width = 8, height = 8, dpi = 600,
       bg = 'white')

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

if(file.exists('models/predictions/diffusion-preds_mu_years.rds')) {
  preds_mu_years <- readRDS('models/predictions/diffusion-preds_mu_years.rds')
} else {
  preds_mu_years <-
    gammals_mean(model = m_diffusion, data = newd_years, nsims = 1e4,
                 unconditional = FALSE,
                 exclude =
                   c('s(days_since_aug_1,animal_year)',
                     's.1(days_since_aug_1,animal_year)')) %>%
    group_by(date, sex_treatment, days_since_aug_1, sex, site, study_year)%>%
    summarize(lwr_95 = quantile(mean, 0.025),
              mu = quantile(mean, 0.500),
              upr_95 = quantile(mean, 0.975),
              .groups = 'drop') %>%
    mutate(sex = case_when(sex == 'f' ~ 'Female',
                           sex == 'm' ~ 'Male'),
           site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                            site == 'staten_island' ~ 'Staten Island'))
  saveRDS(preds_mu_years, 'models/predictions/diffusion-preds_mu_years.rds')
}

if(file.exists('models/predictions/diffusion-preds_s2_years.rds')) {
  preds_s2_years <- readRDS('models/predictions/diffusion-preds_s2_years.rds')
} else {
  preds_s2_years <-
    gammals_var(model = m_diffusion, data = newd_years, nsims = 1e4,
                unconditional = FALSE,
                exclude =
                  c('s(days_since_aug_1,animal_year)',
                    's.1(days_since_aug_1,animal_year)')) %>%
    group_by(date, sex_treatment, days_since_aug_1, sex, site, study_year)%>%
    summarize(lwr_95 = quantile(variance, 0.025),
              s2 = quantile(variance, 0.500),
              upr_95 = quantile(variance, 0.975),
              .groups = 'drop') %>%
    mutate(sex = case_when(sex == 'f' ~ 'Female',
                           sex == 'm' ~ 'Male'),
           site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                            site == 'staten_island' ~ 'Staten Island'))
  saveRDS(preds_s2_years, 'models/predictions/diffusion-preds_s2_years.rds')
}

# mean diffusion
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
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab(expression(bold('Mean diffusion'~(km^2/day)))) +
  theme(legend.position = 'top'); p_mu_y

ggsave('figures/diffusion-mean-years.png',
       p_mu_y, width = 16, height = 8, dpi = 600, bg = 'white')

# variance in diffusion ---
p_s_y <-
  ggplot(preds_s2_years) +
  facet_grid(sex ~ paste('Year', study_year)) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = sqrt(lwr_95), ymax = sqrt(upr_95),
                  fill = site), alpha = 0.5) +
  geom_line(aes(date, sqrt(s2), color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab(expression(bold('SD in diffusion'~(km^2/day)))) +
  theme(legend.position = 'top'); p_s_y

ggsave('figures/diffusion-sd-years.png', p_s_y, width = 16, height = 8,
       dpi = 600, bg = 'white')
