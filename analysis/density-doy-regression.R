library('ctmm')      # for movement modeling
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
library('mgcv')      # for modeling
library('gratia')    # for ggplot-based diagnostics
source('analysis/ref_dates.R')
source('analysis/figures/default-theme.R')
source('functions/betals.r')
source('functions/betals-variance-sims.R')

d <-
  # find density values using the models fit to the full telemetry
  readRDS('models/full-telemetry-movement-models-2024-04-20.rds') %>%
  transmute(animal,
            hr_est_95,
            tel = map(tel, \(.t) {
              data.frame(.t) %>%
                mutate(date = as.Date(as.POSIXct(timestamp))) %>%
                # calculate mean daily density
                group_by(date) %>%
                summarise(density = mean(density), .groups = 'drop') %>%
                # drop first 7 days
                filter(date > min(date) + 7)
            }
            ),
            study_year = factor(study.year)) %>%
  # add moving window information
  left_join(readRDS('data/years-1-and-2-data-no-akde.rds') %>%
              group_by(animal, study_year) %>%
              slice(1) %>%
              select(animal, study_year, study_year2, study_site, sex, Age,
                     sex_treatment, has_fawn) %>%
              ungroup(),
            by = c('animal', 'study_year')) %>%
  mutate(animal = factor(animal),
         animal_year = factor(paste(animal, study_year))) %>%
  unnest(tel) %>%
  as_tibble() %>%
  mutate(days_since_aug_1 =
           if_else(
             # if in Aug, Sept, Oct, Nov, or Dec
             month(date) >= 8,
             # then calculate days since Aug 1
             date - as.Date(paste0(year(date), '-08-01')),
             # otherwise calculate days since Aug 1 of previous year 
             date - as.Date(paste0(year(date) - 1, '-08-01')))%>%
           as.numeric()) # convert difftime to numeric

# checks
if(FALSE) {
  d
  
  # ensure years were added correctly (should be an empty tibble)
  d %>%
    mutate(date_year = year(date)) %>%
    filter(study_year == 1 & date_year == 2023) %>%
    group_by(animal_year) %>%
    slice(1)
  
  # density is much more constrained than the HR estimates
  layout(t(1:2))
  hist(d$density, breaks = 50)
  hist(d$hr_est_95, breaks = 50)
  layout(1)
  quantile(d$hr_est_95, probs = c(0.5, 0.95, 0.97, 0.98, 0.99, 1))
  
  ggplot(d) +
    facet_grid(sex ~ study_site) +
    geom_smooth(aes(days_since_aug_1, density))
  
  ggplot(d) +
    facet_grid(sex ~ study_site + study_year) +
    geom_line(aes(days_since_aug_1, density, group = animal_year),
              alpha = 0.1)
}

# plot the data
# see SI females in year 2
p_density <-
  mutate(d,
         sex = if_else(sex == 'f', 'females', 'males'),
         treatment = if_else(study_site == 'rockefeller', 'Rockefeller',
                             'Staten Island'),
         t_s = paste(treatment, sex),
         study_year = paste('Year', study_year)) %>%
  ggplot(aes(date, density)) +
  facet_grid(t_s ~ study_year, scales = 'free') +
  geom_vline(xintercept = REF_DATES, col = 'red') +
  geom_line(aes(group = animal, color = has_fawn), alpha = 0.3) +
  scale_color_manual('Confirmed fawn', values = c('black', 'darkorange2'),
                     labels = c('No', 'Yes')) +
  labs(x = NULL,
       y = expression(bold('Estimated 7-day space-use requirements'~
                             (km^2)))) +
  theme(legend.position = 'top')
p_density

ggsave('figures/density-estimates.png',
       p_density, width = 8, height = 8, dpi = 600, bg = 'white')

# fit a Hierarchical Generalized Additive Model for Location and Scale ----
if(file.exists('models/m_density-hgamls.rds')) {
  m_density <- readRDS('models/m_density-hgamls.rds')
} else {
  m_density <- gam(formula = list(
    # linear predictor for the mean
    density ~
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
    
    # linear predictor for the scale (sigma2 = mu * (1-mu) * scale)
    # allows mean-variance relationship to be different between sexes
    # sex- and treatment-level trends over season
    ~ sex_treatment +
      s(days_since_aug_1, by = sex_treatment, k = 15, bs = 'tp') +
      s(days_since_aug_1, study_year, by = sex_treatment, k = 15, bs = 'sz') +
      s(days_since_aug_1, animal_year, k = 15, bs = 'fs',
        xt = list(bs = 'cr'))),
    
    family = betals(),
    data = d,
    method = 'REML',
    control = gam.control(trace = TRUE))
  saveRDS(m_density, paste0('models/m_density-hgamls.rds'))
}

appraise(m_density, point_alpha = 0.1, type = 'pearson')
plot(m_density, pages = 1, scheme = 0)
summary(m_density, re.test = FALSE) # good deviance explained

# plot the estimated common trends between the two years ----
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

if(file.exists('models/predictions/density-preds_mu.rds')) {
  preds_mu <- readRDS('models/predictions/density-preds_mu.rds')
} else {
  preds_mu <- betals_mean(model = m_density, data = newd, nsims = 1e4,
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
  saveRDS(preds_mu, 'models/predictions/density-preds_mu.rds')
}

if(file.exists('models/predictions/density-preds_s2.rds')) {
  preds_s2 <- readRDS('models/predictions/density-preds_s2.rds')
} else {
  preds_s2 <- betals_var(model = m_density, data = newd, nsims = 1e4,
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
  saveRDS(preds_s2, 'models/predictions/density-preds_s2.rds')
}

# mean density
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
  ylab('Mean daily excursivity') +
  theme(legend.position = 'top'); p_mu

ggsave('figures/density-mean.png',
       p_mu, width = 8, height = 8, dpi = 600, bg = 'white')

# standard deviation in density ----
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
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab('SD in daily excursivity') +
  theme(legend.position = 'top'); p_s

ggsave('figures/density-sd.png',
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

if(file.exists('models/predictions/density-preds_mu_years.rds')) {
  preds_mu <- readRDS('models/predictions/density-preds_mu_years.rds')
} else {
  preds_mu_years <-
    betals_mean(model = m_density, data = newd_years, nsims = 1e4,
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
  saveRDS(preds_mu_years, 'models/predictions/density-preds_mu_years.rds')
}

if(file.exists('models/predictions/density-preds_s2_years.rds')) {
  preds_s2_years <- readRDS('models/predictions/density-preds_s2_years.rds')
} else {
  preds_s_years <-
    betals_var(model = m_density, data = newd_years, nsims = 1e4,
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
  saveRDS(preds_s2_years, 'models/predictions/density-preds_s2_years.rds')
}

# mean density
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
  ylab('Mean daily excursivity') +
  theme(legend.position = 'top'); p_mu_y

ggsave('figures/density-mean-years.png',
       p_mu_y, width = 16, height = 8, dpi = 600, bg = 'white')

# variance in density ---
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
  ylab('SD in daily excursivity') +
  theme(legend.position = 'top'); p_s_y

ggsave('figures/density-sd-years.png',
       p_s_y, width = 16, height = 8, dpi = 600, bg = 'white')
