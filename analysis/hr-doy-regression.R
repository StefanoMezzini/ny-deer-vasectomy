library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('mgcv')      # for modeling
library('lubridate') # for working with dates
library('ggplot2')    # for fancy plots
library('gratia')    # for predicting from models
theme_set(theme_bw() + theme(text = element_text(face = 'bold')))
source('functions/gammals-variance-simulation-and-derivatives.R')

d <- readRDS('data/years-1-and-2-data.rds') %>%
  mutate(animal = factor(animal),
         study_year = factor(study_year),
         sex_treatment = factor(paste(sex, study_site)),
         s_t_y = factor(paste(sex, study_site, study_year))) %>%
  # drop first 10 days for each animal to remove odd behaviors
  group_by(animal, study_year) %>%
  filter(date >= min(date) + 10)

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
              as.numeric(), # convert difftime to numeric
            weight = 1 / (hr_upr_95 - hr_lwr_95)) %>% # account for HR uncertainty
  group_by(animal) %>%
  mutate(weight = weight / mean(weight)) %>% # normalize weights
  ungroup()

hist(d$days_since_aug_1) # now without breaks

# ensure weights are correct
d %>%
  group_by(animal) %>%
  summarize(weights = round(sum(weight)),
            count = n()) %>%
  mutate(check = weights == count) %>%
  pull(check) %>%
  all()

# plot the data ----
# plot overall raw data
ggplot(d, aes(days_since_aug_1, hr_est_95, group = animal)) +
  facet_wrap(~ sex_treatment) +
  # coord_cartesian(ylim = c(0, 15)) +
  geom_line()

# using free scales to view cycles more clearly
ggplot(d, aes(days_since_aug_1, hr_est_95)) +
  facet_wrap(~ sex_treatment, scales = 'free_y') +
  geom_vline(xintercept = 100, col = 'red') +
  geom_line(aes(group = animal))

# find odd male from S Island
filter(d, sex_treatment == 'm staten_island', hr_est_95 > 15) %>%
  select(study_year, animal, Age)

# find males with large HRs on Rockefeller
filter(d, sex_treatment == 'm rockefeller', hr_est_95 > 7.5) %>%
  select(study_year, animal, Age)

# start by fitting a Hierarchical Generalized Additive Model ----
# not using cyclic cubic splines because the gap is too big (66 days)
# for year 1 the gap is even bigger
range(d$days_since_aug_1) # not close to 0 to 365
365 - (317 - 18)

# location-scale model ----
if(FALSE) {
  m_hr <- gam(formula = list(
    #' not using fixed effects of `sex_treatment` because they are too
    #' uncertain due to the large HRs in May.
    #' the REs of `sex_treatment` are much more constrained.
    
    # linear predictor for the mean
    hr_est_95 ~
      # random intercept of sex and treatment
      s(sex_treatment, bs = 're') +
      # temporal sex- and treatment-level trends with different
      s(days_since_aug_1, sex_treatment, k = 10, bs = 'fs') +
      # accounts for differences in trends between years
      s(days_since_aug_1, s_t_y, k = 10, bs = 'fs') +
      # accounts for differences between individuals
      s(animal, bs = 're'),
    
    # linear predictor for the scale (sigma2 = mu^2 * scale)
    # allows mean-variance relationship to be different between sexes
    # sex- and treatment-level trends over season
    ~ s(days_since_aug_1, sex_treatment, k = 5, bs = 'fs')),
    
    family = gammals(),
    data = d,
    weights = weight,
    method = 'REML',
    control = gam.control(trace = TRUE))
  
  saveRDS(m_hr, paste0('models/m_hr-hgamls-', Sys.Date(), '.rds'))
} else {
  m_hr <- readRDS('models/m_hr-hgamls-2024-03-15.rds')
}

draw(m_hr, scales = 'fixed')

# check the model diagnostics ----
if(FALSE) {
  appraise(m_hr, point_alpha = 0.1) # overdispersed residuals
  
  # check trends in the residuals
  d$e <- residuals(m_hr)
  
  # no visible trends in the residuals over day of year
  ggplot(d, aes(days_since_aug_1, e)) +
    facet_grid(sex ~ study_site) +
    geom_point(alpha = 0.1) +
    geom_smooth()
  
  # check residuals across groups
  ggplot(d, aes(e)) +
    facet_grid(sex ~ study_site) +
    geom_vline(xintercept = 0) +
    geom_density(fill = 'grey', bw = 0.25) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.63))
}

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
                         exclude = 's(days_since_aug_1,s_t_y)') %>%
  group_by(date, sex_treatment, days_since_aug_1, sex, site) %>%
  summarize(lwr_95 = quantile(mean, 0.025),
            mu = quantile(mean, 0.500),
            upr_95 = quantile(mean, 0.975),
            .groups = 'drop') %>%
  mutate(sex = case_when(sex == 'f' ~ 'Female',
                             sex == 'm' ~ 'Male'),
         site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                              site == 'staten_island' ~ 'Staten Island'))

preds_s2 <- gammals_var(model = m_hr, data = newd, nsims = 1e4,
                        unconditional = FALSE,
                        # not excluding the term results in NA
                        exclude = 's(days_since_aug_1,s_t_y)') %>%
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
p_mu <-
  ggplot(preds_mu, aes(group = sex_treatment)) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = as.Date('2021-11-09'), col = 'red')+
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.2) +
  geom_line(aes(date, mu, color = site)) +
  scale_color_brewer('Site', type = 'qual', palette = 2,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  labs(x = NULL,
       y = expression(bold('Mean 7-day home range size'~(km^2)))) +
  theme(legend.position = 'top'); p_mu

ggsave('figures/hr-mean.png', p_mu, width = 8, height = 8, dpi = 600,
       bg = 'white')

# variance in HR ---
p_s <-
  ggplot(preds_s2) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = as.Date('2021-11-09'), col = 'red')+
  geom_ribbon(aes(date, ymin = sqrt(lwr_95), ymax = sqrt(upr_95),
                  fill = site), alpha = 0.3) +
  geom_line(aes(date, sqrt(s2), color = site)) +
  scale_color_brewer('Site', type = 'qual', palette = 2,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  labs(x = NULL,
       y = expression(bold('SD in 7-day home range size'~(km^2)))) +
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
         site_ = case_when(site == 'rockefeller' ~ 'Rockefeller',
                          site == 'staten_island' ~ 'Staten Island'))

preds_s2_years <-
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
p_mu <-
  ggplot(preds_mu_years) +
  facet_grid(sex ~ study_year) +
  geom_vline(xintercept = as.Date('2021-11-09'), col = 'red')+
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.2) +
  geom_line(aes(date, mu, color = site)) +
  scale_color_brewer('Site', type = 'qual', palette = 2,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  labs(x = NULL,
       y = expression(bold('Mean 7-day home range size'~(km^2)))) +
  theme(legend.position = 'top'); p_mu

ggsave('figures/hr-mean-years.png', p_mu, width = 8, height = 8, dpi = 600,
       bg = 'white')

# variance in HR ---
p_s <-
  ggplot(preds_s2_years, aes(group = sex_treatment)) +
  facet_grid(sex ~ study_year) +
  geom_vline(xintercept = as.Date('2021-11-09'), col = 'red')+
  geom_ribbon(aes(date, ymin = sqrt(lwr_95), ymax = sqrt(upr_95),
                  fill = site), alpha = 0.3) +
  geom_line(aes(date, sqrt(s2), color = site)) +
  scale_color_brewer('Site', type = 'qual', palette = 2,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  labs(x = NULL,
       y = expression(bold('SD in 7-day home range size'~(km^2)))) +
  theme(legend.position = 'top'); p_s

ggsave('figures/hr-sd-years.png', p_s, width = 8, height = 8, dpi = 600,
       bg = 'white')
