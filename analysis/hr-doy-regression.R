library('dplyr')     # for data wrangling
library('mgcv')      # for modeling
library('lubridate') # for working with dates
library('ggplot2')    # for fancy plots
library('gratia')    # for predicting from models
theme_set(theme_bw())

d <- readRDS('data/years-1-and-2-data.rds') %>%
  mutate(animal = factor(animal),
         study_year = factor(study_year),
         sex_treatment = factor(paste(sex, study_site)),
         s_t_y = factor(paste(sex, study_site, study_year)))

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
  coord_cartesian(ylim = c(0, 15)) +
  geom_line()

# using free scales to view cycles more clearly
ggplot(d, aes(days_since_aug_1, hr_est_95)) +
  facet_wrap(~ sex_treatment, scales = 'free_y') +
  geom_vline(xintercept = 100, col = 'red') +
  geom_line(aes(group = animal)) +
  geom_smooth()

# scaling by max home range
d %>%
  group_by(animal) %>%
  mutate(hr_rel = hr_est_95 / max(hr_est_95)) %>%
  ungroup() %>%
  ggplot(aes(days_since_aug_1, hr_rel, group = animal)) +
  facet_wrap(~ sex_treatment, scales = 'free_y') +
  geom_vline(xintercept = 100, col = 'red') +
  geom_line()

# find odd male from S Island
filter(d, sex_treatment == 'm staten_island', hr_est_95 > 15) %>%
  select(study_year, animal, Age)

# find males with large HRs on Rockefeller
filter(d, sex_treatment == 'm rockefeller', hr_est_95 > 7.5) %>%
  select(study_year, animal, Age)

# start by fitting a Hierarchical Generalized Additive Model ----
m_hr <- bam(hr_est_95 ~
              # temporal sex- and treatment-level trends
              s(days_since_aug_1, sex_treatment, k = 25, bs = 'fs') +
              # accounts for differences between years
              s(days_since_aug_1, study_year, k = 5, bs = 'fs') +
              # random intercept for individual
              s(animal, bs = 're'),
            family = Gamma(link = 'log'),
            data = d,
            weights = weight,
            method = 'fREML',
            discrete = TRUE,
            control = gam.control(trace = TRUE))

summary(m_hr)

plot(m_hr, pages = 1, trans = exp)

preds <- expand.grid(date = seq(as.Date('2021-09-01'),
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
                       nchar(as.character(sex_treatment)))) %>%
  bind_cols(.,
            predict(object = m_hr, newdata = ., type = 'link',
                    se.fit = TRUE, discrete = FALSE) %>%
              as.data.frame()) %>%
  mutate(mu = exp(fit),
         lwr = exp(fit - 1.96 * se.fit),
         upr = exp(fit + 1.96 * se.fit))

ggplot(preds, aes(group = sex_treatment)) +
  facet_grid(study_year ~ .) +
  geom_vline(xintercept = as.Date('2021-11-09'), col = 'red') +
  geom_ribbon(aes(date, ymin = lwr, ymax = upr,
                  fill = sex), alpha = 0.2) +
  geom_line(aes(date, mu, color = sex, lty = site)) +
  scale_color_brewer('Sex', type = 'qual', palette = 6,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site name') +
  labs(x = NULL, y = expression('7-day home range size'~(km^2))) +
  theme(legend.position = 'top')

# location-scale model ----
m_hr_2 <- gam(list(
  
  # linear predictor for the mean
  hr_est_95 ~
    # temporal sex- and treatment-level trends
    s(days_since_aug_1, sex_treatment, k = 25, bs = 'fs') +
    # accounts for differences between years
    s(days_since_aug_1, study_year, k = 5, bs = 'fs') +
    # random intercept for individual
    s(animal, bs = 're'),
  
  # linear predictor for the scale (sigma2 = mu^2 * scale)
  ~ sex_treatment), # sex- and treatment-level intercepts
  family = gammals(),
  data = d,
  weights = weight,
  method = 'REML',
  control = gam.control(trace = TRUE))

saveRDS(m_hr_2, paste0('models/m_hr_2-', Sys.Date(), '.rds'))

plot(m_hr_2, pages = 1, all.terms = TRUE)

preds_2 <-
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
                       nchar(as.character(sex_treatment)))) %>%
  bind_cols(.,
            predict(object = m_hr_2, newdata = ., type = 'response',
                    se.fit = TRUE) %>%
              as.data.frame()) %>%
  mutate(mu = fit.1,
         sigma2 = fit.1^2 * exp(fit.2))

#' **NOTE:** still need to add credible intervals
# mean HR
ggplot(preds_2, aes(group = sex_treatment)) +
  coord_cartesian(ylim = c(0, 5)) +
  facet_grid(study_year ~ ., ) +
  geom_vline(xintercept = as.Date('2021-11-09'), col = 'red', alpha = 0.3)+
  geom_line(aes(date, mu, color = sex, lty = site)) +
  scale_color_brewer('Sex', type = 'qual', palette = 6,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  labs(x = NULL, y = expression('Mean 7-day home range size'~(km^2))) +
  theme(legend.position = 'top')

# SD in HR
ggplot(preds_2, aes(group = sex_treatment)) +
  coord_cartesian(ylim = c(0, 5)) +
  facet_grid(study_year ~ .) +
  geom_vline(xintercept = as.Date('2021-11-09'), col = 'red', alpha = 0.3)+
  geom_line(aes(date, sqrt(sigma2), color = sex, lty = site)) +
  scale_color_brewer('Sex', type = 'qual', palette = 6,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  labs(x = NULL, y = expression('SD in 7-day home range size'~(km^2))) +
  theme(legend.position = 'top')
