library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for working with dates
library('ggplot2')   # for fancy figures
library('khroma')    # for colorblind-friendly color palettes
library('cowplot')   # for fancy multi-panel figures
library('mgcv')      # for modeling
library('gratia')    # for ggplot-based model figures
source('analysis/figures/default-theme.R')
source('analysis/ref_dates.R')

# model number of transitions/day based on TOD, DOY, sex, location, etc.
d <- bind_rows(map_dfr(list.files('data/Year 1/transition-day-data/',
                                  pattern = '.xlsx', full.names = TRUE),
                       readxl::read_xlsx) %>%
                 mutate(study_year = '1'),
               map_dfr(list.files('data/Year 2/transition-day-data/',
                                  pattern = '.xlsx', full.names = TRUE),
                       readxl::read_xlsx) %>%
                 mutate(study_year = '2')) %>%
  rename(n_transitions = n_of_transitions,
         animal = `individual-local-identifier`) %>%
  mutate(animal_year = factor(paste(animal, study_year))) %>%
  filter(
    ! ((animal_year == '1123 2' & date >= as.Date('2023-01-11')) |
         (animal_year == '170 1'  & date >= as.Date('2022-04-20')) |
         (animal_year == '25 2' & date >= as.Date('2022-11-01')) |
         (animal_year == '28 1'  & date >= as.Date('2022-01-27')) |
         (animal_year == '44 1'  & date >= as.Date('2022-02-26')) |
         (animal_year == '5b 1'  & date >= as.Date('2021-10-24')))) %>%
  # drop first and last rows since they may have less than 24 hours of data
  group_by(animal_year) %>%
  arrange(date) %>%
  slice(- c(1:2, n())) %>%
  ungroup() %>%
  # add other variables
  mutate(study_year = factor(study_year),
         sex_treatment = factor(paste(sex, study_site)),
         days_since_aug_1 =
           if_else(
             # if in Aug, Sept, Oct, Nov, or Dec
             month(date) >= 8,
             # then calculate days since Aug 1
             as.Date(date) - as.Date(paste0(year(date), '-08-01')),
             # otherwise calculate days since Aug 1 of previous year 
             as.Date(date) - as.Date(paste0(year(date) - 1, '-08-01'))) %>%
           as.numeric()) # convert difftime to numeric

# deer 91 in year 2 likely has a malfunctioning collar
filter(d, n_transitions < 50 * 24)

ggplot(d, aes(days_since_aug_1, n_transitions / 24)) +
  facet_grid(study_site ~ sex + study_year) +
  geom_point(aes(color = animal_year == '91 2'), alpha = 0.2) +
  geom_hline(yintercept = 50) +
  scale_color_manual('Deer 91 in year 2', values = 1:2) +
  labs(x = expression(Days~since~August~1^{st}),
       y = 'Number of transitions per hour') +
  theme(legend.position = 'top')

d <- filter(d, animal_year != '91 2')

# histogram of n transitions by day
ggplot(d, aes(n_transitions / 24)) +
  facet_grid(study_site ~ sex + study_year) +
  geom_histogram(fill = 'grey', color = 'black') +
  labs(x = 'Transitions per hour', y = 'Count',
       title = 'Values < 50 should be ok')

ggplot(d, aes(days_since_aug_1, n_transitions / 24)) +
  facet_grid(study_site ~ sex + study_year) +
  geom_point(alpha = 0.2) +
  geom_smooth() +
  labs(x = expression(Days~since~August~1^{st}),
       y = 'Number of transitions per hour') +
  theme(legend.position = 'top')

#' `n_transitions` is highly overdispersed for Pois because `V(Y) >> E(Y)`
if(file.exists('models/m_n_transitions-hgam-2024-06-11.rds')) {
  m_n_transitions <- readRDS('models/m_n_transitions-hgam-2024-06-11.rds')
} else {
  m_n_transitions <- bam(
    n_transitions ~
      # by smooths require a separate explicit intercept for each group
      sex_treatment +
      # temporal sex- and treatment-level trends with different smoothness
      s(days_since_aug_1, by = sex_treatment, k = 10, bs = 'tp') +
      # accounts for deviation from average between years
      s(days_since_aug_1, study_year, by = sex_treatment, k = 10, bs = 'sz') +
      # invidual- and year-level deviations from the average
      s(days_since_aug_1, animal_year, k = 10, bs = 'fs',
        xt = list(bs = 'cr')),
    family = nb(link = 'log'), # data is clearly over-dispersed
    data = d,
    method = 'fREML',
    discrete = TRUE,
    control = gam.control(trace = TRUE))
  
  saveRDS(m_n_transitions, paste0('models/m_n_transitions-hgam-', Sys.Date(), '.rds'))
}

appraise(m_n_transitions, point_alpha = 0.01)
draw(m_n_transitions, ncol = 4, scales = 'fixed')

# check the spread of the residuals
d %>%
  mutate(e = residuals(m_n_transitions)) %>%
  group_by(sex, study_site) %>%
  summarize(x_bar_e = mean(e),
            s2_e = var(e),
            .groups = 'drop')

d %>%
  mutate(e = residuals(m_n_transitions)) %>%
  ggplot(aes(e)) +
  facet_grid(study_year ~ sex) +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_density(aes(color = study_site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1)

# plot the estimated common trends between the two years ----
newd <-
  expand_grid(date = seq(as.Date('2021-09-01'),
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
                  nchar(as.character(sex_treatment))),
    sex = case_when(sex == 'f' ~ 'Female',
                    sex == 'm' ~ 'Male'),
    site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                     site == 'staten_island' ~ 'Staten Island'))

preds <- bind_cols(
  newd,
  predict(m_n_transitions, newdata = newd, type = 'link',
          se.fit = TRUE, unconditional = FALSE,
          discrete = FALSE,
          exclude = c(paste0('s(days_since_aug_1,study_year):sex_treatment',
                             unique(d$sex_treatment)),
                      's(days_since_aug_1,animal_year)'))) %>%
  mutate(mu = exp(fit),
         lwr = exp(fit - 1.96 * se.fit),
         upr = exp(fit + 1.96 * se.fit))

# create the plots
DATES <- as.Date(c('2021-10-01', '2021-12-01', '2022-02-01', '2022-04-01'))
LABS <- format(DATES, '%B 1')

p_mean <-
  ggplot(preds, aes(group = sex_treatment)) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr, ymax = upr, fill = site),
              alpha = 0.5) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab('Transitions between states in a day') +
  theme(legend.position = 'top'); p_mean

ggsave('figures/n-transitions.png',
       p_mean, width = 8, height = 8, dpi = 600, bg = 'white')

# plot the estimated trends for each year ----
newd_years <- bind_rows(mutate(newd, study_year = 1),
                        mutate(newd, study_year = 2))

preds_years <- bind_cols(
  newd_years,
  predict(m_n_transitions, newdata = newd_years, type = 'link',
          se.fit = TRUE, unconditional = FALSE,
          discrete = FALSE, exclude = 's(days_since_aug_1,animal_year)')) %>%
  mutate(mu = exp(fit),
         lwr = exp(fit - 1.96 * se.fit),
         upr = exp(fit + 1.96 * se.fit))

p_years <-
  ggplot(preds_years, aes(group = sex_treatment)) +
  facet_grid(sex ~ study_year) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr, ymax = upr, fill = site),
              alpha = 0.5) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab('Transitions between states in a day') +
  theme(legend.position = 'top'); p_years

ggsave('figures/n-transitions-years.png',
       p_years, width = 16, height = 8, dpi = 600, bg = 'white')
