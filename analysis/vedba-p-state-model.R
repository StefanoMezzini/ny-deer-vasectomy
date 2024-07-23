library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for working with dates
library('ggplot2')   # for fancy figures
library('cowplot')   # for fancy multi-panel figures
library('mgcv')      # for modeling
library('gratia')    # for ggplot-based model figures
library('khroma')    # for colorblind-friendly color palettes
source('analysis/figures/default-theme.R')
source('analysis/ref_dates.R')

pal <- c('#BBBBBB', unname(color('high contrast')(3)))
plot_scheme(pal)

# model P(state) based on TOD, DOY, sex, location, etc.
d <- bind_rows(map_dfr(list.files('data/Year 1/time-state-per-day-data/',
                                  pattern = '.xlsx', full.names = TRUE),
                       readxl::read_xlsx) %>%
                 mutate(study_year = '1'),
               map_dfr(list.files('data/Year 2/time-state-per-day-data/',
                                  pattern = '.xlsx', full.names = TRUE),
                       readxl::read_xlsx) %>%
                 mutate(study_year = '2')) %>%
  rename(animal = `individual-local-identifier`,
         n = total_n_movements,
         percent_n = `%_total_movements`,
         percent_day = `%_of_time_in_each_state`) %>%
  mutate(weight = n * percent_n,
         p_day = percent_day / 100,
         animal_year = factor(paste(animal, study_year))) %>%
  # remove periods with problematically high proportion of no activity
  filter(
    ! ((animal_year == '1123 2' & date >= as.Date('2023-01-11')) |
         (animal_year == '170 1'  & date >= as.Date('2022-04-20')) |
         (animal_year == '25 2' & date >= as.Date('2022-11-01')) |
         (animal_year == '28 1'  & date >= as.Date('2022-01-27')) |
         (animal_year == '44 1'  & date >= as.Date('2022-02-26')) |
         (animal_year == '5b 1'  & date >= as.Date('2021-10-24')))) %>%
  # drop first 2 and last rows since they may have < 24 h and odd behavior
  group_by(animal_year, state) %>%
  arrange(date) %>%
  slice(- c(1:2, n())) %>%
  ungroup() %>%
  mutate(study_year = factor(study_year),
         state = factor(state, levels = c('no activity', 'low', 'medium',
                                          'high')),
         state_int = case_when(state == 'no activity' ~ 1,
                               state == 'low' ~ 2,
                               state == 'medium' ~ 3,
                               state == 'high' ~ 4),
         sex_treatment = factor(paste(sex, study_site)),
         days_since_aug_1 =
           if_else(
             # if in Aug, Sept, Oct, Nov, or Dec
             month(date) >= 8,
             # then calculate days since Aug 1
             as.Date(date) - as.Date(paste0(year(date), '-08-01')),
             # otherwise calculate days since Aug 1 of previous year 
             as.Date(date) - as.Date(paste0(year(date) - 1, '-08-01'))) %>%
           as.numeric()) %>% # convert difftime to numeric
  # drop days with less than 99% of the data for a day
  group_by(animal_year, days_since_aug_1) %>%
  filter(sum(p_day) > 0.99) %>%
  ungroup()

# deer 91 in year 2 likely has a malfunctioning collar
# (see script for number of tansitions per day)
d <- filter(d, animal_year != '91 2')

# remaining proportions of high activity are ok
d %>%
  filter(state == 'no activity') %>%
  filter(percent_day > 10) %>%
  group_by(animal, study_year) %>%
  summarise(n_days = n(), which_days = paste(date, collapse = ', '))

# some days don't have rows for no or high activity
d %>%
  group_by(sex_treatment, state) %>%
  summarise(n = n()) %>%
  mutate(state = state,
         missing = max(n) - n) %>%
  filter(missing != 0)

# only one state is missing for any row
d %>%
  group_by(animal_year, days_since_aug_1) %>%
  summarize(n_misssing = 4 - n_distinct(state),
            .groups = 'drop') %>%
  filter(n_misssing != 0) %>%
  pull(n_misssing) %>%
  unique()

# proportion of "no activity" state is quite low
# proportion of "high" state is negligible
# HGAMs of "no activity" state struggle to converge
#' zeros and 1s are problematic for beta distributions; see `?betar`
d %>%
  group_by(state, sex_treatment) %>%
  summarise(p = mean(p_day),
            p_zero = mean(p_day == 0))

ggplot(d, aes(p_day, fill = state, color = state)) +
  facet_grid(state ~ sex_treatment, scales = 'free') +
  geom_density(alpha = 0.5, show.legend = FALSE) +
  scale_color_manual('State', values = pal, aesthetics = c('color', 'fill'))

# group no activity into low (and high into medium) to only have 2 states
# personally, I (Stefano) don't believe that VeDBA = 0 is a state separate
# from the "low" state because it doesn't happen often enough to accurately
# show when there is no activity (e.g., rest, sleep)
# additionally, deer are in a "high" state even less often, so estimating
# P(state = high) would be too hard and inaccurate
d <- d %>%
  filter(state %in% c('no activity', 'low')) %>%
  select(c(date, animal, sex, study_site, study_year, animal_year,
           sex_treatment, days_since_aug_1, state, p_day)) %>%
  pivot_wider(names_from = state, values_from = p_day) %>%
  #' all days add to a proportion of 0.99 or higher, `NA`s should be ~0
  mutate(`no activity` = if_else(is.na(`no activity`), 0, `no activity`),
         p_low = `no activity` + low) %>%
  select(- c(`no activity`, low))

if(file.exists('models/m_low-hgam.rds')) {
  m_low <- readRDS('models/m_low-hgam.rds')
} else {
  # fits in ~8 minutes
  # k = 15 increases fitting time to ~28 minutes with excessive wiggliness
  m_low <- bam(
    p_low ~
      # by smooths require a separate explicit intercept for each group
      sex_treatment +
      # temporal sex- and treatment-level trends with different smoothness
      s(days_since_aug_1, by = sex_treatment, k = 10, bs = 'tp') +
      # accounts for deviation from average between years
      s(days_since_aug_1, study_year, by = sex_treatment, k = 10, bs = 'sz') +
      # invidual- and year-level deviations from the average
      s(days_since_aug_1, animal_year, k = 10, bs = 'fs', xt = list(bs = 'cr')),
    family = betar(link = 'logit'),
    data = d,
    method = 'fREML',
    discrete = TRUE,
    control = gam.control(trace = TRUE))
  
  draw(m_low, rug = FALSE)
  beepr::beep()
  summary(m_low)
  appraise(m_low, point_alpha = 0.1, type = 'pearson')
  saveRDS(m_low, 'models/m_low-hgam.rds')
}

# make figures ----
# create figures
DATES <- as.Date(c('2021-10-01', '2021-12-01', '2022-02-01', '2022-04-01'))
DATES_AUG_1 <- if_else(
  # if in Aug, Sept, Oct, Nov, or Dec
  month(DATES) >= 8,
  # then calculate days since Aug 1
  DATES - as.Date(paste0(year(DATES), '-08-01')),
  # otherwise calculate days since Aug 1 of previous year
  DATES - as.Date(paste0(year(DATES) - 1, '-08-01'))) %>%
  as.numeric()
LABS <- format(DATES, '%B 1')

# overall mean
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
      as.numeric())

preds <- bind_cols(
  newd,
  predict(m_low, newdata = newd, type = 'link', se.fit = TRUE,
          unconditional = FALSE, discrete = FALSE,
          terms = c('(Intercept)', 'sex_treatment',
                    paste0('s(days_since_aug_1):sex_treatment',
                           unique(d$sex_treatment)))) %>%
    bind_cols()) %>%
  # calculate predictions on probability scale
  mutate(mu = m_low$family$linkinv(fit),
         lwr = m_low$family$linkinv(fit - 1.96 * se.fit),
         upr = m_low$family$linkinv(fit + 1.96 * se.fit)) %>%
  # add variables for faceting
  mutate(sex = if_else(substr(sex_treatment, 1, 1) == 'f',
                       'Female', 'Male'),
         study_site = if_else(grepl('rockefeller', sex_treatment),
                              'Rockefeller', 'Staten Island'))

d <- mutate(d,
                   sex = if_else(substr(sex_treatment, 1, 1) == 'f',
                                 'Female', 'Male'),
                   study_site = if_else(grepl('rockefeller', sex_treatment),
                                        'Rockefeller', 'Staten Island'))

p <- ggplot(preds, aes(date)) +
  facet_grid(sex ~ .) +
  geom_hline(yintercept = 0.5, color = 'grey', lty = 'dashed') +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = study_site), alpha = 0.5) +
  geom_line(aes(y = mu, color = study_site), lwd = 1) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  scale_y_continuous('Proportion of a day in low-activity state',
                     limits = c(0, 1), expand = c(0, 0)) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = 'top', panel.spacing.y = unit(10, 'points'))

ggsave('figures/p-low.png',
       p, width = 8, height = 8, dpi = 600, bg = 'white')

# differences between years
newd_y <-
  expand_grid(date = seq(as.Date('2021-09-01'),
                         as.Date('2022-05-30'),
                         length.out = 400),
              sex_treatment = unique(d$sex_treatment),
              study_year = unique(d$study_year),
              animal_year = 'new animal') %>%
  mutate(days_since_aug_1 = if_else(
    # if in Aug, Sept, Oct, Nov, or Dec
    month(date) >= 8,
    # then calculate days since Aug 1
    date - as.Date(paste0(year(date), '-08-01')),
    # otherwise calculate days since Aug 1 of previous year
    date - as.Date(paste0(year(date) - 1, '-08-01'))) %>%
      as.numeric())

preds_y <- bind_cols(
  newd_y,
  predict(m_low, newdata = newd_y, type = 'link', se.fit = TRUE,
          unconditional = FALSE, discrete = FALSE,
          terms = c('(Intercept)', 'sex_treatment',
                    paste0('s(days_since_aug_1,study_year):sex_treatment',
                           unique(d$sex_treatment)),
                    paste0('s(days_since_aug_1):sex_treatment',
                           unique(d$sex_treatment)))) %>%
    bind_cols()) %>%
  # calculate predictions on probability scale
  mutate(mu = m_low$family$linkinv(fit),
         lwr = m_low$family$linkinv(fit - 1.96 * se.fit),
         upr = m_low$family$linkinv(fit + 1.96 * se.fit)) %>%
  # add variables for faceting
  mutate(sex = if_else(substr(sex_treatment, 1, 1) == 'f',
                       'Female', 'Male'),
         study_site = if_else(grepl('rockefeller', sex_treatment),
                              'Rockefeller', 'Staten Island'))

p_y <- ggplot(preds_y, aes(date)) +
  facet_grid(sex ~ paste('Year', study_year)) +
  geom_hline(yintercept = 0.5, color = 'grey', lty = 'dashed') +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = study_site), alpha = 0.5) +
  geom_line(aes(y = mu, color = study_site), lwd = 1) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  scale_y_continuous('Proportion of a day in low-activity state',
                     limits = c(0, 1), expand = c(0, 0)) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = 'top', panel.spacing.y = unit(10, 'points'))

ggsave('figures/p-low-years.png',
       p_y, width = 16, height = 8, dpi = 600, bg = 'white')
