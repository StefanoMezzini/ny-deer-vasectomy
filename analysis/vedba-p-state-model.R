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
         p_day = percent_day  / 100,
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
  group_by(animal_year) %>%
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
           as.numeric()) # convert difftime to numeric

# remaining proportions of high activity are ok
d %>%
  filter(state == 'no activity') %>%
  filter(percent_day > 10) %>%
  group_by(animal, study_year) %>%
  summarise(n_days = n(), which_days = paste(date, collapse = ', '))

ggplot(d, aes(days_since_aug_1, percent_day)) +
  facet_grid(study_site ~ sex + study_year) +
  geom_point(aes(color = state), alpha = 0.2) +
  geom_smooth(aes(group = state, fill = state), color = 'black', lwd = 2) +
  geom_smooth(aes(color = state), se = FALSE) +
  scale_color_manual('State', values = pal,
                     aesthetics = c('color', 'fill')) +
  labs(x = expression(Days~since~August~1^{st}),
       y = 'Time spent per day (%)') +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = 'top')

ggplot(d, aes(p_day, fill = state, color = state)) +
  facet_grid(state ~ sex_treatment, scales = 'free') +
  geom_density(alpha = 0.5, show.legend = FALSE) +
  scale_color_manual('State', values = pal, aesthetics = c('color', 'fill'))

min(d$p_day)

d %>%
  group_by(sex_treatment, state) %>%
  summarise(p_zero = mean(p_day == 0))

min(filter(d, p_day > 0)$p_day) # keep EPS < min(p | p > 0)
EPS <- 1e-6
min(filter(d, p_day > 0)$p_day) / EPS

if(file.exists('models/m_active-hgam.rds')) {
  m_active <- readRDS('models/m_active-hgam.rds')
} else {
  m_active <- bam(
    p_day ~
      # by smooths require a separate explicit intercept for each group
      sex_treatment +
      # temporal sex- and treatment-level trends with different smoothness
      s(days_since_aug_1, by = sex_treatment, k = 10, bs = 'tp') +
      # accounts for deviation from average between years
      s(days_since_aug_1, study_year, by = sex_treatment, k = 10, bs = 'sz') +
      # invidual- and year-level deviations from the average
      s(days_since_aug_1, animal_year, by = sex_treatment, k = 10, bs = 'fs',
        xt = list(bs = 'cr')),
    family = betar(link = 'logit', eps = EPS), #' `eps` because of 0s
    data = d,
    subset = state == 'no activity',
    method = 'fREML',
    discrete = TRUE,
    control = gam.control(trace = TRUE))
  
  draw(m_active, rug = FALSE, parametric = TRUE)
  beepr::beep()
  summary(m_active)
  saveRDS(m_active, paste0('models/m_active-hgam-', Sys.Date(), '.rds'))  
}

# group no activity into low (and high into medium) to only have 2 states
d_joined <- d %>%
  filter(state %in% c('no activity', 'low')) %>%
  pivot_wider(names_from = state, values_from = p_day) %>%
  mutate(p_low = `no activity` + low) %>%
  select(- c(`no activity`, low))

if(file.exists('models/m_low-hgam.rds')) {
  m_low <- readRDS('models/m_low-hgam.rds')
} else {
  m_low <- bam(
    p_day ~
      # by smooths require a separate explicit intercept for each group
      sex_treatment +
      # temporal sex- and treatment-level trends with different smoothness
      s(days_since_aug_1, by = sex_treatment, k = 10, bs = 'tp') +
      # accounts for deviation from average between years
      s(days_since_aug_1, study_year, by = sex_treatment, k = 10, bs = 'sz') +
      # invidual- and year-level deviations from the average
      s(days_since_aug_1, animal_year, by = sex_treatment, k = 10, bs = 'fs',
        xt = list(bs = 'cr')),
    family = betar(link = 'logit'), #' higher `eps` because of 0s
    data = d,
    subset = state == 'low',
    method = 'fREML',
    discrete = TRUE,
    control = gam.control(trace = TRUE))
  
  draw(m_low, rug = FALSE)
  beepr::beep()
  summary(m_low)
  saveRDS(m_low, paste0('models/m_low-hgam-', Sys.Date(), '.rds'))
}

# newd <-
#   expand.grid(date = seq(as.Date('2021-09-01'),
#                          as.Date('2022-05-30'),
#                          length.out = 400),
#               state = unique(d$state),
#               sex_treatment = unique(d$sex_treatment),
#               study_year = 'new year',
#               animal_year = 'new animal') %>%
#   mutate(days_since_aug_1 = if_else(
#     # if in Aug, Sept, Oct, Nov, or Dec
#     month(date) >= 8,
#     # then calculate days since Aug 1
#     date - as.Date(paste0(year(date), '-08-01')),
#     # otherwise calculate days since Aug 1 of previous year 
#     date - as.Date(paste0(year(date) - 1, '-08-01'))) %>%
#       as.numeric(),
#     sex = substr(sex_treatment, 1, 1),
#     site = substr(sex_treatment, 3,
#                   nchar(as.character(sex_treatment))))
# 
# preds <- bind_cols(newd,
#                    predict(m, newdata = newd, type = 'link',
#                            se.fit = TRUE, unconditional = FALSE) %>%
#                      bind_cols()) %>%
#   # calculate predictions on probability scale
#   mutate(mu = m$family$linkinv(fit) * 100,
#          lwr = m$family$linkinv(fit - 1.96 * se.fit) * 100,
#          upr = m$family$linkinv(fit + 1.96 * se.fit) * 100) %>%
#   # make sure probabilities add up to 1
#   group_by(sex_treatment, days_since_aug_1) %>%
#   mutate(mu = mu / sum(mu) * 100,
#          lwr = lwr / sum(lwr) * 100,
#          upr = upr / sum(upr) * 100) %>%
#   ungroup() %>%
#   # add variables for faceting
#   mutate(sex = substr(sex_treatment, 1, 1),
#          study_site = substr(sex_treatment, 3,
#                              nchar(as.character(sex_treatment))))
# 
# # create figures
# DATES <- as.Date(c('2021-10-01', '2021-12-01', '2022-02-01', '2022-04-01'))
# DATES_AUG_1 <- if_else(
#   # if in Aug, Sept, Oct, Nov, or Dec
#   month(DATES) >= 8,
#   # then calculate days since Aug 1
#   DATES - as.Date(paste0(year(DATES), '-08-01')),
#   # otherwise calculate days since Aug 1 of previous year 
#   DATES - as.Date(paste0(year(DATES) - 1, '-08-01'))) %>%
#   as.numeric()
# LABS <- format(DATES, '%B 1')
# 
# ggplot(d, aes(days_since_aug_1)) +
#   facet_grid(study_site ~ sex) +
#   geom_point(aes(y = percent_day, color = state), alpha = 0.1) +
#   geom_ribbon(aes(ymin = lwr, ymax = upr, group = state), preds,
#               alpha = 0.3) +
#   geom_line(aes(y = mu, group = state), preds, color = 'black', lwd = 1.5) +
#   geom_line(aes(y = mu, color = state), preds) +
#   scale_color_manual('State', values = pal,
#                      aesthetics = c('color', 'fill')) +
#   scale_x_continuous(NULL, breaks = DATES_AUG_1, labels = LABS,
#                      limits = c(61, 272)) +
#   ylab('Time spent per day (%)') +
#   guides(colour = guide_legend(override.aes = list(alpha = 1))) +
#   theme(legend.position = 'top')
# 
# ggplot() +
#   facet_grid(study_site ~ sex) +
#   geom_area(aes(date, mu, color = state, fill = state), preds,
#             position = position_stack(reverse = TRUE), alpha = 0.5, lwd = 1) +
#   scale_fill_manual('State', values = pal,
#                     aesthetics = c('color', 'fill')) +
#   guides(colour = guide_legend(override.aes = list(alpha = 1))) +
#   scale_x_continuous(NULL, breaks = DATES, labels = LABS,
#                      limits = as.Date(c('2021-10-01', '2022-04-30'))) +
#   scale_y_continuous('Time spent per day (%)', expand = c(0, 0)) +
#   theme(legend.position = 'top')
