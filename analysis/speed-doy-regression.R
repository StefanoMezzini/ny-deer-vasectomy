library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('mgcv')      # for modeling
library('lubridate') # for working with dates
library('ggplot2')    # for fancy plots
library('gratia')    # for predicting from models
library('ctmm')  # for movement modeling
source('functions/gammals-variance-simulation-cis.R')
source('analysis/ref_dates.R')
source('analysis/figures/default-theme.R')

d <- readRDS('data/years-1-and-2-data-no-akde.rds') %>%
  filter(is.finite(speed_est)) %>% # loosing ~40% of the windows
  mutate(dof_speed = map_dbl(model, ctmm:::DOF.speed))

range(d$dof_speed) # range may be ok
hist(d$dof_speed)

# plot the data ----
p_speed <-
  mutate(d,
         sex = if_else(sex == 'f', 'females', 'males'),
         treatment = if_else(study_site == 'rockefeller', 'Rockefeller',
                             'Staten Island'),
         t_s = paste(treatment, sex),
         study_year = paste('Year', study_year)) %>%
  ggplot(aes(date, speed_est, group = animal)) +
  facet_grid(t_s ~ study_year, scales = 'free_x') +
  geom_vline(xintercept = REF_DATES, col = 'red') +
  geom_line() +
  labs(x = NULL,
       y = expression(bold('Estimated weekly speed'~(km/day))))
p_speed

ggsave('figures/speed-estimates.png', p_speed, width = 8,
       height = 8, dpi = 600, bg = 'white')

# fit a Hierarchical Generalized Additive Model for Location and Scale ----
# not using cyclic cubic splines because the gap is too big
# for year 1 the gap is even bigger
range(d$days_since_aug_1) # not close to 0 to 365
365 - diff(range(d$days_since_aug_1))

# there are some high speed values, but they are well explained by the
# model. The slight overdispersion is caused by low speeds more than high
# speeds 

# location-scale model ----
#' using `bs = 'ad'` worsens AIC and BIC by >1500, even with `k = 30`
if(FALSE) {
  m_speed <- gam(formula = list(
    # linear predictor for the mean
    speed_est ~
      # temporal sex- and treatment-level trends with different
      s(days_since_aug_1, by = sex_treatment, k = 10) +
      # accounts for differences in trends between years
      s(days_since_aug_1, by = sex_treatment, study_year, k = 10, bs = 'sz') +
      # accounts for differences between individuals
      s(animal, bs = 're'),
    
    # linear predictor for the scale (sigma2 = mu^2 * scale)
    # allows mean-variance relationship to be different between sexes
    # sex- and treatment-level trends over season
    ~ s(days_since_aug_1, by = sex_treatment, k = 10)),
    
    family = gammals(),
    data = d,
    method = 'REML',
    control = gam.control(trace = TRUE))
  
  appraise(m_speed, method = 'simulate', n_bins = 30, point_alpha = 0.1)
  plot(m_speed, pages = 1) #' `gratia::draw()` doesn't support `bs = 'sz'`
  
  saveRDS(m_speed, paste0('models/m_speed-hgamls-', Sys.Date(), '.rds'))
} else {
  m_speed <- readRDS('models/m_speed-hgamls-2024-04-04.rds')
}

appraise(m_speed, method = 'simulate')
plot(m_speed, pages = 1)
summary(m_speed)

# check what groups cause oddly large outliers
d %>%
  ungroup() %>%
  mutate(sex = if_else(sex == 'f', 'Females', 'Males'),
       treatment = if_else(study_site == 'rockefeller', 'Rockefeller',
                           'Staten Island'),
       fitted = m_speed$fitted.values[, 1],
       large = resid(m_speed) > 3) %>%
  ggplot() +
  facet_grid(treatment ~ sex) +
  geom_point(aes(fitted, speed_est, color = large, alpha = large)) +
  labs(x = 'Fitted values', y = 'Observed values') +
  scale_color_manual('Deviance residuals > 3', values = 1:2) +
  scale_alpha_manual('Deviance residuals > 3', values = c(0.3, 1)) +
  theme(legend.position = 'top')

ggsave('figures/speed-model-obs-fitted.png', width = 8, height = 8,
       dpi = 600, bg = 'white')

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

preds_mu <- gammals_mean(model = m_speed, data = newd, nsims = 1e4,
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

preds_s <- gammals_var(model = m_speed, data = newd, nsims = 1e4,
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

# mean speed
DATES <- as.Date(c('2021-09-15', '2021-11-15', '2022-01-15', '2022-03-15',
                   '2022-05-15'))
LABS <- format(DATES, '%B 15')
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
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab('Mean distance travelled (km/day)') +
  theme(legend.position = 'top'); p_mu

ggsave('figures/speed-mean.png', p_mu, width = 8, height = 8, dpi = 600,
       bg = 'white')

# standard deviation in speed ---
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
  ylab('SD in distance travelled (km/day)') +
  theme(legend.position = 'top'); p_s

ggsave('figures/speed-sd.png', p_s, width = 8, height = 8, dpi = 600,
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
  gammals_mean(model = m_speed, data = newd_years, nsims = 1e4,
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
  gammals_var(model = m_speed, data = newd_years, nsims = 1e4,
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

# mean speed
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
  ylab('Mean distance travelled (km/day)') +
  theme(legend.position = 'top'); p_mu_y

ggsave('figures/speed-mean-years.png',
       p_mu_y, width = 16, height = 8, dpi = 600, bg = 'white')

# variance in speed ---
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
  ylab('SD in distance travelled (km/day)') +
  theme(legend.position = 'top'); p_s

ggsave('figures/speed-sd-years.png', p_s_y, width = 16, height = 8, dpi = 600,
       bg = 'white')
