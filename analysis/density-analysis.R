library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
library('mgcv')      # for modeling
library('gratia')    # for ggplot-based diagnostics
source('analysis/ref_dates.R')
source('analysis/figures/default-theme.R')
source('../Directed-Study/functions/betals.r')

d <-
  readRDS('models/full-telemetry-movement-models-2024-04-20.rds') %>%
  transmute(animal,
            hr_est_95,
            tel = map(tel, \(.t) {
              data.frame(.t) %>%
                mutate(date = as.Date(as.POSIXct(timestamp))) %>%
                # calculate median daily density
                group_by(date) %>%
                summarise(density = mean(density), .groups = 'drop') %>%
                # drop first 7 days
                filter(date > min(date) + 7)
            }
            )) %>%
  left_join(readRDS('data/years-1-and-2-data-no-akde.rds') %>%
              group_by(animal) %>%
              slice(1) %>%
              select(animal, study_year, study_year2, study_site, sex, Age,
                     sex_treatment, has_fawn) %>%
              ungroup(),
            by = 'animal') %>%
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

d

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
  geom_line(aes(days_since_aug_1, density, group = animal_year), alpha = 0.1)

# fit the model ----
# 40 minutes with s(by=s_t) + s(s_t, year), s(animal)
# takes ~ 10 seconds
# study_year smooth widens CIs substantially
# fits in ~ 35 minutes
if(file.exists('models/m_density-hgam-2024-05-07.rds')) {
  m_density <- readRDS('models/m_density-hgam-2024-05-07.rds')
} else {
  m_density <- bam(
    density ~
      0 + sex_treatment +
      s(days_since_aug_1, by = sex_treatment, k = 15) +
      s(days_since_aug_1, study_year, by = sex_treatment, k = 15, bs = 'fs') +
      s(days_since_aug_1, animal_year, k = 15, bs = 'fs',
        xt = list(bs = 'cr')),
    family = betar(link = 'logit'), # since CDF is in (0, 1)
    data = d,
    method = 'fREML',
    discrete = TRUE,
    control = gam.control(trace = TRUE))
  saveRDS(m_density, paste0('models/m_density-hgam-', Sys.Date(), '.rds'))
  
  m_density <- gam(formula = list(
    # linear predictor for the mean
    density ~
      0 + sex_treatment +
      # temporal sex- and treatment-level trends with different
      s(days_since_aug_1, by = sex_treatment, k = 15, bs = 'tp') +
      # accounts for differences in trends between years
      s(days_since_aug_1, by = sex_treatment, study_year, k = 15, bs = 'sz') +
      # accounts for differences between individuals
      s(animal, bs = 're'),
    
    # linear predictor for the scale (sigma2 = mu^2 * scale)
    # allows mean-variance relationship to be different between sexes
    # sex- and treatment-level trends over season
    ~ 0 + sex_treatment +
      s(days_since_aug_1, by = sex_treatment, k = 15, bs = 'tp') +
      # accounts for differences between individuals
      s(animal, bs = 're')),
    
    family = gammals(),
    data = d,
    method = 'REML',
    control = gam.control(trace = TRUE))
}

appraise(m_density, point_alpha = 0.05, type = 'pearson')
plot(m_density, pages = 1, scheme = 0)
summary(m_density, re.test = FALSE) # good deviance explained

mutate(d,
       e = resid(m_density)) %>%
  ggplot(aes(days_since_aug_1, e, group = animal)) +
  facet_wrap(study_site ~ sex + study_year) +
  geom_line(alpha = 0.1)

# residuals by group look good
ggplot(mutate(d, e = resid(m_density, type = 'pearson'))) +
  facet_grid(sex ~ study_site, scales = 'free') +
  geom_density(aes(e, fill = study_year, color = study_year), alpha = 0.4) +
  geom_vline(xintercept = 0, color = 'grey20') +
  scale_fill_brewer('Study year', type = 'qual', palette = 6,
                    aesthetics = c('fill', 'color'))

# no problematic trends in the residuals over time
ggplot(mutate(d, e = resid(m_density, type = 'pearson'))) +
  facet_grid(sex ~ study_site + study_year, scales = 'free') +
  geom_point(aes(days_since_aug_1, e), alpha = 0.2) +
  geom_smooth(aes(days_since_aug_1, e))

# create figures ----
newd <-
  expand.grid(date = seq(as.Date('2021-09-01'),
                         as.Date('2022-05-30'),
                         length.out = 400),
              sex_treatment = unique(d$sex_treatment),
              study_year = 'null',
              animal_year = 'null') %>%
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

preds <- bind_cols(
  newd,
  predict(m_density, newdata = newd, unconditional = FALSE, se.fit = TRUE,
          discrete = FALSE,
          exclude = c('s(days_since_aug_1,animal_year)',
                      paste0('s(days_since_aug_1,study_year):sex_treatment',
                             unique(d$sex_treatment)))) %>%
    as.data.frame()) %>%
  mutate(lwr_95 = m_density$family$linkinv(fit - se.fit * 1.96),
         mu = m_density$family$linkinv(fit),
         upr_95 = m_density$family$linkinv(fit + se.fit * 1.96),
         sex = case_when(sex == 'f' ~ 'Female',
                         sex == 'm' ~ 'Male'),
         site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                          site == 'staten_island' ~ 'Staten Island'))

preds_y <-
  bind_rows(mutate(newd, study_year = 1),
            mutate(newd, study_year = 2)) %>%
  bind_cols(
    .,
    predict(m_density, newdata = ., unconditional = FALSE, se.fit = TRUE,
            discrete = FALSE,
            exclude = c('s(days_since_aug_1,animal_year)')) %>%
      as.data.frame()) %>%
  mutate(lwr_95 = m_density$family$linkinv(fit - se.fit * 1.96),
         mu = m_density$family$linkinv(fit),
         upr_95 = m_density$family$linkinv(fit + se.fit * 1.96),
         sex = case_when(sex == 'f' ~ 'Female',
                         sex == 'm' ~ 'Male'),
         site = case_when(site == 'rockefeller' ~ 'Rockefeller',
                          site == 'staten_island' ~ 'Staten Island'))

# create the figures ----
DATES <- as.Date(c('2021-09-15', '2021-11-15', '2022-01-15', '2022-03-15',
                   '2022-05-15'))
LABS <- format(DATES, '%B 15')

p_mu <-
  ggplot(preds, aes(group = sex_treatment)) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab('Mean daily density') +
  theme(legend.position = 'top'); p_mu

ggsave('figures/density-mean.png', p_mu, width = 8, height = 8, dpi = 600,
       bg = 'white')

p_mu_y <-
  ggplot(preds_y, aes(group = sex_treatment)) +
  facet_grid(sex ~ study_year) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_linetype('Site') +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS) +
  ylab('Mean daily density') +
  theme(legend.position = 'top'); p_mu_y

ggsave('figures/density-mean-years.png', p_mu_y, width = 8, height = 8,
       dpi = 600, bg = 'white')
