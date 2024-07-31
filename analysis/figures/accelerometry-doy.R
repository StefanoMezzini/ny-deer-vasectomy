library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
library('cowplot')   # for fancy plots in grids
source('analysis/ref_dates.R')
source('analysis/figures/default-theme.R')

# date breaks and labels
DATES <- as.Date(c('2021-10-01', '2021-12-01', '2022-02-01', '2022-04-01'))
LABS <- format(DATES, '%B 1')

DATES_AUG_1 <- if_else(
  # if in Aug, Sept, Oct, Nov, or Dec
  month(DATES) >= 8,
  # then calculate days since Aug 1
  DATES - as.Date(paste0(year(DATES), '-08-01')),
  # otherwise calculate days since Aug 1 of previous year
  DATES - as.Date(paste0(year(DATES) - 1, '-08-01'))) %>%
  as.numeric()

# predictions
m_low <- readRDS('models/m_low-hgam.rds')
m_n <- readRDS('models/m_n_transitions-hgam-2024-06-11.rds')

newd <-
  expand_grid(date = seq(as.Date('2021-10-01'),
                         as.Date('2022-05-30'),
                         length.out = 400),
              sex_treatment = unique(m_n$model$sex_treatment),
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
  # P(low)
  predict(m_low, newdata = newd, type = 'link', se.fit = TRUE,
          unconditional = FALSE, discrete = FALSE,
          terms = c('(Intercept)', 'sex_treatment',
                    paste0('s(days_since_aug_1):sex_treatment',
                           unique(newd$sex_treatment)))) %>%
    bind_cols() %>%
    transmute(p_low_mu = m_low$family$linkinv(fit),
              p_low_lwr = m_low$family$linkinv(fit - 1.96 * se.fit),
              p_low_upr = m_low$family$linkinv(fit + 1.96 * se.fit)),
  # n transitions
  predict(m_n, newdata = newd, type = 'link', se.fit = TRUE,
          unconditional = FALSE, discrete = FALSE,
          terms = c('(Intercept)', 'sex_treatment',
                    paste0('s(days_since_aug_1):sex_treatment',
                           unique(newd$sex_treatment)))) %>%
    bind_cols() %>%
    transmute(n_mu = m_n$family$linkinv(fit),
              n_lwr = m_n$family$linkinv(fit - 1.96 * se.fit),
              n_upr = m_n$family$linkinv(fit + 1.96 * se.fit))) %>%
  # add variables for faceting
  mutate(sex = if_else(substr(sex_treatment, 1, 1) == 'f',
                       'Female', 'Male'),
         study_site = if_else(grepl('rockefeller', sex_treatment),
                              'Rockefeller', 'Staten Island'))

p_p_low <-
  ggplot(preds, aes(date)) +
  facet_grid(sex ~ .) +
  geom_hline(yintercept = 0.5, color = 'grey', lty = 'dashed') +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(ymin = p_low_lwr, ymax = p_low_upr, fill = study_site),
              alpha = 0.5) +
  geom_line(aes(y = p_low_mu, color = study_site), lwd = 1) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  scale_y_continuous('Proportion of a day in low-activity state',
                     limits = c(0, 1), expand = c(0, 0)) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = 'none', panel.spacing.y = unit(10, 'points'))

p_n <-
  ggplot(preds, aes(date)) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(ymin = n_lwr, ymax = n_upr, fill = study_site),
              alpha = 0.5) +
  geom_line(aes(y = n_mu, color = study_site), lwd = 1) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  scale_y_continuous('Number of transitions per hour') +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = 'none', panel.spacing.y = unit(10, 'points'))

plot_grid(get_plot_component(p_n + theme(legend.position = 'top'),
                             'guide-box-top'),
          plot_grid(p_p_low, p_n, labels = 'AUTO', nrow = 1),
          ncol = 1, rel_heights = c(1, 10))
ggsave('figures/mean-accelerometry-doy.png',
       width = 16, height = 6 * 1.2, dpi = 600, bg = 'white')
