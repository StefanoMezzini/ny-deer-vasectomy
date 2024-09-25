library('dplyr')     # for data wrangling
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
library('cowplot')   # for fancy plots in grids
source('analysis/ref_dates.R')
source('analysis/figures/default-theme.R')

# date breaks and labels
DATES <- as.Date(c('2021-10-01', '2021-12-01', '2022-02-01', '2022-04-01'))
LABS <- format(DATES, '%B 1')

d <- readRDS('data/years-1-and-2-data-no-akde.rds') %>%
  mutate(sex = if_else(sex == 'f', 'Female', 'Male'),
         site = if_else(study_site == 'staten_island', 'Staten Island',
                        'Rockefeller'))

# without data ----
p_hr <-
  readRDS('models/predictions/hr-preds_mu.rds') %>%
  ggplot() +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill'),
                     labels = c('Control site', 'Treatment site')) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab(expression(bold('Mean 7-day HR size'~(km^2)))) +
  theme(legend.position = 'none')

p_speed <-
  readRDS('models/predictions/speed-preds_mu.rds') %>%
  ggplot() +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab('Mean distance traveled (km/day)') +
  theme(legend.position = 'none')

p_diff <-
  readRDS('models/predictions/diffusion-preds_mu.rds') %>%
  ggplot() +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab(expression(bold('Mean diffusion'~(km^2/day)))) +
  theme(legend.position = 'none')

p_exc <-
  readRDS('models/predictions/density-preds_mu.rds') %>%
  ggplot() +
  facet_grid(sex ~ .) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab('Mean daily excursivity') +
  theme(legend.position = 'none')

plot_grid(get_plot_component(p_hr + theme(legend.position = 'top'),
                             'guide-box-top'),
          plot_grid(p_hr, p_speed, p_diff, p_exc, labels = 'AUTO',
                    ncol = 2),
          ncol = 1, rel_heights = c(1, 20))

ggsave('figures/mean-movement-parameters-doy-faceted.png',
       width = 12, height = 9 * 1.05, dpi = 600, bg = 'white')

# with data ----
d_exc <-
  # find density values using the models fit to the full telemetry
  readRDS('models/full-telemetry-movement-models-2024-04-20.rds') %>%
  transmute(animal,
            hr_est_95,
            tel = purrr::map(tel, \(.t) {
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
  tidyr::unnest(tel) %>%
  as_tibble() %>%
  mutate(days_since_aug_1 =
           if_else(
             # if in Aug, Sept, Oct, Nov, or Dec
             month(date) >= 8,
             # then calculate days since Aug 1
             date - as.Date(paste0(year(date), '-08-01')),
             # otherwise calculate days since Aug 1 of previous year 
             date - as.Date(paste0(year(date) - 1, '-08-01')))%>%
           as.numeric()) %>% # convert difftime to numeric
  mutate(sex = if_else(sex == 'f', 'Female', 'Male'),
         site = if_else(study_site == 'staten_island', 'Staten Island',
                        'Rockefeller'))

p_hr <-
  readRDS('models/predictions/hr-preds_mu.rds') %>%
  ggplot() +
  facet_grid(sex ~ site) +
  geom_point(aes(as_date('2021-08-01') + days_since_aug_1, hr_est_95),
             alpha = 0.1, filter(d, hr_est_95 <= 10)) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill'),
                     labels = c('Control site', 'Treatment site')) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab(expression(bold('Mean 7-day HR size'~(km^2)))) +
  theme(legend.position = 'none')

p_speed <-
  readRDS('models/predictions/speed-preds_mu.rds') %>%
  ggplot() +
  facet_grid(sex ~ site) +
  geom_point(aes(as_date('2021-08-01') + days_since_aug_1, speed_est),
             alpha = 0.1, d) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab('Mean distance traveled (km/day)') +
  theme(legend.position = 'none')

p_diff <-
  readRDS('models/predictions/diffusion-preds_mu.rds') %>%
  ggplot() +
  facet_grid(sex ~ site) +
  geom_point(aes(as_date('2021-08-01') + days_since_aug_1, diffusion_est),
             alpha = 0.1, d) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab(expression(bold('Mean diffusion'~(km^2/day)))) +
  theme(legend.position = 'none')

p_exc <-
  readRDS('models/predictions/density-preds_mu.rds') %>%
  ggplot() +
  facet_grid(sex ~ site) +
  geom_point(aes(as_date('2021-08-01') + days_since_aug_1, density),
             alpha = 0.1, d_exc) +
  geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
  geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site),
              alpha = 0.3) +
  geom_line(aes(date, mu, color = site), lwd = 1) +
  scale_color_brewer('Site', type = 'qual', palette = 1,
                     aesthetics = c('color', 'fill')) +
  scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                     limits = as.Date(c('2021-10-01', '2022-04-30'))) +
  ylab('Mean daily excursivity') +
  theme(legend.position = 'none')

plot_grid(get_plot_component(p_hr + theme(legend.position = 'top'),
                             'guide-box-top'),
          plot_grid(p_hr, p_speed, p_diff, p_exc, labels = 'AUTO',
                    ncol = 2),
          ncol = 1, rel_heights = c(1, 20))

ggsave('figures/mean-movement-parameters-doy-faceted-with-data.png',
       width = 24, height = 9 * 1.05, dpi = 600, bg = 'white')

# single figure for all parameters with no faceting for ease of comparison ----
make_plot <- function(filename) {
  readRDS(filename) %>%
    mutate(
      sex = factor(sex, levels = c('Female', 'Male'))) %>%
    ggplot() +
    geom_vline(xintercept = REF_DATES[c(1, 3)], col = 'red') +
    geom_ribbon(aes(date, ymin = lwr_95, ymax = upr_95, fill = site,
                    group = sex_treatment), alpha = 0.3) +
    geom_line(aes(date, mu, color = site, lty = sex), lwd = 1) +
    scale_color_brewer('Site', type = 'qual', palette = 1,
                       aesthetics = c('color', 'fill'),
                       labels = c('Control site', 'Treatment site')) +
    theme(legend.position = 'top') +
    scale_x_continuous(NULL, breaks = DATES, labels = LABS,
                       limits = as.Date(c('2021-10-01', '2022-04-30'))) +
    scale_linetype_manual('Sex', values = c(2, 1)) +
    theme(legend.position = 'none')
}

plot_grid(
  get_plot_component(
    make_plot('models/predictions/hr-preds_mu.rds') +
      theme(legend.position = 'top',
            legend.key.width = rel(2)),
    pattern = 'guide-box-top', return_all = TRUE),
  plot_grid(
    make_plot('models/predictions/hr-preds_mu.rds') +
      ylab(expression(bold('Mean 7-day HR size'~(km^2)))),
    make_plot('models/predictions/speed-preds_mu.rds') +
      ylab('Mean distance traveled (km/day)'),
    make_plot('models/predictions/diffusion-preds_mu.rds') +
      ylab(expression(bold('Mean diffusion'~(km^2/day)))),
    make_plot('models/predictions/density-preds_mu.rds') +
      ylab('Mean daily excursivity'),
    labels = 'AUTO'),
  ncol = 1, rel_heights = c(1, 20))

ggsave('figures/mean-movement-parameters-doy.png',
       width = 12, height = 6 * 1.05, dpi = 600, bg = 'white')
