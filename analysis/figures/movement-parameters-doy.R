library('dplyr')     # for data wrangling
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
library('cowplot')   # for fancy plots in grids
source('analysis/ref_dates.R')
source('analysis/figures/default-theme.R')

# date breaks and labels
DATES <- as.Date(c('2021-10-01', '2021-12-01', '2022-02-01', '2022-04-01'))
LABS <- format(DATES, '%B 1')

p_hr <-
  readRDS('models/predictions/hr-preds_mu.rds') %>%
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
       width = 16, height = 12 * 1.05, dpi = 600, bg = 'white')

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
                       aesthetics = c('color', 'fill')) +
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
       width = 16, height = 8 * 1.05, dpi = 600, bg = 'white')
