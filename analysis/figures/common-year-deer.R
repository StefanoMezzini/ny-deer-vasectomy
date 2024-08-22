library('dplyr')   # for data wrangling
library('ggplot2') # for fancy plots
library('cowplot') # for fancy multi-panel plots
source('analysis/figures/default-theme.R')

d <- readRDS('data/years-1-and-2-data-no-akde.rds')

# find deer that were collared twice
common <- d %>%
  group_by(study_year2, animal) %>%
  slice(1) %>%
  ungroup() %>%
  filter(duplicated(animal)) %>%
  pull(animal)

d <- filter(d, animal %in% common) %>%
  arrange(study_site, sex, animal) %>%
  mutate(sex = if_else(sex == 'f', 'Female', 'Male'),
         study_site = if_else(study_site == 'rockefeller', 'Rockefeller',
                              'Staten Island'),
         id = paste(study_site, sex, animal, sep = ', ') %>%
           factor(., levels = unique(.)),
         study_year2 = factor(study_year2))

p_hr <-
  ggplot(d, aes(days_since_aug_1, hr_est_95, group = study_year2,
                color = study_year2)) +
  facet_wrap(~ id, scales = 'free_y', drop = TRUE, nrow = 2) +
  geom_line() +
  labs(x = expression(bold('Days since August 1')^st),
       y = expression(bold('Estimated 7-day home-range size'~
                             (km^2)))) +
  scale_color_brewer('Study year', type = 'qual', palette = 2,
                     aesthetics = c('color', 'fill')) +
  theme(legend.position = 'none')

p_speed <-
  ggplot(d, aes(days_since_aug_1, speed_est, group = study_year2,
                color = study_year2)) +
  facet_wrap(~ id, scales = 'free_y', drop = TRUE, nrow = 2) +
  geom_line() +
  labs(x = expression(bold('Days since August 1')^st),
       y = 'Estimated distance traveled (km/day)') +
  scale_color_brewer('Study year', type = 'qual', palette = 2,
                     aesthetics = c('color', 'fill')) +
  theme(legend.position = 'none')

p_diffusion <-
  ggplot(d, aes(days_since_aug_1, diffusion_est, group = study_year2,
                color = study_year2)) +
  facet_wrap(~ id, scales = 'free_y', drop = TRUE, nrow = 2) +
  geom_line() +
  labs(x = expression(bold('Days since August 1')^st),
       y = expression(bold('Estimated diffusion'~(km^2/day)))) +
  scale_color_brewer('Study year', type = 'qual', palette = 2,
                     aesthetics = c('color', 'fill')) +
  theme(legend.position = 'none')

p <-
  plot_grid(get_legend(p_hr + theme(legend.position = 'top')),
          p_hr, p_speed, p_diffusion, ncol = 1, rel_heights = c(0.2, 1, 1, 1),
          labels = c('', 'A', 'B', 'C'))

ggsave('figures/twice-collared-deer-movement-data.png', p,
       width = 16, height = 16, units = 'in', dpi = 600, bg = 'white')
