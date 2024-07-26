library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
source('analysis/figures/default-theme.R')

d <- readRDS('models/full-telemetry-movement-models-2024-04-20.rds') %>%
  transmute(animal,
            Sex = if_else(animal.sex == 'f', 'Female', 'Male'),
            site = if_else(study.area == 'rockefeller', 'Rockefeller', 
                           'Staten Island'),
            study_year = paste('Year', study.year),
            tel = map(tel, data.frame)) %>%
  unnest(tel) %>%
  group_by(animal, Sex, site, study_year) %>%
  summarize(dt = median(diff(timestamp), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(dt = as.numeric(dt) / (1 %#% 'hour'))
d

ggplot(d, aes(animal, dt)) +
  facet_grid(study_year ~ .) +
  geom_text(aes(label = animal), filter(d, dt > 2),
            nudge_y = -0.2, nudge_x = 0.5) +
  geom_vline(aes(xintercept = animal),
             filter(d, duplicated(animal)) %>%
               select(! study_year), alpha = 0.1) +
  geom_point(aes(shape = Sex, color = site)) +
  labs(x = NULL, y = 'Median sampling interval (hours)') +
  scale_color_brewer('Site', type = 'qual', palette = 1) +
  theme(axis.text.x.bottom = element_text(angle = 90, face = 'plain',
                                          size = 7),
        legend.position = 'top')

ggsave('figures/sampling-intervals.png', width = 12, height = 4,
       units = 'in', dpi = 600, bg = 'white')
