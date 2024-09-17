library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
library('ctmm')      #' for movement modeling and `%#%`
source('analysis/figures/default-theme.R')

d <- readRDS('models/full-telemetry-movement-models-2024-04-20.rds') %>%
  transmute(animal,
            Sex = if_else(animal.sex == 'f', 'Female', 'Male'),
            site = if_else(study.area == 'rockefeller', 'Rockefeller', 
                           'Staten Island'),
            study_year = paste('Year', study.year),
            tel = map(tel, data.frame)) %>%
  unnest(tel) %>%
  # create a column of days since August 1st
  mutate(date = as.Date(timestamp),
         days_since_aug_1 =
           if_else(
             # if in Aug, Sept, Oct, Nov, or Dec
             month(date) >= 8,
             # then calculate days since Aug 1
             date - as.Date(paste0(year(date), '-08-01')),
             # otherwise calculate days since Aug 1 of previous year 
             date - as.Date(paste0(year(date) - 1, '-08-01'))) %>%
           as.numeric()) # convert difftime to numeric
d

d_dt <- d %>%
  group_by(animal, Sex, site, study_year) %>%
  summarize(dt = median(diff(timestamp), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(dt = as.numeric(dt) / (1 %#% 'hour'))

# make figure of doy and days since August 1st
cowplot::plot_grid(
  ggplot(d_dt, aes(animal, dt)) +
    facet_grid(study_year ~ .) +
    geom_vline(aes(xintercept = animal), alpha = 0.1,
               d %>%
                 mutate(animal_year = paste(animal, study_year)) %>%
                 group_by(animal_year) %>%
                 slice(1) %>%
                 ungroup() %>%
                 filter(duplicated(animal)) %>%
                 select(! study_year)) +
    geom_point(aes(shape = Sex, color = site)) +
    geom_text(aes(label = animal), filter(d_dt, dt > 2),
              nudge_y = -0.2, nudge_x = 0.5) +
    labs(x = NULL, y = 'Median sampling interval (hours)') +
    scale_color_brewer('Site', type = 'qual', palette = 1) +
    theme(axis.text.x.bottom = element_text(angle = 90, face = 'plain',
                                            size = 7),
          legend.position = 'top'),
  ggplot() +
    geom_histogram(aes(yday(date)), d, fill = 'grey', color = 'black',
                   center = 5, binwidth = 10) +
    labs(x = 'Day of year', y = 'Count'),
  ggplot() +
    geom_histogram(aes(days_since_aug_1), d, fill = 'grey', color = 'black',
                   center = 5, binwidth = 10) +
    labs(x = expression(bold(Days~since~August~1^bold(st))), y = 'Count') +
    xlim(c(0, 370)),
  ncol = 1, labels = 'AUTO', rel_heights = c(1.5, 1, 1))

ggsave('figures/dt-and-yday-days-since-august-first-hist.png',
       width = 12, height = 8, dpi = 600, bg = 'white')
