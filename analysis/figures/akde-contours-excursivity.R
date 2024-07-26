library('ctmm')   # for movement models
library('dplyr')  # for data wrangling
library('tidyr')  # for data wrangling
library('purrr')  # for functional programming
library('sf')     # for spatial data
library('terra')  # to work with rasters
library('ggplot2') # for fancy figures
source('analysis/figures/default-theme.R')

tel <- d <- readRDS('models/full-telemetry-movement-models-2024-04-20.rds') %>%
  pull(tel) %>%
  first() %>%
  data.frame() %>%
  mutate(date = as.Date(timestamp)) %>%
  select(x, y, date)

a <- readRDS('models/full-telemetry-movement-models-2024-04-20.rds') %>%
  pull(akde) %>%
  first() %>%
  raster(DF = 'CDF') %>%
  as.data.frame(xy = TRUE)

ggplot() +
  geom_raster(aes(x / 1e3, y / 1e3, fill = layer), a) +
  geom_point(aes(x / 1e3, y / 1e3), tel, alpha = 0.5, pch = '.') +
  geom_path(aes(x / 1e3, y / 1e3), tel, alpha = 0.02) +
  geom_contour(aes(x / 1e3, y / 1e3, z = layer), a, color = 'darkorange',
               breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous('x (km)', expand = c(0, 0)) +
  scale_y_continuous('y (km)', expand = c(0, 0)) +
  scale_fill_gradient('Excursivity', low = 'blue', high = 'white',
                      limits = c(0, 1), aesthetics = c('color', 'fill')) +
  theme(legend.position = 'inside', legend.position.inside = c(0.15, 0.8))

ggsave('figures/ud-density-example.png', width = 7, height = 5, units = 'in',
       dpi = 600, bg = 'white')
