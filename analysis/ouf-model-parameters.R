library('dplyr')   # for data wrangling
library('ggplot2') # for fancy plots
library('ggExtra') # for marginal histograms
library('ctmm')    # for %#% operator
library('purrr')   # for
theme_set(theme_minimal())

d <- readxl::read_xlsx('data/OUF-parameters-year-1.xlsx') %>%
  filter(! is.na(ind.id)) %>% # last row has Na ID and only t_p non-NA
  mutate(
    # tau_p
    # clean up the units column
    t_position.units = substr(x = t_position.units,
                              start = nchar('τ[position] (x'),
                              stop = nchar(t_position.units) - 1),
    # convert all tau_p values to days
    across(c(t_position_low, t_position_mean, t_position_high),
           \(x) {imap_dbl(x, \(..x, ..i) {
             'days' %#% ..x %#% t_position.units[..i]
           })
           }),
    t_position.units = 'days', # correct units
    # tau_v
    # column is imported as character
    t_velocity_high = as.numeric(t_velocity_high),
    # clean up the units column
    t_velocity.units = substr(x = t_velocity.units,
                              start = nchar('τ[velocity] (x'),
                              stop = nchar(t_velocity.units) - 1),
    # convert all tau_p values to hours
    across(c(t_velocity_low, t_velocity_mean, t_velocity_high),
           \(x) {imap_dbl(x, \(..x, ..i) {
             'hours' %#% ..x %#% t_velocity.units[..i]
           })
           }),
    t_velocity.units = 'hours')

# tau_p
ggplot(d) +
  facet_wrap(~ period) +
  geom_density(aes(t_position_mean, color = sex, fill = sex),
               alpha = 0.3, adjust = 5) +
  labs(x = expression(tau[p]~(days)), y = 'Density') +
  scale_color_brewer('Sex', palette = 6, type = 'qual',
                     aesthetics = c('color', 'fill'))

d %>%
  group_by(sex, period, t_position.units) %>%
  summarise(t_p = median(t_position_mean, na.rm = TRUE))

# tau_v
ggplot(d) +
  facet_wrap(~ period) +
  geom_density(aes(t_velocity_mean, color = sex, fill = sex),
               alpha = 0.3, adjust = 2) +
  labs(x = expression(tau[v]~(hours)), y = 'Density') +
  scale_color_brewer('Sex', palette = 6, type = 'qual',
                     aesthetics = c('color', 'fill'))

d %>%
  group_by(sex, period, t_velocity.units) %>%
  summarise(t_v = median(t_velocity_mean, na.rm = TRUE))
