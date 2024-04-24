library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
library('mgcv')      # for modeling
source('../ndvi-stochasticity/functions/betals.r')

# polys <- map(x$akde[-c(17, 47, 61, 73, 85)],
#              \(.x) SpatialPolygonsDataFrame.UD(.x) %>%
#                st_as_sf() %>%
#                st_transform('EPSG:4326') %>%
#                st_union() %>%
#                st_as_sf()) %>%
#   bind_rows()

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
  mutate(animal = factor(animal)) %>%
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

hist(d$density)

hist(d$hr_est_95)
quantile(d$hr_est_95, probs = c(0.5, 0.95, 0.97, 0.98, 0.99, 1))

d %>%
  group_by(animal) %>%
  mutate(density = density / mean(density)) %>%
  ungroup() %>%
  ggplot() +
  facet_wrap(sex ~ study_site) +
  geom_smooth(aes(days_since_aug_1, density))

# 40 minutes with s(by=s_t) + s(s_t, year), s(animal)

d$animal_year <- factor(paste(d$animal, d$study_year))

# takes ~ 10 seconds
# study_year smooth widens CIs substantially
m <- bam(
  density ~
    s(days_since_aug_1, by = sex_treatment, k = 10) +
    # s(days_since_aug_1, study_year, by = sex_treatment, bs = 'fs', k = 10) +
    s(animal, bs = 're'),
  family = betar(link = 'logit'), # since CDF is in (0, 1)
  data = d,
  method = 'fREML',
  discrete = TRUE,
  control = gam.control(trace = TRUE))

# m <- bam(
#   density ~
#     s(days_since_aug_1, by = sex_treatment, k = 10) +
#     s(days_since_aug_1, animal_year, bs = 'fs', k = 10),
#   family = betar(link = 'logit'), # since CDF is in (0, 1)
#   data = d,
#   method = 'fREML',
#   discrete = TRUE,
#   control = gam.control(trace = TRUE))

summary(m)
plot(m, pages = 1, scheme = 1)

layout(matrix(1:4, ncol = 2))
gam.check(m, type = 'pearson')
abline(a = 0, b = 1, col = 'red3')
layout(1)

# all four groups have a long right tail
ggplot(mutate(d, e = resid(m, type = 'pearson'))) +
  facet_grid(sex ~ study_site, scales = 'free') +
  geom_density(aes(e), fill = 'grey') +
  geom_vline(xintercept = 0, color = 'grey20')

# no problematic trends in the residuals over time
ggplot(mutate(d, e = resid(m, type = 'pearson'))) +
  facet_grid(sex ~ study_site, scales = 'free') +
  geom_point(aes(days_since_aug_1, e)) +
  geom_smooth(aes(days_since_aug_1, e))
