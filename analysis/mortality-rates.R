library('dplyr')

# rates for both years at once
tibble(sex = c('all', 'male', 'female', 'all', 'male', 'female'),
       group = c('treatment', 'treatment', 'treatment',
                 'control', 'control', 'control'),
       deaths = c(5, 3, 2, 8, 2, 6),
       n = c(64, 21, 43, 60, 19, 41), # collared twice appear twice
       perc = deaths / n * 100, # percent deaths
       perc_se = sqrt(perc * (100 - perc) / n)) %>%
  mutate(perc = round(perc, 1),
         perc_se = round(perc_se, 1)) %>%
  as.data.frame() # to avoid rounding

# Fisher's exact test ----
# overall
matrix(c(5, 64 - 5,
         8, 60 - 8), ncol=2) %>%
  fisher.test()

# males
matrix(c(3, 21 - 3,
         2, 19 - 2), ncol=2) %>%
  fisher.test()

# females
matrix(c(2, 43 - 2,
         6, 41 - 6), ncol=2) %>%
  fisher.test()
