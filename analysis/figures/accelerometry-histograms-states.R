library('purrr')   # for functional programming
library('dplyr')   # for data wrangling
library('ggplot2') # for fancy plots
source('analysis/figures/default-theme.R')

pal <- c('no activity' = 'grey',
         'low activity' = '#004488',
         'medium activity' = '#DDAA33',
         'high activity' = '#BB5566')

d <-
  map_dfr(list.files('data/18-vedba-data/', pattern = 'acc_data_18_.*.csv',
                     full.names = TRUE),
          \(fn) {
            data.table::fread(fn) %>%
              transmute(VeDBA,
                        log_VeDBA,
                        month = format(timestamp, '%b %Y') %>%
                          as.character()) %>%
              return()
          }) %>%
  data.frame() %>%
  as_tibble() %>%
  mutate(State = case_when(VeDBA == 0 ~ 'no activity',
                           log(VeDBA) < 1.31 ~ 'low activity',
                           log(VeDBA) < 4.17 ~ 'medium activity',
                           TRUE ~ 'high activity') %>%
           factor(levels = c('no activity',
                             'low activity',
                             'medium activity',
                             'high activity')),
         # place -Inf log(VeDBA) to min - 3 rather than - Inf
         log_VeDBA = if_else(VeDBA == 0,
                             log(min(VeDBA[which(VeDBA > 0)])) - 3,
                             log_VeDBA),
         month = factor(month, levels = unique(month)))

min_log_VeDBA <- log(min(d$VeDBA[which(d$VeDBA > 0)]))

p <-
  d %>%
  ggplot(aes(log_VeDBA, fill = State)) +
  geom_histogram(bins = 300) +
  geom_vline(xintercept = 1.31, lty = 'dashed') +
  geom_vline(xintercept = 4.17, lty = 'dashed') +
  labs(x = expression(bold(log(VedBA))), y = 'Count (square-root scale)') +
  scale_fill_manual(values = pal) +
  scale_y_continuous(transform = 'sqrt') +
  theme(legend.position = 'top'); p

ggsave('figures/vedba-histogram.png',
       plot = p, width = 6, height = 4, units = 'in', dpi = 600,
       bg = 'white')

p_month <-
  p +
  facet_wrap(~ month) +
  theme(legend.position = 'top')
ggsave('figures/vedba-histogram-months.png',
       plot = p_month, width = 6, height = 4, units = 'in', dpi = 600,
       bg = 'white', scale = 1.5)
