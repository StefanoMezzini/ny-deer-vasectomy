library('purrr') #' for functional programming (`map()`)
library('dplyr') #' for data wrangling (`mutate()`, `%>%`, ...)

setwd('../ny-deer/')
files <- list.files('data/movement-state-data', '.zip', full.names = TRUE)

map(files, \(.file) unzip(zipfile = .file, exdir = 'data/movement-state-data'))
beepr::beep()

files <- list.files('data/movement-state-data/Classified_states', '*',
                    full.names = FALSE, recursive = TRUE)

file.copy(paste0('data/movement-state-data/Classified_states/', files),
          files %>%
            gsub(pattern = '/[^/]+/', replacement = '/', perl = TRUE) %>%
            gsub(pattern = '/classified',
                 replacement = '/classified-states/classified') %>%
            paste0('data/', .),
          overwrite = TRUE)
beepr::beep()

# need to split by year
for(y in 1:2){
  files <- list.files(path = paste0('data/Year 1/classified-states/'),
                      pattern = '*', full.names = TRUE, recursive = TRUE)
  
  head(data.table::fread(files[1]))
  
  d <- map_dfr(files,
               .f = \(.f) data.table::fread(.f, verbose = FALSE,
                                            showProgress = FALSE) %>%
                 transmute(
                   timestamp,
                   n = as.integer(n),
                   animal = gsub(pattern = '.*classified_states_',
                                 replacement = '', x = .f,
                                 perl = TRUE) %>%
                     gsub(pattern = '.csv', replacement = '', x = .,
                          perl = TRUE),
                   study_year = gsub(pattern = '.*Year ',
                                     replacement = '',
                                     x = .f, perl = TRUE) %>%
                     gsub(pattern = '/.*', replacement = '', x = .,
                          perl = TRUE),
                   state, sex, study_site),
               .progress = TRUE)
  
  saveRDS(d, paste0('data/movement-states-year-', y, '.rds'))
  rm(d)
  gc()
}

d <- bind_rows(readRDS('data/movement-states-year-1.rds'),
               readRDS('data/movement-states-year-2.rds'))

saveRDS(d, 'data/movement-states-years-1-and-2.rds')

plot(factor(d$state, levels = c('high', 'medium', 'low', 'no activity')) ~
       d$timestamp, subset = d$animal == d$animal[1],
     col = c('grey', '#004488', '#DDAA33', '#BB5566'),
     xlab = '', ylab = 'Proportion for each state')
