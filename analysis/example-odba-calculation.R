library('dplyr') # for data wrangling

#' the there must be a space somewhere that causes `read.csv()` to see the
#' dataset as having too many columns to name them VX, where X is the
#' column number
odba_data <- read.csv('data/odba-data-2021-11-25.csv', sep = ',')

odba_data <- read.csv('data/odba-data-2021-11-25.csv',
                      skip = 1) # skip metadata row
colnames(odba_data)
colnames(odba_data) <- c('timestamp_utc', 'milliseconds', 'acc_x_g',
                         'acc_y_g', 'acc_z_g', 'temperature_c')
head(odba_data)

accel_mat <- with(odba_data, cbind(acc_x_g, acc_y_g, acc_z_g))

# I don't know what the values should be, so I'm just guessing...
tictoc::tic()
odba_data$odba <- tagtools::odba(A = accel_mat, # N by 3 matrix of values
                                 sampling_rate = 1,
                                 fh = 1 / 5) # sampling rate in Hz
tictoc::toc()
