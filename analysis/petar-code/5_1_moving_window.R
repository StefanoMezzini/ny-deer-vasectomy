library(ctmm)
library(tidyverse)
library(dplyr)
library(purrr)
library(ggplot2)
library(magrittr)

library(doParallel)
library(foreach)

source("analysis/petar-code/5_1_moving_window_functions.R")

# Read the datasets
# ------------------------------------------------------------------------------

data.y1 <- data.table::fread("data/Year_1_processed_data.csv",
                             stringsAsFactors = FALSE,
                             header = TRUE) %>%
  as.data.frame()
data.y2 <- data.table::fread("data/Year_2_processed_data.csv",
                             stringsAsFactors = FALSE,
                             header = TRUE) %>%
  as.data.frame()

head(data.y1)
head(data.y2)


# Read the calibration error model
# ------------------------------------------------------------------------------
calibration_model <- readRDS("Data/moving_window/calibration.error.model.inventa_no_class.rds")
calibration_model

rownames(calibration_model$UERE) <- 'gps'
rownames(calibration_model$DOF) <- 'gps'


# Path to manually flagged outliers
outliers.files.path <- "Data/moving_window/outliers/"


# Set variables/paths

window = 7 %#% 'day' # size of the moving window
dt = 3 %#% 'day' # step size for the moving window - slide
cores = 6


# Year 1
for(ind in 46:length(unique(data.y1$`individual-local-identifier`))){ # 1:length(unique(data.y1$`individual-local-identifier`))
  year.id <- 1
  individual.id <- unique(data.y1$`individual-local-identifier`)[ind]
  
  moving_window(individual.id = individual.id, 
                year.id = year.id,
                data.y1 = data.y1, 
                data.y2 = data.y2, 
                calibration_model, 
                outliers.files.path = outliers.files.path, 
                window = window, 
                dt = dt, 
                fig_path = "Data/moving_window/results/Year 1/", 
                cores = cores, 
                rds_path = "Data/moving_window/results/Year 1/")

}



# Year 2
for(ind in 58:length(unique(data.y2$`individual-local-identifier`))){
  year.id <- 2
  individual.id <- unique(data.y2$`individual-local-identifier`)[ind]
  
  moving_window(individual.id = individual.id, 
                year.id = year.id,
                data.y1 = data.y1, 
                data.y2 = data.y2, 
                calibration_model, 
                outliers.files.path = outliers.files.path, 
                window = window, 
                dt = dt, 
                fig_path = "Data/moving_window/results/Year 2/", 
                cores = cores, 
                rds_path = "Data/moving_window/results/Year 2/")
  
}
