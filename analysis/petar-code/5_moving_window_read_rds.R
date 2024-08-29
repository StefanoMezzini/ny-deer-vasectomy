# Author: Petar Bursac (bursac.petar3@gmail.com)
# Date: 05.02.2024.

# Load major packages
library(ctmm)
library(tidyverse)
library(dplyr)
library(purrr)
library(ggplot2)
library(magrittr)

# Script and functions to open and read rds files as results from moving window analysis
# ------------------------------------------------------------------------------

# Function: make_graphics_mw
# Author: Petar Bursac (bursac.petar3@gmail.com)
# Date: 05.02.2024.
# Description: Function to make additional plots per id for mowing window results

# Parameters: 
# ind.id - individual id of interest
# year.id - study year of interest

make_graphics_mw <- function(ind.id, year.id) {
  
  ind.path <- paste0("Data/MoveBank_both_years/moving_window/Year ", year.id, "/", ind.id, "-window-7-days-dt-3-days.rds")
  data.mw <- readRDS(file = ind.path)
  
  # Example - plot the variogram and model fitted for each moving window - LARGE number of models
  
  # variogram and model fitted
  
  if(year.id == 2) {
    png(paste0("Data/MoveBank_both_years/moving_window/Year ", year.id, "/", ind.id, "-window-7-days-dt-3-days_var_mf.png"), width = 65, height = 65, units='cm', res = 600)
    
    par(mfrow = c(10, 10), mar = c(1,1,1,1))
  } else {
    png(paste0("Data/MoveBank_both_years/moving_window/Year ", year.id, "/", ind.id, "-window-7-days-dt-3-days_var_mf.png"), width = 60, height = 60, units='cm', res = 600)
    
    par(mfrow = c(9, 10), mar = c(1,1,1,1))
  }
  
  for(i in 1:dim(data.mw)[1]){
    if(length(data.mw[i, "model"][[1]][[1]]) == 1){
      plot.new()
    } else {
      plot(data.mw[i, "variogram"][[1]], data.mw[i, "model"][[1]])
    }
  }
  dev.off()
  par(mfrow = c(1, 1))
  
  
  # akde - home range
  
  if(year.id == 2) {
    png(paste0("Data/MoveBank_both_years/moving_window/Year ", year.id, "/", ind.id, "-window-7-days-dt-3-days_akdes.png"), width = 65, height = 65, units='cm', res = 600)
    
    par(mfrow = c(10, 10), mar = c(1,1,1,1))
  } else {
    png(paste0("Data/MoveBank_both_years/moving_window/Year ", year.id, "/", ind.id, "-window-7-days-dt-3-days_akdes.png"), width = 60, height = 60, units='cm', res = 600)
    
    par(mfrow = c(9, 10), mar = c(1,1,1,1))
  }
  
  for(i in 1:dim(data.mw)[1]){
    
    if(length(data.mw[i, "model"][[1]][[1]]) == 1) {
      plot.new()
    } else {
      plot(data.mw[i, "akde"][[1]])
    }
  }
  dev.off()
  par(mfrow = c(1, 1))
  
  
  # plot diffusion over time with CIs
  
  # plot results
  
  theme_set(theme_bw() + theme(legend.position = 'none'))
  
  plt_b <-
    ggplot(data.mw) +
    
    # 95% CIs for 95% diffusion estimates
    geom_ribbon(aes(date, ymin = diffusion_lwr, ymax = diffusion_upr), alpha = 0.3) +
    
    # 95% diffusion estimates
    geom_line(aes(date, diffusion_est), linewidth = 1.25) +
    geom_line(aes(date, diffusion_est, color = posixct)) +
    
    # scale_x_date(NULL, date_labels = '%b %Y') +
    scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 week",
                 date_labels = "%b %Y") +
    #scale_color_viridis_c() +
    labs(y = expression(Diffusion~(km^2/day))) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
  
  
  # Save figure as a png using the animal's name
  ggsave(filename = paste0("Data/MoveBank_both_years/moving_window/Year ", year.id, "/", ind.id, "-window-7-days-dt-3-days_diffusion.png"),
         plot = plt_b, units = 'in', width = 9, height = 3, dpi = 600, bg = 'white')
  
  print(paste0("Saved graphics for individual: ", ind.id))
  
  return("Done!")
}

# Example
make_graphics_mw(ind.id = "148", 
                 year.id = 1)



# Make graphics and export results for all individuals from Year 1
# ------------------------------------------------------------------------------

data.y1 <- data.table::fread("Data/MoveBank_both_years/Year_1_processed_data.csv",
                             stringsAsFactors = FALSE,
                             header = TRUE) %>%
  as.data.frame()

names.id.y1 <- unique(data.y1$`individual-local-identifier`)

for(i in 45:length(names.id.y1)){
  
  ind.id <- names.id.y1[i]
  
  print(paste0("Making graphics for individual: ", ind.id))
  
  make_graphics_mw(ind.id = ind.id, 
                   year.id = 1)
  
  
}


# Make graphics and export results for all individuals from Year 2
# ------------------------------------------------------------------------------

data.y2 <- data.table::fread("Data/MoveBank_both_years/Year_2_processed_data.csv",
                             stringsAsFactors = FALSE,
                             header = TRUE) %>%
  as.data.frame()

names.id.y2 <- unique(data.y2$`individual-local-identifier`)

for(i in 22:length(names.id.y2)){
  
  ind.id <- names.id.y2[i]
  
  print(paste0("Making graphics for individual: ", ind.id))
  
  make_graphics_mw(ind.id = ind.id, 
                   year.id = 2)
  
  print(paste0("Done: ", i, " / ", length(names.id.y2)))
}











