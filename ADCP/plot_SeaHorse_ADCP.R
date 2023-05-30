# plot_SeaHorse_ADCP.R
# Blair Greenan
# Fisheries and Oceans Canada
# 30 May 2023
#
# Description: this script generates a plot of the ADCP
# data collected by the SeaHorse mooring on Ross Bank. 
#
# load libraries
library(oce)
library(R.matlab)
library(R.oo)
library(tidyverse)
library(lubridate)
library(cmocean)
library(interp)
library(patchwork)

# Magnetic Declination
# Date 23 Jan 2012
# Latitude 76.6601S
# Longitude 179.2532 W
# https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml
# Declination 117.5904E 

# Extract data from the downward looking ADCP on the SeaHorse mooring
DPL5 <- read.oce(file="DPL5_000.000")
# Serial number 9184 is downward looking as per logsheet - the mix of upward and downward is because there are profiles before and after the mooring was deployed (i.e., on deck) that are included in the data file.
DPL5<-oceSetMetadata(DPL5,'orientation','downward')
# Apply the magnetic declination adjustment
magField <- magneticField(-179.2531, -76.6601, 2012.06, version = 13) # 23-Jan-2012 = 2012.06
DPL5 <- applyMagneticDeclination(object = DPL5, declination = magField$declination, debug = 1)
# Extract the eastward velocity component
DPL5_east <- DPL5@data$v[,,1]
DPL5_tibble_east <- as_tibble(DPL5_east)
# Extract the northward velocity component
DPL5_north <- DPL5@data$v[,,2]
DPL5_tibble_north <- as_tibble(DPL5_north)
# Compute the depth-averaged mean for each profile
DPL5_rmean_east <- rowMeans(DPL5_tibble_east, na.rm = TRUE)
DPL5_rmean_north <- rowMeans(DPL5_tibble_north, na.rm = TRUE)
# Extract the time vector
DPL5_time <- DPL5@data$time
DPL5_tibble <- bind_cols(DPL5_time, DPL5_rmean_east, DPL5_rmean_north)
DPL5_tibble <- rename(DPL5_tibble,Time=1,East=2,North=3)

# Extract data from the upward looking ADCP on the SeaHorse mooring
DPL6 <- read.oce(file="DPL6_000.000")
# Serial number 14074 is downward looking as per logsheet - the mix of upward and downward is because there are profiles before and after the mooring was deployed (i.e., on deck) that are included in the data file.
DPL6<-oceSetMetadata(DPL6,'orientation','upward')
# Apply the magnetic declination adjustment
magField <- magneticField(-179.2531, -76.6601, 2012.06, version = 13) # 23-Jan-2012 = 2012.06
DPL6 <- applyMagneticDeclination(object = DPL6, declination = magField$declination, debug = 1)
# Extract the eastward velocity component
DPL6_east <- DPL6@data$v[,,1]
DPL6_tibble_east <- as_tibble(DPL6_east)
# Extract the northward velocity component
DPL6_north <- DPL6@data$v[,,2]
DPL6_tibble_north <- as_tibble(DPL6_north)
# Compute the depth-averaged mean for each profile
DPL6_rmean_east <- rowMeans(DPL6_tibble_east, na.rm = TRUE)
DPL6_rmean_north <- rowMeans(DPL6_tibble_north, na.rm = TRUE)
# Extract the time vector
DPL6_time <- DPL6@data$time
DPL6_tibble <- bind_cols(DPL6_time, DPL6_rmean_east, DPL6_rmean_north)
DPL6_tibble <- rename(DPL6_tibble,Time=1,East=2,North=3)


dev.new()
ggplot() + geom_line(data=DPL5_tibble, aes(Time, East, color="red")) + 
  geom_line(data=DPL6_tibble, aes(Time, East, color="blue")) +
  labs(x=NULL,y="East (m/s)") + 
  scale_color_discrete(name=NULL, labels=c("Upper","Lower"))
# save plot
ggsave(filename = "ADCP_East.png", device = "png", scale = 1.5, width = 6, height = 10, units = "in", dpi = 1200)


dev.new()
ggplot() + geom_line(data=DPL5_tibble, aes(Time, North, color="red")) + 
  geom_line(data=DPL6_tibble, aes(Time, North, color="blue")) + 
  labs(x=NULL,y="North (m/s)") + 
  scale_color_discrete(name=NULL, labels=c("Upper","Lower"))
# save plot
ggsave(filename = "ADCP_NOrth.png", device = "png", scale = 1.5, width = 6, height = 10, units = "in", dpi = 1200)







