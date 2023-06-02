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

# Estimate the depth-averaged currents
DPL_east <- (DPL5_tibble$East+DPL6_tibble$East)/2
DPL_north <- (DPL5_tibble$North+DPL6_tibble$North)/2

# Tidal model reults
tide_u <- readMat("ross_bank_top_u.mat")
tide_v <- readMat("ross_bank_top_v.mat")
Tme_posixct_u <- as.POSIXct((as.numeric(tide_u$SerialDay) - 719529)*86400, origin = "1970-01-01", tz = "UTC")
tide_u_df <- cbind(Tme_posixct_u, tide_u$TimeSeries[1,])
tide_u_tibble <- as_tibble(tide_u_df)
tide_u_tibble <- rename(tide_u_tibble,Time=1,East=2)
tide_u_tibble$Time <- as.POSIXct(tide_u_tibble$Time)
tide_u_tibble$East <- (tide_u_tibble$East)/100  # convert from cm/s to m/s
model_east <- ggplot(tide_u_tibble, aes(Time, East)) + 
  geom_line() + 
  scale_x_datetime(limits = c(as.POSIXct("2012-01-21 00:00:00"),as.POSIXct("2012-01-27 00:00:00")))
# model_east
Tme_posixct_v <- as.POSIXct((as.numeric(tide_v$SerialDay) - 719529)*86400, origin = "1970-01-01", tz = "UTC")
tide_v_df <- cbind(Tme_posixct_v, tide_v$TimeSeries[1,])
tide_v_tibble <- as_tibble(tide_v_df)
tide_v_tibble <- rename(tide_v_tibble,Time=1,North=2)
tide_v_tibble$Time <- as.POSIXct(tide_v_tibble$Time)
tide_v_tibble$North <- (tide_v_tibble$North)/100  # convert from cm/s to m/s
model_north <- ggplot(tide_v_tibble, aes(Time, North)) + 
  geom_line() + 
  scale_x_datetime(limits = c(as.POSIXct("2012-01-21 00:00:00"),as.POSIXct("2012-01-27 00:00:00")))
# model_north

####### Plots of the data ##################
dev.new()
ggplot() + geom_line(data=DPL5_tibble, aes(Time, East, color="red")) + 
  geom_line(data=DPL6_tibble, aes(Time, East, color="blue")) +
  labs(x=NULL,y="East (m/s)") + 
  scale_color_discrete(name=NULL, labels=c("Lower","Upper"))
# save plot
ggsave(filename = "ADCP_East.png", device = "png", scale = 1.5, width = 6, height = 10, units = "in", dpi = 1200)


dev.new()
ggplot() + geom_line(data=DPL5_tibble, aes(Time, North, color="red")) + 
  geom_line(data=DPL6_tibble, aes(Time, North, color="blue")) + 
  labs(x=NULL,y="North (m/s)") + 
  scale_color_discrete(name=NULL, labels=c("Lower","Upper"))
# save plot
ggsave(filename = "ADCP_NOrth.png", device = "png", scale = 1.5, width = 6, height = 10, units = "in", dpi = 1200)

dev.new()
ggplot() + geom_line(data=DPL_east, aes(Time, East, color="red")) + 
  geom_line(data=tide_u_tibble, aes(Time, East, color="blue")) +
  labs(x=NULL,y="East (m/s)") + 
  scale_color_discrete(name=NULL, labels=c("Model","ADCP")) +
  scale_x_datetime(limits = c(as.POSIXct("2012-01-21 00:00:00"),as.POSIXct("2012-01-27 00:00:00")))

dev.new()
ggplot() + geom_line(data=DPL_north, aes(Time, East, color="red")) + 
  geom_line(data=tide_u_tibble, aes(Time, East, color="blue")) +
  labs(x=NULL,y="East (m/s)") + 
  scale_color_discrete(name=NULL, labels=c("Model","ADCP")) +
  scale_x_datetime(limits = c(as.POSIXct("2012-01-21 00:00:00"),as.POSIXct("2012-01-27 00:00:00")))

# Plot tiled images of downward and upward looking ADCP results
# East component
dev.new()
par(mfrow=c(2,1))
plot(DPL5,which=1,zlim=c(-0.6,0.6))
plot(DPL6,which=1,zlim=c(-0.6,0.6))
# North component
dev.new()
par(mfrow=c(2,1))
plot(DPL5,which=2,zlim=c(-0.6,0.6))
plot(DPL6,which=2,zlim=c(-0.6,0.6))
# Vertical component
dev.new()
par(mfrow=c(2,1))
plot(DPL5,which=3,zlim=c(-0.1,0.1))
plot(DPL6,which=3,zlim=c(-0.1,0.1))






