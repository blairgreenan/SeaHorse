# plot_SeaHorse_ADCP.R
# Blair Greenan
# Fisheries and Oceans Canada
# 5 June 2023
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

# Settings to play with the declination to see if we can get better alignment of the ADCPs
#decl5 <- 117.59
#decl6 <- 117.59

setwd("C:/Science Projects/Ross Sea/Documents/Papers/Ross Bank/Figures/Mooring/ADCP")
# Extract data from the downward looking ADCP on the SeaHorse mooring using the oce package
DPL5 <- read.oce(file="DPL5_000.000")
# Serial number 9184 is downward looking as per logsheet - the mix of upward and downward is because there are profiles before and after the mooring was deployed (i.e., on deck) that are included in the data file.
DPL5 <- oceSetMetadata(DPL5,'orientation','downward')
# Apply the magnetic declination adjustment
magField <- magneticField(-179.2531, -76.6601, 2012.06, version = 13) # 23-Jan-2012 = 2012.06
DPL5 <- applyMagneticDeclination(object = DPL5, declination = magField$declination, debug = 1)
#DPL5 <- applyMagneticDeclination(object = DPL5, declination = decl5, debug = 1)
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

# Extract data from the upward looking ADCP on the SeaHorse mooring using the oce package
DPL6 <- read.oce(file="DPL6_000.000")
# Serial number 14074 is upward looking as per logsheet - the mix of upward and downward is because there are profiles before and after the mooring was deployed (i.e., on deck) that are included in the data file.
DPL6<-oceSetMetadata(DPL6,'orientation','upward')
# Apply the magnetic declination adjustment
magField <- magneticField(-179.2531, -76.6601, 2012.06, version = 13) # 23-Jan-2012 = 2012.06
DPL6 <- applyMagneticDeclination(object = DPL6, declination = magField$declination, debug = 1)
#DPL6 <- applyMagneticDeclination(object = DPL6, declination = decl6, debug = 1)
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


# Plot tiled images of downward and upward looking ADCP results using oce defaults
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


# Create a merged data set with the upward and downward-looking ADCPs on the SeaHorse
# mooring. Given the physical distance between the two ADCPs on the mooring wire
# of about 200m, we will use Bins 1:24 from each ADCP and this will approximately 
# create a merged data set.

# Eastward velocity component
# flip the depth bin order of the top ADCP since it is downward-looking
rev_DPL5_east <- DPL5@data$v[,27:1,1]
# bind the columns of the two ADCP usng the first 24 bins of each
bind_east <- cbind(DPL6@data$v[,1:24,1],rev_DPL5_east[,4:27])
# We are going to ignore that the top ADCP sampled at 5 minutes to the hour and 
# the bottom ADCP sampled 5 after the hour to avoid side lobe interference with the
# SeaHorse. We will add 5 minutes to the top (DPL5) time vector
DPL_time <- DPL5_time + (5*60)
# there are some empty profiles from the ADCPs while they are on deck before and
# after deployment.  Clipping the data to remove these empty profiles.
bind_east_clipped <- bind_east[c(-1,-2,-3,-4,-5,-145,-146),]
DPL_time_clipped <- DPL_time[c(-1,-2,-3,-4,-5,-145,-146)]
# bind the time vector to the 2D velocity matrix
tidy_east_clipped <- bind_cols(DPL_time_clipped,bind_east_clipped)
# Add name for the time column
names(tidy_east_clipped)[1] <- "DateTime"
# Add names for the depth bins columns (which are just the actual depths)
names(tidy_east_clipped)[2:49] <- seq(from=197,to=9,by=-4)
# Use pivot_longer to turn this data into a tidy data set
tidy_east_long_clipped <- pivot_longer(data=tidy_east_clipped,cols=-DateTime,names_to="Depth",values_to="East")
tidy_east_long_clipped$Depth <- as.numeric(tidy_east_long_clipped$Depth)
# This results in the following: 6672 obs. of 3 variables
# head(tidy_east_long_clipped)
# A tibble: 6 × 3
#DateTime            Depth East_component
#<dttm>              <dbl>          <dbl>
# 1 2012-01-21 03:00:00   197         -0.172
# 2 2012-01-21 03:00:00   193         -0.171
# 3 2012-01-21 03:00:00   189         -0.174
# 4 2012-01-21 03:00:00   185         -0.172
# 5 2012-01-21 03:00:00   181         -0.172
# 6 2012-01-21 03:00:00   177         -0.170

# Plot the result
tidy_east_plot <- ggplot() + 
  geom_tile(data=tidy_east_long_clipped, aes(x=DateTime,y=Depth,fill=East), na.rm = TRUE) + 
  scale_fill_cmocean(name = "balance", limits=c(-0.65,0.65)) + 
  scale_y_reverse() + 
  xlab(NULL) + 
  ylab("Depth (m)")

# Northward velocity component
# flip the depth bin order of the top ADCP since it is downward-looking
rev_DPL5_north <- DPL5@data$v[,27:1,2]
# bind the columns of the two ADCP usng the first 24 bins of each
bind_north <- cbind(DPL6@data$v[,1:24,2],rev_DPL5_north[,4:27])
# We are going to ignore that the top ADCP sampled at 5 minutes to the hour and 
# the bottom ADCP sampled 5 after the hour to avoid side lobe interference with the
# SeaHorse. We will add 5 minutes to the top (DPL5) time vector
# DPL_time already done for east component
# DPL_time <- DPL5_time + (5*60)
# there are some empty profiles from the ADCPs while they are on deck before and
# after deployment.  Clipping the data to remove these empty profiles.
bind_north_clipped <- bind_north[c(-1,-2,-3,-4,-5,-145,-146),]
# DPL_time_clipped already done for east component
# DPL_time_clipped <- DPL_time[c(-1,-2,-3,-4,-5,-145,-146)]
# bind the time vector to the 2D velocity matrix
tidy_north_clipped <- bind_cols(DPL_time_clipped,bind_north_clipped)
# Add name for the time column
names(tidy_north_clipped)[1] <- "DateTime"
# Add names for the depth bins columns (which are just the actual depths)
names(tidy_north_clipped)[2:49] <- seq(from=197,to=9,by=-4)
# Use pivot_longer to turn this data into a tidy data set
tidy_north_long_clipped <- pivot_longer(data=tidy_north_clipped,cols=-DateTime,names_to="Depth",values_to="North")
tidy_north_long_clipped$Depth <- as.numeric(tidy_north_long_clipped$Depth)
# This results in the following: 6672 obs. of 3 variables
# head(tidy_north_long_clipped)
# A tibble: 6 × 3
#DateTime            Depth East_component
#<dttm>              <dbl>          <dbl>
# 1 2012-01-21 03:00:00   197 -0.348
# 2 2012-01-21 03:00:00   193 -0.353
# 3 2012-01-21 03:00:00   189 -0.363
# 4 2012-01-21 03:00:00   185 -0.369
# 5 2012-01-21 03:00:00   181 -0.363
# 6 2012-01-21 03:00:00   177 -0.377

# Plot the result
tidy_north_plot <- ggplot() + 
  geom_tile(data=tidy_north_long_clipped, aes(x=DateTime,y=Depth,fill=North), na.rm = TRUE) + 
  scale_fill_cmocean(name = "balance", limits=c(-0.65,0.65)) + 
  scale_y_reverse() + 
  xlab(NULL) + 
  ylab("Depth (m)")


##### Plot depth-average East and North components and compare to tidal model

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

# Use dplyr::summarise to compute the depth-averaged mean velocity components
tidy_east_summary <- tidy_east_long_clipped %>% group_by(DateTime) %>% summarise(avg = mean(East, na.rm = TRUE))
tidy_north_summary <- tidy_north_long_clipped %>% group_by(DateTime) %>% summarise(avg = mean(North, na.rm = TRUE))

# Add the ROMS model output at the SeaHorse mooring site
date_start <- ymd_hms("2012-01-21 04:00:00")  #UTC
ROMS_vel <- read_csv("ROMS_velocity.csv", col_names = c("time_inc","East","North"))
ROMS_vel$time_inc <- date_start + seconds(ROMS_vel$time_inc)


####### Plots of the depth-averaged velocity data ##################
compare_east <- ggplot() + 
  geom_line(data=tidy_east_summary, aes(DateTime,avg, color="red")) + 
  geom_line(data=tide_u_tibble, aes(Time, East, color="blue")) + 
  geom_line(data=ROMS_vel, aes(time_inc, East, color="black")) +
  labs(x=NULL,y=expression(paste("East (m ", s^-1, ")"))) +
  ylim(-0.5, 0.5) +
  scale_color_manual(name=NULL, labels=c("ROMS","Tidal Model","ADCP"), values = c("black", "red", "blue")) + 
  scale_x_datetime(limits = c(as.POSIXct("2012-01-21 00:00:00"),as.POSIXct("2012-01-27 00:00:00"))) +
  theme(axis.text.x = element_blank())  # remove labels on tick marks

compare_north <- ggplot() + 
  geom_line(data=tidy_north_summary, aes(DateTime,avg, color="red")) + 
  geom_line(data=tide_v_tibble, aes(Time, North, color="blue")) + 
  geom_line(data=ROMS_vel, aes(time_inc, North, color="black")) + 
  labs(x=NULL,y=expression(paste("North (m ", s^-1, ")"))) +
  ylim(-0.5, 0.5) +
  scale_color_manual(name=NULL, labels=c("ROMS","Tidal Model","ADCP"), values = c("black", "red", "blue")) + 
  scale_x_datetime(limits = c(as.POSIXct("2012-01-21 00:00:00"),as.POSIXct("2012-01-27 00:00:00")))

# Use patchwork to plot results
dev.new()
compare_east/compare_north
ggsave(filename = "SH_ADCP_ROMS_17Jul2024.png", device = "png", scale = 1, width = 6, height = 3, units = "in", dpi = 1200)
dev.off()

############# Need to add a section on vertical shear

# Calculate the 4-m bin difference in the east and north velocity components
# Results in an matrix with one less row
# Need to do a transpose of the matrix for the diff function to work on the depth bins
# The transport after the diff operation to return the matrix to its orignal form
# Divide by 4 because that is the bin depth
bind_east_shear <- t(diff(t(bind_east_clipped)))/4
bind_north_shear <- t(diff(t(bind_north_clipped)))/4
shear <- sqrt(bind_east_shear^2 + bind_north_shear^2)
# Remove the bins where the two ADCPs overlap because these values distort the data
shear[,23:25] <- NA
# Just checking to see the range of the shear values
max(shear, na.rm=TRUE)
min(shear, na.rm=TRUE)
# Plot the results
# This does not provide any pattern for vertical shear so I think we should
# not pursue this any further. The upward-looking ADCP is Likely to0 far off
# the seabed to measure shear in the bottom boundary layer
image(shear, zlim = c(0,0.01))

