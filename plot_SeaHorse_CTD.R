# plot_SeaHorse_CTD.R
# Blair Greenan
# Fisheries and Oceans Canada
# 24 May 2023
#
# Description: this script generates a filled contour plot with profile of CTD
# data collected by the SeaHorse mooring on Ross Bank
#
# load libraries
library(oce)
library(R.matlab)
library(R.oo)
library(tidyverse)
library(lubridate)
library(cmocean)

# load CTD Data from the .mat file
SeaHorse_CTD <- readMat("sh_ctd_deployment3.mat")

# Replace NAN with NA which represents Not Available in R
nanTime <- is.nan(SeaHorse_CTD$Time)
SeaHorse_CTD$Time[nanTime] <- NA
# Create a time vector
Tme <- SeaHorse_CTD$Time
Tme_vec <- Tme[1,]

# Extract the temperature matrix
T <- SeaHorse_CTD$Temp
# Replace NAN with NA which represents Not Available in R
nanTemp <- is.nan(T)
T[nanTemp] <- NA

# Remove the first column which doesn't seem to have any data
T <- T[,2:139]
# Bind a vector of 1:250 in the first column to represent the depth bins
T <- cbind(c(1:250), T)

# Convert to a tibble
T <- as_tibble(T)

# Set the names of the columns to be Depth and the profile time
names(T)[1] <- "Depth"
names(T)[2:139] <- as.character(Tme_vec[2:139])

T2 <- pivot_longer(T, cols = 2:139, names_to = "Time", values_to = "Temperature")
# Convert to Matlab datenum time to R date/time
T2$Time <- as.POSIXct((as.numeric(T2$Time) - 719529)*86400, origin = "1970-01-01", tz = "UTC")

# Plot the data
dev.new()
dev.new()
ggplot(T2) + geom_contour_filled(aes(x=Time, y=Depth, z=Temperature), na.rm = TRUE, bins = 10) + scale_y_reverse() + scale_fill_cmocean(name="thermal", discrete = TRUE) + labs(x=NULL, y="Depth (m)")




