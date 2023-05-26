# plot_SeaHorse_CTD.R
# Blair Greenan
# Fisheries and Oceans Canada
# 24 May 2023
#
# Description: this script generates a filled contour plot with profile of CTD
# data collected by the SeaHorse mooring on Ross Bank. This version uses the 
# Matlab mat file the was processed by Susanne Craig post-cruise in 2012.
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

# load CTD Data from the .mat file
SeaHorse_CTD2 <- readMat("SeaHorse_Deployment_3.mat")

for (ii in 1:138) {
  Tme <- as.vector(SeaHorse_CTD2$ctd[[1+((ii-1)*20)]])
  Cond <- as.vector(SeaHorse_CTD2$ctd[[3+((ii-1)*20)]])
  SigmaT <- as.vector(SeaHorse_CTD2$ctd[[4+((ii-1)*20)]])
  DescentRate <- as.vector(SeaHorse_CTD2$ctd[[5+((ii-1)*20)]])
  Fluor <- as.vector(SeaHorse_CTD2$ctd[[6+((ii-1)*20)]])
  O2PercentSat <- as.vector(SeaHorse_CTD2$ctd[[7+((ii-1)*20)]])
  O2Conc <- as.vector(SeaHorse_CTD2$ctd[[8+((ii-1)*20)]])
  PotTemp <- as.vector(SeaHorse_CTD2$ctd[[9+((ii-1)*20)]])
  Press <- as.vector(SeaHorse_CTD2$ctd[[10+((ii-1)*20)]])
  Dpth <- as.vector(SeaHorse_CTD2$ctd[[11+((ii-1)*20)]])
  Sal <- as.vector(SeaHorse_CTD2$ctd[[12+((ii-1)*20)]])
  SoundVel <- as.vector(SeaHorse_CTD2$ctd[[13+((ii-1)*20)]])
  Tmp <- as.vector(SeaHorse_CTD2$ctd[[14+((ii-1)*20)]])
  V0 <- as.vector(SeaHorse_CTD2$ctd[[15+((ii-1)*20)]])
  V1 <- as.vector(SeaHorse_CTD2$ctd[[16+((ii-1)*20)]])
  BF2 <- as.vector(SeaHorse_CTD2$ctd[[17+((ii-1)*20)]])
  BF <- as.vector(SeaHorse_CTD2$ctd[[18+((ii-1)*20)]])
  Flag <- as.vector(SeaHorse_CTD2$ctd[[19+((ii-1)*20)]])
  SHUnits <- SeaHorse_CTD2$ctd[[20+((ii-1)*20)]]

  # Convert Matlab datenum to POSIX compliant date/time
  Tme_posixct <- as.POSIXct((as.numeric(Tme) - 719529)*86400, origin = "1970-01-01", tz = "UTC")
  
  if (ii == 1) {
    tidy_SH <- bind_cols(Tme_posixct,Cond,SigmaT,DescentRate,Fluor,O2PercentSat,O2Conc,PotTemp,Press,Dpth,Sal,SoundVel,Tmp,V0,V1,BF2,BF)
  } else {
  tidy_tmp <- bind_cols(Tme_posixct,Cond,SigmaT,DescentRate,Fluor,O2PercentSat,O2Conc,PotTemp,Press,Dpth,Sal,SoundVel,Tmp,V0,V1,BF2,BF)
  tidy_SH <- bind_rows(tidy_SH, tidy_tmp)
  }
  
}
names(tidy_SH)[1] <- "Time"
names(tidy_SH)[2] <- "Conductivity"
names(tidy_SH)[3] <- "Density"
names(tidy_SH)[4] <- "Descent Rate"
names(tidy_SH)[5] <- "Fluorescence"
names(tidy_SH)[6] <- "O2 Percent Sat"
names(tidy_SH)[7] <- "O2 Concentration"
names(tidy_SH)[8] <- "Potential Temperature"
names(tidy_SH)[9] <- "Pressure"
names(tidy_SH)[10] <- "Depth"
names(tidy_SH)[11] <- "Salinity"
names(tidy_SH)[12] <- "Sound Velocity"
names(tidy_SH)[13] <- "Temperature"
names(tidy_SH)[14] <- "V0"
names(tidy_SH)[15] <- "V1"
names(tidy_SH)[16] <- "Buoyancy Freq 2"
names(tidy_SH)[17] <- "Buoyancy Freq"

# Temperature plot
p1 <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = Temperature)) + 
  scale_fill_cmocean(name = "thermal") + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Temperature (\u00B0C)"))) + 
  scale_x_continuous(labels = NULL)

# Conductivity plot
p2 <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = Conductivity)) + 
  scale_fill_cmocean(name = "haline") + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Conductivity (S m"^"-1"*")"))) + 
  scale_x_continuous(labels = NULL)





# save(tidy_SH, file = "SeaHorse_CTD_data.RData")

########################
# SAVE DATA in a Rdat file
########################



