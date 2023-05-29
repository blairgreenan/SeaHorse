# plot_SeaHorse_CTD2.R
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

# Toggle whether to load data from scratch (i.e., from the Matlab .mat file)
From_Scratch <- FALSE

if (From_Scratch){
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
    
    # Compute the mixed layer depth based on temperature and density criterion
    # first subset the profile to remove the upper and lower 10m of the water column
    # The following algorithm is based on D. Kelley (2018) p. 126-7
    Press_subset <- subset(Press, Press>10 & Press<190)
    Tmp_subset <- subset(Tmp, Press>10 & Press<190)
    SigmaT_subset <- subset(SigmaT, Press>10 & Press<190)
    criterion_T <- 0.1
    criterion_rho <- 0.125
    inMLD_T <- abs(Tmp_subset[1]-Tmp_subset) < criterion_T
    MLDindex_T <- which.min(inMLD_T)
    MLDpressure_T <- Press_subset[MLDindex_T]
    inMLD_rho <- abs(SigmaT_subset[1]-SigmaT_subset) < criterion_rho
    MLDindex_rho <- which.min(inMLD_rho)
    MLDpressure_rho <- Press_subset[MLDindex_rho]
    # there seem to be cases where the density criterion is not met through the whole water column
    if (MLDpressure_rho>150){
      MLDpressure_rho <- NA
    }
    
    if (ii == 1) {
      tidy_SH <- bind_cols(Tme_posixct,Cond,SigmaT,DescentRate,Fluor,O2PercentSat,O2Conc,PotTemp,Press,Dpth,Sal,SoundVel,Tmp,V0,V1,BF2,BF)
      tidy_MLD <- bind_cols(Tme_posixct,MLDpressure_T,MLDpressure_rho)
    } else {
      tidy_tmp <- bind_cols(Tme_posixct,Cond,SigmaT,DescentRate,Fluor,O2PercentSat,O2Conc,PotTemp,Press,Dpth,Sal,SoundVel,Tmp,V0,V1,BF2,BF)
      tidy_SH <- bind_rows(tidy_SH, tidy_tmp)
      tidy_MLD_tmp <- bind_cols(Tme_posixct,MLDpressure_T,MLDpressure_rho)
      tidy_MLD <- bind_rows(tidy_MLD, tidy_MLD_tmp)
    }
    
  }
  names(tidy_SH)[1] <- "Time"
  names(tidy_SH)[2] <- "Conductivity"
  names(tidy_SH)[3] <- "Density"
  names(tidy_SH)[4] <- "Descent_Rate"
  names(tidy_SH)[5] <- "Fluorescence"
  names(tidy_SH)[6] <- "O2_Percent_Sat"
  names(tidy_SH)[7] <- "O2_Concentration"
  names(tidy_SH)[8] <- "Potential_Temperature"
  names(tidy_SH)[9] <- "Pressure"
  names(tidy_SH)[10] <- "Depth"
  names(tidy_SH)[11] <- "Salinity"
  names(tidy_SH)[12] <- "Sound_Velocity"
  names(tidy_SH)[13] <- "Temperature"
  names(tidy_SH)[14] <- "V0"
  names(tidy_SH)[15] <- "V1"
  names(tidy_SH)[16] <- "Buoyancy_Freq_2"
  names(tidy_SH)[17] <- "Buoyancy_Freq"
  
  names(tidy_MLD)[1] <- "Time"
  names(tidy_MLD)[2] <- "Pressure_T"
  names(tidy_MLD)[3] <- "Pressure_rho"
  
  
  # Save the data to a RData file
  save(tidy_SH, tidy_MLD, file = "SeaHorse_CTD_data.RData")
  
} else {
  load("SeaHorse_CTD_data.RData")
}

# Estimate the mixed layer depth

###### Plotting Section ###############

# Temperature plot
# Set flag for plotting ML depth based on temperature criterion
# Note that the T criterion does not seem to work as well as the density
# criterion, so I am not planning to present this
Tflag <- FALSE

if (Tflag) {
  pTmp <- ggplot() +
    geom_tile(data=tidy_SH, aes(Time, Pressure, fill = Temperature)) +
    scale_fill_cmocean(name = "thermal") +
    scale_y_reverse() +
    labs(x=NULL,y="Depth (m)") +
    guides(fill = guide_colourbar(title=expression("Temperature (\u00B0C)"))) +
    scale_x_continuous(labels = NULL) +
    geom_smooth(data=tidy_MLD, aes(Time,Pressure_T), method="loess", se=FALSE, span=0.1, colour="white") +
    geom_point(data=tidy_MLD, aes(Time,Pressure_T))
  
} else {
  pTmp <- ggplot(tidy_SH, aes(Time, Pressure)) + 
    geom_tile(aes(fill = Temperature)) + 
    scale_fill_cmocean(name = "thermal") + 
    scale_y_reverse() + 
    labs(x=NULL,y="Depth (m)") + 
    guides(fill = guide_colourbar(title=expression("Temperature (\u00B0C)"))) + 
    scale_x_continuous(labels = NULL)
  
}


# Conductivity plot
pCond <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = Conductivity)) + 
  scale_fill_cmocean(name = "haline") + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Conductivity (S m"^"-1"*")"))) + 
  scale_x_continuous(labels = NULL)

# Density plot
pDensity <- ggplot() + 
  geom_tile(data=tidy_SH, aes(Time, Pressure, fill = Density)) + 
  scale_fill_cmocean(name = "dense") + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Density (kg m"^"-3"*")"))) + 
  scale_x_continuous(labels = NULL) +
  geom_smooth(data=tidy_MLD, aes(Time,Pressure_rho), method="loess", se=FALSE, span=0.1, colour="white") 
# could add dots for individual profiles but this clutters the plot
# +
#  geom_point(data=tidy_MLD, aes(Time,Pressure_rho))


# Descent Rate plot
pDescent_Rate <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = Descent_Rate)) + 
  scale_fill_cmocean(name = "speed", limits=c(-0.3,0)) + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Descent Rate (m s"^"-1"*")"))) + 
  scale_x_continuous(labels = NULL)

# Fluorescence plot
pFluorescence <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = Fluorescence)) + 
  scale_fill_cmocean(name = "algae") + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Chlorophyll (mg m"^"-3"*")"))) + 
  scale_x_continuous(labels = NULL)

# O2_Percent_Sat plot
pO2_Percent_Sat <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = O2_Percent_Sat)) + 
  scale_fill_cmocean(name = "dense", limits=c(65,105)) + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("O2_Percent_Sat")))

# O2_Concentration plot
pO2_Concentration <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = O2_Concentration)) + 
  scale_fill_cmocean(name = "dense", limits=c(5.5,8.5)) + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Oxygen (mL L"^"-1"*")")))

# Potential_Temperature plot
pPotential_Temperature <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = Potential_Temperature)) + 
  scale_fill_cmocean(name = "thermal") + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Potential_Temperature (\u00B0C)")))

# Salinity plot
pSalinity <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = Salinity)) + 
  scale_fill_cmocean(name = "haline") + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Salinity        "))) + 
  scale_x_continuous(labels = NULL)

# Sound_Velocity plot
pSound_Velocity <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = Sound_Velocity)) + 
  scale_fill_cmocean(name = "speed") + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Sound_Velocity (m s"^"-1"*")"))) + 
  scale_x_continuous(labels = NULL)

# Buoyancy_Freq_2 plot
pBuoyancy_Freq_2 <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = Buoyancy_Freq_2)) + 
  scale_fill_cmocean(name = "dense") + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Buoyancy_Freq_Sqrd (rad"^2*"s"^"-1"*")")))

# Buoyancy_Freq plot
pBuoyancy_Freq <- ggplot(tidy_SH, aes(Time, Pressure)) + 
  geom_tile(aes(fill = Buoyancy_Freq)) + 
  scale_fill_cmocean(name = "dense") + 
  scale_y_reverse() + 
  labs(x=NULL,y="Depth (m)") + 
  guides(fill = guide_colourbar(title=expression("Buoyancy_Freq (cycles h"^"-1"*")")))

# facet plot using patchwork package
pTmp/pSalinity/pDensity/pFluorescence/pO2_Concentration

# save plot
ggsave(filename = "SH_CTD.png", device = "png", scale = 1.5, width = 6, height = 10, units = "in", dpi = 1200)


