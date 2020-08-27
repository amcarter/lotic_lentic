########################################
# Helpers for running Lake Metabolizer on Stream Data

# a carter
# 2020-07-06
library(tidyverse)


format_data_for_lake_metabolizer <- function(dat){
  dat <- dat %>%
    select(datetime=DateTime_UTC, do.sat=DO.sat, doobs=DO_mgL, par=light, airt=AirTemp_C, wtr_tmp=WaterTemp_C, z.mix=level_m, discharge=Discharge_m3s)
  
  dat$datetime <- ymd_hms(dat$datetime)
  #dat$wnd_10 <- NA
  dat$day_night <- 1
  dat$day_night[dat$par<1]<- 0
  #dat$k.gas <- 0.5
  dat <- as.data.frame(dat)
  return(dat)
}



calc_raymond_K600_eq7 <- function(dat, site_deets, site){
  dat$S <- site_deets$slope[site_deets$sitecode==site]
  dat$V = dat$discharge/(dat$z.mix*site_deets$width_m[site_deets$sitecode==site])
  dat$D = dat$z.mix
  dat$Q = dat$discharge
  
  #get min and max for all terms in raymond eqn 7 (from table 2)
  dat$t1a = (dat$V * dat$S)^(0.86 + 0.016)
  dat$t1b = (dat$V * dat$S)^(0.86 - 0.016)
  dat$t2a = dat$Q^(-0.14 + 0.012)
  dat$t2b = dat$Q^(-0.14 - 0.012)
  dat$t3a = dat$D^(0.66 + 0.029)
  dat$t3b = dat$D^(0.66 - 0.029)
  dat$t1min = apply(dat[,c("t1a", "t1b")], 1, min)
  dat$t1max = apply(dat[,c("t1a", "t1b")], 1, max)
  dat$t2min = apply(dat[,c("t2a", "t2b")], 1, min)
  dat$t2max = apply(dat[,c("t2a", "t2b")], 1, max)
  dat$t3min = apply(dat[,c("t3a", "t3b")], 1, min)
  dat$t3max = apply(dat[,c("t3a", "t3b")], 1, max)
  
  #calculate k600 (m/d)
  dat$raymondk600min = (4725 - 445) * dat$t1min * dat$t2min * dat$t3min
  dat$raymondk600max = (4725 + 445) * dat$t1max * dat$t2max * dat$t3max
  dat$raymondk600mean = 4725 * (dat$V * dat$S)^(0.86) * dat$Q^(-0.14) * dat$D^(0.66)
  
  #convert to K600 (1/d)
  K600s = data.frame(
    apply(select(dat, starts_with('raymond')), 2, function(x) x / dat$D))
  colnames(K600s) = sub('k', 'K', colnames(K600s))
  dat = select(dat, datetime, raymondk600max,k600 =raymondk600mean, raymondk600min)
  return(dat)
}
