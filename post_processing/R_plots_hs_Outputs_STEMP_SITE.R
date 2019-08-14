#---------------------------------------------------------------------------------------
# Import observed temperature data and heat source hourly temperature output,
# compare to observed data (hourly and 7DADM), make plots, and generate model fit summary statistics
# Ryan Michie
#---------------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(hydroGOF)

sim_name <- "Yach_s1_2"

sim_file <- "Temp_H2O"
obs_file <- "ObservedData.csv"

out_dir <-"F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_2/outputs/"
sim_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_2/outputs/"
obs_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/"
fun_dir <- "E:/GitHub/Rscripts/heatsource"

source(paste0(fun_dir,"/","R_functions_hs_read_and_format.R"))

setwd(out_dir)

sim.hour <- read.hs.outputs(output_dir=sim_dir, file_name="Temp_H2O", 
                            constituent_name="Temperature (deg-C)",
                            statistic_name="hourly",
                            sim_name=sim_name, hs_ver=9, sheet_name = NULL)

sim.7dadm <- calc.7dadm(sim.hour)

obs.raw <- read.obs(obs_dir=obs_dir, file_name=obs_file)

obs <- obs.raw[obs.raw$constituentCode == "STEMPH",]

obs$hour <-as.integer(format(obs$Datetime, "%H"))

# Drop these cols
#obs$notes <- NULL
#obs$constituentCode <- NULL
#obs$siteName <-NULL

# data frame of obkm and simkm with site name
km <- obs2simkm(obs_dir=obs_dir, file_name="ObservedData.csv", 
                constituentCode="STEMPH",
                simkm=as.numeric(unique(sim.hour$Stream_km)))

# Drop these cols
obs$notes <- NULL
obs$constituentCode <- NULL
obs$constituent <- NULL
obs$statistic <- NULL
obs$siteName <-NULL

obs <- merge(km,obs,by.x="obs.km",by.y="Stream_km")

obs$Stream_km <- obs$sim.km
obs$obs.km <- NULL
obs$sim.km <- NULL

obs.wide  <- obs.wide <- obs %>%
  dplyr::select(Stream_km, siteName) %>%
  dplyr::distinct()

# calc 7DADM
obs$sim <- "Observations"
obs.7dadm <- calc.7dadm(df=obs)

# Only use observed data for the same time period as the model output
obs <- obs[obs$Datetime >= min(sim.hour$Datetime) & obs$Datetime <= max(sim.hour$Datetime),]

obs$siteName <- NULL
obs$Constituent <- NULL

sim.hour$statistic <- NULL
sim.hour$constituent <- NULL
sim.7dadm$constituent <- NULL

sim.hour$sim <-"Predictions"
sim.7dadm$sim <- "Predictions"

# vector of co-located model and obs stream km
obs.km <- as.numeric(unique(obs.wide$Stream_km))

# Add posixct date back into 7dadm
obs.7dadm$Date <- as.POSIXct(obs.7dadm$Date,format="%m/%d/%Y")
sim.7dadm$Date <- as.POSIXct(sim.7dadm$Date,format="%m/%d/%Y")

#---------------------------------------------------------------------------------------
# PLOT HOURLY SIMULATED VS. OBSERVED TEMPERATURES #

df.hour <- merge(sim.hour,obs, by=c("Date","Stream_km", "hour"))

# set y axis plot limits
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

ymin <- round_any(min(sim.hour$value, obs$value, na.rm=TRUE), 5, floor)
ymax <- round_any(max(sim.hour$value, obs$value, na.rm=TRUE), 5, ceiling)

i <-1

for(i in 1:length(obs.km)) {

  # The actual model stream kilometer
  skm <- obs.km[i]
  
  obs.i <- obs[obs$Stream_km == skm,]
  sim.i <- sim.hour[sim.hour$Stream_km == skm,]
  df.stats <- df.hour[df.hour$Stream_km == skm,]

  df.i <- rbind(sim.i,obs.i)
  
  obs.wide$hourly_me[i] <- round(me(df.stats$value.x, df.stats$value.y,na.rm = TRUE), digits = 2)
  obs.wide$hourly_mae[i] <- round(mae(df.stats$value.x, df.stats$value.y,na.rm = TRUE), digits = 2)
  obs.wide$hourly_rmse[i] <- round(rmse(df.stats$value.x, df.stats$value.y,na.rm = TRUE), digits = 2)
  obs.wide$hourly_ns[i] <- round(mNSE(df.stats$value.x, df.stats$value.y, j=2, na.rm = TRUE), digits = 2)
  obs.wide$hourly_n[i] <- nrow(na.omit(df.stats[,c("value.x","value.y")]))

  p1 <- ggplot(data=df.i, aes(x=Datetime, y=value, size=sim, 
                              linetype=sim, colour=sim)) +
    geom_line() +
    scale_color_manual(values=c("black", "red")) +
    scale_size_manual(values=c(1, 0.5)) +
    scale_linetype_manual(values=c("solid","solid")) +
    #scale_y_continuous(limits = c(5, 20)) +
    scale_y_continuous(limits = c(ymin, ymax)) +
    xlab("Date") +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y") +
    ylab("Hourly Temperature (C)") + 
    labs(title=unique(obs.wide$SiteName[i]),
         subtitle=paste0("Model Kilometer ",skm)) +
    theme(text=element_text(size=10)) +
    theme(legend.position="bottom") +
    theme(legend.title=element_blank())
  
  ggsave(file=paste0(out_dir,sim_name,"_HOURLY_km_",skm,".png"),
         plot=p1,
         height=3.5,
         width=6.75,
         units="in")
}

#---------------------------------------------------------------------------------------
# PLOT 7DADM SIMULATED VS. OBSERVED TEMPERATURES 

df.7dadm <- merge(sim.7dadm, obs.7dadm, by=c("Date","Stream_km","statistic"))

i <- 2

for(i in 1:length(obs.km)) {

  # The actual model stream kilometer
  skm <- obs.km[i]
  
  obs.i <- obs.7dadm[obs.7dadm$Stream_km == skm & obs.7dadm$statistic== "7DADM Temperature",]
  sim.i <- sim.7dadm[sim.7dadm$Stream_km == skm & sim.7dadm$statistic== "7DADM Temperature",] 
  df.stats <- df.7dadm[df.7dadm$Stream_km == skm & df.7dadm$statistic== "7DADM Temperature",] 

  df.i <- rbind(obs.i,sim.i)

  obs.wide$sdadm_me[i] <- round(me(df.stats$value.x, df.stats$value.y,na.rm = TRUE), digits = 2)
  obs.wide$sdadm_mae[i] <- round(mae(df.stats$value.x, df.stats$value.y,na.rm = TRUE), digits = 2)
  obs.wide$sdadm_rmse[i] <- round(rmse(df.stats$value.x, df.stats$value.y,na.rm = TRUE), digits = 2)
  obs.wide$sdadm_ns[i] <- round(mNSE(df.stats$value.x, df.stats$value.y, j=2, na.rm = TRUE), digits = 2)
  obs.wide$sdadm_n[i] <- nrow(na.omit(df.stats[,c("value.x","value.y")]))
  
  p1 <- ggplot(data=df.i, aes(x=Date, y=value, colour=sim, size=sim, linetype=sim)) +
    geom_line() +
    scale_color_manual(values=c("black", "red")) +
    scale_size_manual(values=c(1, 0.5)) +
    scale_linetype_manual(values=c("solid","solid")) +
    #scale_y_continuous(limits = c(5, 20)) +
    scale_y_continuous(limits = c(ymin, ymax)) +
    xlab("Date") +
    ylab("7DADM Temperature (C)") + 
    labs(title=unique(obs.wide$SiteName[i]),
         subtitle=paste0("Model Kilometer ",skm)) +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y") +
    theme(text=element_text(size=10)) +
    theme(legend.position="bottom") +
    theme(legend.title=element_blank())
  
  ggsave(file=paste0(out_dir,sim_name,"_7DADM_km_",skm,".png"),
         plot=p1,
         height=3.5,
         width=6.75,
         units="in")
}

#---------------------------------------------------------------------------------------
# STATS FOR ALL SITES

# add row for all sites
new.row <- data.frame("All Sites", -99, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
colnames(new.row) <-colnames(obs.wide)
obs.wide <- rbind(obs.wide,new.row)

i <- nrow(obs.wide)

obs.wide$hourly_me[i] <- round(me(df.hour$value.x, df.hour$value.y,na.rm = TRUE), digits = 2)
obs.wide$hourly_mae[i] <- round(mae(df.hour$value.x, df.hour$value.y,na.rm = TRUE), digits = 2)
obs.wide$hourly_rmse[i] <- round(rmse(df.hour$value.x, df.hour$value.y,na.rm = TRUE), digits = 2)
obs.wide$hourly_ns[i] <- round(mNSE(df.hour$value.x, df.hour$value.y, j=2, na.rm = TRUE), digits = 2)
obs.wide$hourly_n[i] <- nrow(na.omit(df.hour[,c("value.x","value.y")]))

obs.wide$sdadm_me[i] <- round(me(df.7dadm$value.x, df.7dadm$value.y,na.rm = TRUE), digits = 2)
obs.wide$sdadm_mae[i] <- round(mae(df.7dadm$value.x, df.7dadm$value.y,na.rm = TRUE), digits = 2)
obs.wide$sdadm_rmse[i] <- round(rmse(df.7dadm$value.x, df.7dadm$value.y,na.rm = TRUE), digits = 2)
obs.wide$sdadm_ns[i] <- round(mNSE(df.7dadm$value.x, df.7dadm$value.y, j=2, na.rm = TRUE), digits = 2)
obs.wide$sdadm_n[i] <- nrow(na.omit(df.7dadm[,c("value.x","value.y")]))

# make sure it is sorted
obs.wide <- obs.wide[with(obs.wide, order(-Stream_km)), ]

## write data to ouput
write.csv(obs.wide, file=paste(out_dir,sim_name,"_STATS.csv",sep=""),row.names = FALSE)

setwd(fun_dir)
