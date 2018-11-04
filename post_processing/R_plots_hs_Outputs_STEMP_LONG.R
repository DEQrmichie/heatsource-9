#---------------------------------------------------------------------------------------
# Plot heat source output hourly and 7DADM stream temperatures and observed temps longitudinally by hour
# for selected dates. Plot output model period min, max, and mean 7DADM stream temps longitudnally by stream km,
#---------------------------------------------------------------------------------------

library (reshape2)
library (ggplot2)
library(plyr)

sim_name <- "Yach_s1_2"

plot.dates <- c("07/01/2005","07/15/2005","08/01/2005","08/15/2005","09/01/2005")

out_dir <-"F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_2/outputs/"
sim_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_2/outputs/"
obs_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/"
fun_dir <- "E:/GitHub/Rscripts/heatsource"

source(paste0(fun_dir,"/","R_functions_hs_read_and_format.R"))

setwd(out_dir)

#--  Read temps ------------------------------------
sim.hour <- read.hs.outputs(output_dir=out_dir, file_name="Temp_H2O", 
                          constituent_name="Hourly Stream Temperature", 
                          sim_name=sim_name, hs8=FALSE)

sim.7dadm <- calc.7dadm(sim.hour)

# calc summary 7DADM
sim.7dadm.summ <- calc.summary(sim.7dadm)

# Read observation data
obs <- read.obs(obs_dir=obs_dir, file_name="ObservedData.csv")
obs <- obs[obs$ConstituentCode == "STEMPH",]
obs$hour <-as.integer(format(obs$Datetime, "%H"))

# calc obs 7DADM
obs$sim <- "Observations"
obs.7dadm <- calc.7dadm(df=obs)

# Just get data on the plot date
sim.hour <- sim.hour[sim.hour$Date %in% plot.dates,]
obs <- obs[obs$Date %in% plot.dates,]

sim.7dadm<- sim.7dadm[sim.7dadm$Date %in% plot.dates,]
obs.7dadm <- obs.7dadm[obs.7dadm$Date %in% plot.dates,]

# Drop these cols
obs$Notes <- NULL
obs$ConstituentCode <- NULL
obs$Constituent <- NULL
obs$SiteName <-NULL
obs$Datetime <- NULL

# data frame of obkm and simkm with site name
km <- obs2simkm(obs_dir=obs_dir, file_name="ObservedData.csv", 
                ConstituentCode="STEMPH",
                simkm=as.numeric(unique(sim.hour$Stream_km)))

obs <- merge(km,obs,by.x="obs.km",by.y="Stream_km")

obs$Stream_km <- obs$sim.km
obs$obs.km <- NULL
obs$sim.km <- NULL
obs$SiteName <- NULL

df <- merge(sim.hour,obs, by=c("Date","hour","Stream_km"),all.x=TRUE)

colnames(df) <- c("Date","hour","Stream_km","Datetime","sim","constituent","pred","obs")
df$Date <- format(df$Datetime, "%m/%d/%Y")
df$hour <- format(df$Datetime, "%H:%M")

obs.7dadm$sim <- NULL
df.7dadm <- merge(sim.7dadm,obs.7dadm, by=c("Date","Stream_km","constituent"),all.x=TRUE)
colnames(df.7dadm) <- c("Date","Stream_km","constituent","sim","pred","obs")

#---------------------------------------------------------------------------------------
# PLOT HOURLY SIMULATED VS. OBSERVED TEMPERATURES #

i<-1

# set y axis plot limits
ymin <- round_any(min(df$pred, df$obs, na.rm=TRUE), 5, floor)
ymax <- round_any(max(df$pred, df$obs, na.rm=TRUE), 5, ceiling)

for (i in 1:length(plot.dates)) {
  
  df.i <- df[df$Date == plot.dates[i],]
  
  p.long <- ggplot(data=df.i, aes(x=Stream_km)) +
    geom_line(aes(y=pred, color="Predictions")) +
    geom_point(aes(y=obs, color="Observations"),size=.75) +
    guides(color=guide_legend(override.aes=list(shape=c(16,NA),linetype=c(0,1))))  +
    scale_colour_manual(values=c("Predictions"="red","Observations"="Black")) +
    theme(legend.position="bottom") +
    theme(legend.title=element_blank()) +
    labs(subtitle=paste0(unique(df.i$Date))) +
    xlab("Model Stream Kilometer") +
    ylab("Temperature (C)") +
    #ylim(5,20) +
    ylim(ymin, ymax) +
    facet_wrap(~hour,nrow = 6)
  
  ggsave(file=paste0(out_dir,sim_name,"_LONG_",gsub("/","_", plot.dates[i]),".png"),
         plot=p.long,
         height=7.2,
         width=7,
         units="in")
}


#---------------------------------------------------------------------------------------
# PLOT 7DADM SIMULATED VS. OBSERVED TEMPERATURES

df.i <- df.7dadm[df.7dadm$constituent== "7DADM Temperature",]

p.long2 <- ggplot(data=df.i, aes(x=Stream_km)) +
    geom_line(aes(y=pred, color="Predictions")) +
    geom_point(aes(y=obs, color="Observations"),size=.75) +
    guides(color=guide_legend(override.aes=list(shape=c(16,NA),linetype=c(0,1))))  +
    scale_colour_manual(values=c("Predictions"="red","Observations"="Black")) +
    theme(legend.position="bottom") +
    theme(legend.title=element_blank()) +
    #labs(subtitle=sim_name) +
    xlab("Model Stream Kilometer") +
    ylab("Temperature (C)") +
    #ylim(5,20) +
    ylim(ymin, ymax) +
    facet_wrap(~Date, nrow = length(plot.dates))
  
ggsave(file=paste0(out_dir,sim_name,"_LONG_7DADM_",gsub("/","_", plot.dates[i]),".png"),
         plot=p.long2,
         height=7.2,
         width=7,
         units="in")

#---------------------------------------------------------------------------------------
# PLOT MIN, MEAN, and MAX 7DADM SIMULATED TEMPERATURES FOR THE ENTIRE MODEL PERIOD

df.summ <- sim.7dadm.summ[sim.7dadm.summ$constituent== "7DADM Temperature",]

p.summ <- ggplot(data=df.summ, aes(x=Stream_km)) +
  geom_ribbon(aes(ymax = max, ymin = min, fill="Min & Max Range"), alpha = 0.6) + 
  geom_line(aes(y=mean, linetype="Mean")) +
  scale_linetype_manual(values=c("Mean"="solid")) +
  scale_fill_manual(values="skyblue") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("7DADM Temperature (C)") + 
  ylim(5,21)

ggsave(file=paste0(out_dir,sim_name,"_LONG_7DADM_SUMM.png"),
       plot=p.summ,
       height=3,
       width=6.75,
       units="in")

setwd(fun_dir)
