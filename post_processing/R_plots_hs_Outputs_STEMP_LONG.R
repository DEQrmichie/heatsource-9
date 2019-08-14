#---------------------------------------------------------------------------------------
# Plot heat source output hourly and 7DADM stream temperatures and observed temps longitudinally by hour
# for selected dates. Plot output model period min, max, and mean 7DADM stream temps longitudnally by stream km,
#---------------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)

sim_name <- "Yach_s1_2"

sim_file <- "Temp_H2O"
obs_file <- "ObservedData.csv"
bbnc <- 18

plot.dates <- c("07/01/2005","07/15/2005","08/01/2005","08/15/2005","09/01/2005")

out_dir <-"F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_2/outputs/"
sim_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_2/outputs/"
obs_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/"
fun_dir <- "E:/GitHub/Rscripts/heatsource"

source(paste0(fun_dir,"/","R_functions_hs_read_and_format.R"))

setwd(out_dir)

#--  Read temps ------------------------------------
sim.hour <- read.hs.outputs(output_dir=sim_dir, file_name="Temp_H2O", 
                            constituent_name="Temperature (deg-C)",
                            statistic_name="hourly",
                            sim_name=sim_name, hs_ver=9, sheet_name = NULL)

sim.7dadm <- calc.7dadm(sim.hour)

# calc summary 7DADM
sim.7dadm.summ <- calc.summary(sim.7dadm)

# Read observation data
obs <- read.obs(obs_dir=obs_dir, file_name=obs_file)
obs <- obs[obs$constituentCode == "STEMPH",]
obs$hour <-as.integer(format(obs$Datetime, "%H"))

# calc obs 7DADM
obs$sim <- "Observations"
obs.7dadm <- calc.7dadm(df=obs)

# Just get data on the plot date
sim.hour <- sim.hour[sim.hour$Date %in% plot.dates,]
obs <- obs[obs$Date %in% plot.dates,]

sim.7dadm<- sim.7dadm[sim.7dadm$Date %in% plot.dates,]
obs.7dadm <- obs.7dadm[obs.7dadm$Date %in% plot.dates,]



# data frame of obkm and simkm with site name
km <- obs2simkm(obs_dir=obs_dir, file_name=obs_file, 
                constituentCode="STEMPH",
                simkm=as.numeric(unique(sim.hour$Stream_km)))

# Drop these cols
obs$notes <- NULL
obs$constituentCode <- NULL
obs$constituent <- NULL
obs$statistic <- NULL
obs$siteName <-NULL
obs$Datetime <- NULL

obs <- merge(km, obs, by.x=c("obs.km"), by.y=c("Stream_km"))

obs$Stream_km <- obs$sim.km
obs$obs.km <- NULL
obs$sim.km <- NULL

df <- merge(sim.hour, obs, by=c("Date","hour","Stream_km"), all.x=TRUE)

df$sim.y <- NULL
df$siteName <- NULL

colnames(df) <- c("Date","hour","Stream_km","Datetime","sim","constituent","statistic", "pred","obs")
df$Date <- format(df$Datetime, "%m/%d/%Y")
df$hour <- format(df$Datetime, "%H:%M")

obs.7dadm <- merge(obs.7dadm, km, by.x="Stream_km",by.y="obs.km")

obs.7dadm$Stream_km <- obs.7dadm$sim.km
obs.7dadm$sim.km <- NULL
obs.7dadm$siteName <- NULL
obs.7dadm$sim <- NULL

df.7dadm <- merge(sim.7dadm, obs.7dadm, by=c("Date","Stream_km","constituent", "statistic"),all.x=TRUE)
colnames(df.7dadm) <- c("Date","Stream_km","constituent","statistic","sim","pred","obs")

# set y axis plot limits
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

ymin <- round_any(min(df$pred, df$obs, na.rm=TRUE), 5, floor)
ymax <- round_any(max(df$pred, df$obs, na.rm=TRUE), 5, ceiling)

#---------------------------------------------------------------------------------------
# PLOT HOURLY SIMULATED VS. OBSERVED TEMPERATURES #

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
    ylim(ymin, ymax) +
    facet_wrap(~hour,nrow = 6)
  
  ggsave(file=paste0(out_dir,sim_name,"_LONG_HOURLY_",gsub("/","_", plot.dates[i]),".png"),
         plot=p.long,
         height=7.2,
         width=7,
         units="in")
}

#---------------------------------------------------------------------------------------
# PLOT 7DADM SIMULATED VS. OBSERVED TEMPERATURES

df.i <- df.7dadm[df.7dadm$statistic== "7DADM Temperature",]

p.long2 <- ggplot(data=df.i, aes(x=Stream_km)) +
  geom_line(aes(y=pred, color="Predictions")) +
  geom_point(aes(y=obs, color="Observations"),size=.75) +
  guides(color=guide_legend(override.aes=list(shape=c(16,NA),linetype=c(0,1))))  +
  scale_colour_manual(values=c("Predictions"="red","Observations"="Black")) +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.title.y =element_blank(),
        plot.title = element_text(size=12, hjust = 0.5)) +
  xlab("Model Stream Kilometer") +
  ylab("Temperature (C)") +
  #ylim(5,20) +
  ylim(ymin, ymax) +
  facet_wrap(~Date, nrow = length(plot.dates))

ggsave(file=paste0(out_dir,sim_name,"_LONG_7DADM.png"),
       plot=p.long2,
       height=7.2,
       width=7,
       units="in")

#---------------------------------------------------------------------------------------
# PLOT MIN, MEAN, and MAX 7DADM SIMULATED TEMPERATURES FOR THE ENTIRE MODEL PERIOD

df.summ <- sim.7dadm.summ[sim.7dadm.summ$statistic== "7DADM Temperature",]

p.summ <- ggplot(data=df.summ, aes(x=Stream_km)) +
  geom_hline(aes(yintercept=bbnc, linetype="Applicable Criteria")) +
  geom_ribbon(aes(ymax = max, ymin = min, fill="Min and Max Range"), alpha = 0.6) + 
  geom_line(aes(y=median, linetype="Median")) +
  scale_linetype_manual(values=c("Applicable Criteria"="dashed","Median"="solid")) +
  scale_fill_manual(values="darkgrey") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.title.y =element_blank(),
        panel.background = element_rect(fill="white", colour = "black"),
        strip.background =element_rect(fill="white", colour = "black"),
        panel.grid.major =element_line(colour = "lightgrey"),
        plot.title = element_text(size=12, hjust = 0.5)) +
  xlab("Model Stream Kilometer") +
  ylab("7DADM Temperature (C)") + 
  ylim(ymin, ymax)

ggsave(file=paste0(out_dir,sim_name,"_LONG_7DADM_SUMM.png"),
       plot=p.summ,
       height=3,
       width=6.75,
       units="in")

setwd(fun_dir)
