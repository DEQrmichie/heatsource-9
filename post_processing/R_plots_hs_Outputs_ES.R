#---------------------------------------------------------------------------------------
# Import heat source effective shade outputs and compare to observed data on a plot
# output a csv file for a certain date that can be used for import into GIS
#---------------------------------------------------------------------------------------

#library(reshape2)
#library(reshape)
library(ggplot2)

name <- "Yach_s1_2"

plot.date <- c("08/15/2005")

out_dir <-"F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_2/outputs/"
sim_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_2/outputs/"
obs_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/"
fun_dir <- "E:/GitHub/Rscripts/heatsource"

source(paste0(fun_dir,"/","R_functions_hs_read_and_format.R"))

setwd(out_dir)

sim <- read.hs.outputs(sim_dir,"shade", "Effective Shade",name, hs8=FALSE)

obs.raw <- read.obs(obs_dir=obs_dir, file_name="ObservedData.csv")

obs <- obs.raw[obs.raw$ConstituentCode == "ES",]

# data frame of obkm and simkm with site name
km <- obs2simkm(obs_dir=obs_dir, file_name="ObservedData.csv", 
                ConstituentCode="ES",
                simkm=as.numeric(unique(sim$Stream_km)))


####################################################################################
# Effective Shade  Data


# convert to percent
obs$value <- obs$value*100
sim$value <- sim$value*100

# Drop these cols
obs$Notes <- NULL
obs$ConstituentCode <- NULL
obs$SiteName <- NULL
obs$Constituent <- NULL
obs$Datetime <- NULL
obs$Date <-NULL
sim$Datetime <- NULL
sim$hour <- NULL

# Build a list of dates when the effective shade data was collected.
#obs.dates <- unique(obs$Date)


#km <- cbind(data.frame(obs.km),data.frame(sim.km))
obs <- merge(km,obs,by.x="obs.km",by.y="Stream_km")

obs$Stream_km <- obs$sim.km
obs$obs.km <- NULL
obs$sim.km <- NULL

# Just get predictions on the plot date
sim.plot <- sim[sim$Date %in% plot.date,]

df <- merge(sim.plot, obs, by=c("Stream_km"),all=TRUE)
colnames(df) <- c("Stream_km", "sim", "constituent", "es_sim", "Date", "SiteName", "es_obs") 

#df <- merge(sim.l,obs, by=c("Date","Stream_km"),all=TRUE)
#colnames(df) <- c("Date", "Stream_km","es_sim", "es_obs") 

df$Date <- factor(df$Date)

es.lm <- lm(es_obs~es_sim,data=df)

# write the fit stats to a text file
sink(paste0(out_dir,name,"_ES2_STATS.txt"))
summary(es.lm)
sink()

p.es <- ggplot(data=df, aes(x=Stream_km)) +
  geom_line(aes(y=es_sim, color="Predictions")) +
  geom_point(aes(y=es_obs, color="Observations"),size= 2.0) +
  guides(color=guide_legend(override.aes=list(shape=c(16,NA),linetype=c(0,1))))  +
  scale_colour_manual(values=c("Predictions"="red","Observations"="Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Effective Shade (%)") + 
  ylim(0,100) +
  facet_wrap(~Date, nrow = length(plot.date))

p2.es <- ggplot(data=df, aes(x=es_sim,y=es_obs)) +
  geom_point(shape=1) +
  geom_smooth(method="lm") +
  ylab("Observed Effective Shade (%)") + 
  xlab("Predicted Effective Shade (%)")

p2.es

ggsave(file=paste0(out_dir,name,"_ES.png"),
       plot=p.es,
       height=3,
       width=6.75,
       units="in")

ggsave(file=paste0(out_dir,name,"_ES2.png"),
       plot=p2.es,
       height=3,
       width=6.75,
       units="in")

# Output csv for GIS
write.csv(sim.plot, file=paste(out_dir,name,"_ES_",gsub("/","_", plot.date),".csv",sep=""),row.names = FALSE)

setwd(fun_dir)
