#---------------------------------------------------------------------------------------
# Import heat source 8 model flows and compare to observed data on a longitudinal plot
#---------------------------------------------------------------------------------------

library(reshape2)
library(reshape)
library(ggplot2)

name <- "Yach_s1_16"

outpath <-"F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_16/outputs/"
sim_path <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_16/outputs/"
obs_path <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/"

setwd(outpath)

sim.raw <- read.table(paste0(sim_path,"Hyd_Flow.txt"),
                      sep="",dec=".",skip=2,header=TRUE, 
                      stringsAsFactors = FALSE, na.strings = "NA")

obs.raw <- read.table(paste0(obs_path,"ObservedData.csv"),
                      sep=",",dec=".",skip=0,header=TRUE, 
                      stringsAsFactors = FALSE, na.strings = "NA")


# Format Datetime
sim.raw$Datetime <- round(as.POSIXct('1899-12-30')+(floor(sim.raw$Datetime)
                                                    *60*60*24)-3600,"mins")
obs.raw$Datetime <- as.POSIXct(obs.raw$Datetime,format="%m/%d/%Y %H:%M") #6/17/2003 0:00


# Convert data from wide to long
sim.l <- melt(sim.raw, id.vars =c("Datetime"),variable_name="Stream_km")
sim.l$Stream_km <- as.numeric(gsub(pattern="X", replacement="", sim.l$Stream_km, ignore.case = FALSE,fixed = FALSE))

# Get the observed flow rate
obs <- obs.raw[obs.raw$ParameterCode == "Q",]


####################################################################################
# Date formatting
obs$Date <- format(obs$Datetime, "%m/%d/%Y")
sim.l$Date <- format(sim.l$Datetime, "%m/%d/%Y")


# Calculate the daily average flow rate
sim <- as.data.frame(aggregate(sim.l$value,by=list(sim.l$Date,sim.l$Stream_km), FUN=mean,na.rm = TRUE))
colnames(sim)[1] <- "Date"
colnames(sim)[2] <- "Stream_km"
colnames(sim)[3] <- "value"

# Only use observed data for the same time period as the model
Tstart <-min(sim.raw$Datetime)
Tend <-max(sim.raw$Datetime)

obs <- obs[obs$Datetime >= Tstart & obs$Datetime <= Tend,]

# convert to cfs
obs$value <- obs$value * 35.3147
sim$value <- sim$value * 35.3147

# Drop these cols
obs$Notes <- NULL
obs$ParameterCode <- NULL
obs$SiteName <- NULL
obs$ParameterName <- NULL
obs$Datetime <- NULL
sim.l$Datetime <- NULL

# Add these cols
#sim.l$ParameterName <- unique(obs$ParameterName)

#obs.wide  <- reshape(obs, idvar = c("ParameterName","Date"), 
#                     direction="wide", drop = c("Datetime","value","group"))


# Build a list of dates when the flow rate data was collected.
obs.dates <- unique(obs$Date)

obs.km <- as.numeric(unique(obs$Stream_km))
simkm <- as.numeric(unique(sim.l$Stream_km))

sim.km <- numeric(0)

for (i in 1:length(obs.km)) {
  # this makes a vector of the model km that is closest to the obs km
  sim.km[i] <-simkm[which(abs(simkm-obs.km[i])==min(abs(simkm-obs.km[i])))[1]]
}

km <- cbind(data.frame(obs.km),data.frame(sim.km))
obs <- merge(km,obs,by.x="obs.km",by.y="Stream_km")

obs$Stream_km <- obs$sim.km
obs$obs.km <- NULL
obs$sim.km <- NULL

# Just get predictions on the observation dates
sim <- sim[sim$Date %in% obs.dates,]

df <- merge(sim,obs, by=c("Date","Stream_km"),all=TRUE)
colnames(df) <- c("Date", "Stream_km","sim", "obs") 

#df <- rbind(obs, sim)
df$Date <- factor(df$Date)

p.flow <- ggplot(data=df, aes(x=Stream_km)) +
  geom_line(aes(y=sim, color="Model")) +
  geom_point(aes(y=obs, color="Observations"),size= 2.0) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,16),linetype=c(1,0))))  +
  scale_colour_manual(values=c("Observations"="Black","Model"="red")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Flow Rate (cfs)") + 
  facet_wrap(~Date, nrow = length(obs.dates))

ggsave(file=paste0(outpath,name,"_FLOW.png"),
       plot=p.flow,
       height=6.75,
       width=6.75,
       units="in")

# for(i in 1:length(obs.dates)) {
#   
#   obs.i <- obs[obs$Date == obs.dates[i],]
#   sim.i <- sim.l[sim.l$Date == obs.dates[i],]
#   data.i <- rbind(obs.i,sim.i)
#   data.i$group <- factor(data.i$group)
#   
#   p1 <- ggplot(data=data.i, aes(x=Stream_km, y=value, colour=group, linetype=group)) +
#     geom_line() +
#     geom_point(aes(shape=group),size=2) +
#     scale_shape_manual(values=c(16,NA)) + 
#     scale_color_manual(values=c("black", "red")) +
#     scale_linetype_manual(values=c(0,1)) +
#     scale_y_continuous(limits = c(0, 100)) +
#     xlab("Model Stream Kilometer") +
#     ylab("Effective Shade (%)") + 
#     theme(text=element_text(size=10)) +
#     theme(legend.position="bottom") +
#     theme(legend.title=element_blank())
#   
#   ggsave(file=paste0(outpath,name,"_ES",i,".png"),
#          plot=p1,
#          height=3,
#          width=6.75,
#          units="in")
# }
