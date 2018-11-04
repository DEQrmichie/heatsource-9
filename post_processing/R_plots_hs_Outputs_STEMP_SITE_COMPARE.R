#---------------------------------------------------------------------------------------
# Import heat source hourly temperature outputs from two simulations, 
# calculate 7DADM, calculate the change at specfic sites, and plot
#---------------------------------------------------------------------------------------
library(ggplot2)
library(plyr)

fun_dir <- "E:/GitHub/Rscripts/heatsource"

sim1_name <- "s1_1"
sim2_name <- "s1_2"

plot.sims <- c(sim1_name,sim2_name)
plot.constituents <-c("7DADM Temperature")

name <- "Yach"
out_dir <-"F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_2/outputs/"
sim1_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_1/outputs/"
sim2_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_2/outputs/"
obs_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/"

#name <- "Example_River"
#out_dir <-"C:/workspace/Heatsource_example/comparisions/"
#sim1_dir <- "C:/workspace/Heatsource_example/outputs_808/"
#sim2_dir <- "C:/workspace/Heatsource_example/hs9/outputs/"
#obs_dir <- "C:/workspace/Heatsource_example/"

source(paste0(fun_dir,"/","R_functions_hs_read_and_format.R"))

setwd(out_dir)

#--  Read temps and calc 7dadm ------------------------------------

data1 <- calc.7dadm2(sim1_dir,sim1_name,hs8 = FALSE)
data2 <- calc.7dadm2(sim2_dir,sim2_name,hs8 = FALSE)

datad <- merge(data1,data2, by=c("Date", "Stream_km", "constituent"))

# find the difference
datad$value <- datad$value.y - datad$value.x

datad <- datad[,c("Date","sim.x","Stream_km","constituent","value")]

datad <- rename(datad, c("sim.x"="sim"))
datad$sim <- "Change"

data0 <- rbind(data1,data2,datad)

# convert the character date back to POSIXct
data0$Datetime <- as.POSIXct(data0$Date,format="%m/%d/%Y") #6/17/2003

#--  Read obs sites ------------------------------------

obs.raw <- read.obs(obs_dir=obs_dir, file_name="ObservedData.csv")

obs <- obs.raw[obs.raw$ConstituentCode == "STEMPH",]

# Drop these cols
obs$Notes <- NULL
obs$ConstituentCode <- NULL
obs$Constituent <- NULL

# just get the sites and stream_km
obs.wide  <- reshape(obs, timevar= "Datetime", idvar = c("Stream_km","SiteName"), 
                     direction="wide", drop = c("Date","value","group"))

#-- Plot 7dadm Temperature at obs sites ---------------

obs.km <- as.numeric(unique(obs$Stream_km))
simkm <- as.numeric(unique(data0$Stream_km))

sim.km <- numeric(0)

for (i in 1:length(obs.km)) {
  # this makes a vector of the model km that is closest to the obs km
  sim.km[i] <-simkm[which(abs(simkm-obs.km[i])==min(abs(simkm-obs.km[i])))[1]]
}

km <- cbind(data.frame(obs.km),data.frame(sim.km))
km <- merge(km,obs.wide,by.x="obs.km",by.y="Stream_km")

df <- data0[data0$Stream_km %in% sim.km & 
             data0$sim %in% c(plot.sims,"Change") &
             data0$constituent %in% plot.constituents,]
df <- merge(df,km, by.x="Stream_km",by.y="sim.km")
#df$constituent <- factor(df$constituent,levels=c("Daily Maximums", "7DADM"))
df$sim <- factor(df$sim,levels=c(plot.sims,"Change"))

# set y axis plot limits
ymin <- round_any(min(df[df$sim %in% plot.sims,]$value, na.rm=TRUE), 5, floor)
ymax <- round_any(max(df[df$sim %in% plot.sims,]$value, na.rm=TRUE), 5, ceiling)

i <- 1

for (i in 1:length(sim.km)) {
  
  df.i <- df[df$Stream_km==sim.km[i] & df$sim %in% plot.sims,]
  
  p.data <- ggplot(data=df.i,aes(x=Datetime,y=value,color=constituent, size=constituent)) +
    geom_point() +
    scale_color_manual(values=c("black","black")) +
    scale_size_manual(values = c(0.5,0.5)) +
    labs(title=unique(df.i$SiteName)) +
    theme(legend.position="bottom") +
    theme(legend.title=element_blank()) +
    theme(plot.title=element_text(size=10,face="bold", color ="black")) +
    xlab("Day") +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y") +
    ylab("Temperature (C)") +
    ylim(ymin,ymax) +
    facet_wrap(~sim+constituent, nrow=2, scales = "fixed")
  
  p.data
  
  df.i <- df[df$Stream_km==sim.km[i] & df$sim %in% c("Change"),]
  
  p.data.change <- ggplot(data=df.i, aes(x=Datetime,y=value,color=constituent, size=constituent)) +
    geom_point() +
    scale_color_manual(values=c("black","black")) +
    scale_size_manual(values = c(0.5,0.5)) +
    labs(title=unique(df.i$SiteName),
         subtitle=paste0(sim2_name," minus ",sim1_name)) +
    theme(legend.position="bottom") +
    theme(legend.title=element_blank()) +
    theme(plot.title=element_text(size=10,face="bold", color ="black")) +
    xlab("Day") +
    scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y") +
    ylab("Temperature (C)") +
    ylim(round_any(min(df.i$value, na.rm=TRUE), 5, floor),round_any(max(df.i$value, na.rm=TRUE), 5, ceiling)) +
    facet_wrap(~sim+constituent, nrow=1)
  
  p.data.change
  
  
  ggsave(file=paste0(out_dir,name,"_7DADM_compare_km_",sim.km[i],".png"),
         plot=p.data,
         height=5,
         width=8,
         units="in")
  
  ggsave(file=paste0(out_dir,name,"_7DADM_change_km_",sim.km[i],".png"),
         plot=p.data.change,
         height=5,
         width=8,
         units="in")
  
}


