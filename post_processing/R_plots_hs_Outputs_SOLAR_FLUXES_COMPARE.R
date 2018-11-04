#---------------------------------------------------------------------------------------
# Import heat source 8 hourly flux outputs from two simulations
# and plot the changes
#---------------------------------------------------------------------------------------

library(reshape2)
library(ggplot2)
library(plyr)
#library(dplyr)

fun_dir <- "E:/GitHub/Rscripts/heatsource"

name <- "Example_River"

sim1_name <- "hs8.0.8"
sim2_name <- "hs9.0.0"

plot.dates <- c("07/01/2001")
#plot.dates <- c("08/01/2005")

plot.sims <- c(sim1_name,sim2_name)

out_dir <-"C:/workspace/Heatsource_example/comparisions/"
sim1_dir <- "C:/workspace/Heatsource_example/outputs_808/"
sim2_dir <- "C:/workspace/Heatsource_example/hs9/outputs/"
obs_dir <- "C:/workspace/Heatsource_example/"

#out_dir <-"C:/workspace/Heatsource_example/comparisions/"
#sim1_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_16/outputs/"
#sim2_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_6/outputs/"
#obs_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/"

source(paste0(fun_dir,"/","R_functions_hs_read_and_format.R"))

# use km from obs sites for plot km
km <- obs2simkm(obs_dir=obs_dir,file_name="ObservedData.csv", 
                ConstituentCode= "STEMPH",
                simkm=as.numeric(unique(flux$Stream_km)))

plot.kms <- km$sim.km

setwd(out_dir)

#--  Read sim fluxes ------------------------------------

flux1 <- read.solar.flux.outputs(sim1_dir,sim1_name,hs8=TRUE)
flux2 <- read.solar.flux.outputs(sim2_dir,sim2_name,hs8=FALSE)

# Add character Date
flux1$Date <- format(flux1$Datetime, "%m/%d/%Y")
flux2$Date <- format(flux2$Datetime, "%m/%d/%Y")

# Add hour
flux1$hour <- flux1$Datetime$hour
flux2$hour <- flux2$Datetime$hour

fluxd <- merge(flux1,flux2, by=c("Date", "hour", "Stream_km", "constituent"))

#find the difference
fluxd$value <- fluxd$value.y - fluxd$value.x

fluxd <- fluxd[,c("Datetime.x","sim.x","Stream_km","constituent","value","Date", "hour")]

fluxd <- rename(fluxd, c("Datetime.x"="Datetime","sim.x"="sim"))
fluxd$sim <- "Change"

flux <- rbind(flux1,flux2,fluxd)

rm(flux1,flux2)

#--  subset data for plotting ------------------------------------

df <- flux[flux$Stream_km %in% plot.kms & flux$Date %in% plot.dates & flux$sim %in% c(plot.sims,"Change"),]
df <- merge(df,km, by.x="Stream_km",by.y="sim.km")

levels.hs8 <- c("SR1: Solar Radiation above Topo", 
                "SR4: Solar Radiation above Stream",
                "SR6: Solar Radiation received by Stream")

levels.hs9 <- c("SR1: Solar Radiation above Topo",
                "SR2: Solar Radiation below Topo",
                "SR3: Solar Radiation below land cover",
                "SR4: Solar Radiation above Stream",
                "SR5: Solar Radiation entering Stream",
                "SR6: Solar Radiation received by Stream")

df$constituent <- factor(df$constituent,levels=levels.hs8)
df$sim <- factor(df$sim,levels=c(plot.sims,"Change"))

#--  plots ------------------------------------
i <- 1

for (i in 1:length(plot.kms)) {
  
  df.i <- df[df$Stream_km==plot.kms[i] & df$sim %in% plot.sims,]
  
  p.flux <- ggplot(data=df.i,aes(x=Datetime,y=value,color=constituent, size=constituent)) +
    geom_line() +
    scale_color_manual(values=c("black","black","black")) +
    scale_size_manual(values = c(0.5,0.5,0.5)) +
    #scale_color_manual(values=c("red","yellow","light blue","dark blue", "green","black")) +
    #scale_size_manual(values = c(0.5,0.5,0.5,0.5,0.5,0.5)) +
    #labs(title=unique(df.i$SiteName),
    labs(title=paste0("Model km ",unique(df.i$Stream_km)),
         subtitle=unique(df.i$Date)) +
    theme(legend.position="bottom") +
    theme(legend.title=element_blank()) +
    theme(plot.title=element_text(size=10,face="bold", color ="black")) +
    xlab("Hour") +
    scale_x_datetime(date_breaks = "4 hours", date_labels = "%H") +
    ylab("Flux (watts/square meter)") +
    ylim(round_any(min(df$value), 100, floor),round_any(max(df$value), 100, ceiling)) +
    facet_wrap(~sim+constituent, nrow=2, scales = "fixed")
  
  p.flux
  
  
  df.i <- df[df$Stream_km==plot.kms[i] & df$sim %in% c("Change"),]
  
  p.flux.change <- ggplot(data=df.i, aes(x=Datetime,y=value,color=constituent, size=constituent)) +
    geom_line() +
    scale_color_manual(values=c("black","black","black")) +
    scale_size_manual(values = c(0.5,0.5,0.5)) +
    #scale_color_manual(values=c("red","yellow","light blue","dark blue", "green","black")) +
    #scale_size_manual(values = c(0.5,0.5,0.5,0.5,0.5,0.5)) +
    #labs(title=unique(df.i$SiteName),
    labs(title=paste0("Model km ",unique(df.i$Stream_km)),
         subtitle=paste0(unique(df.i$Date),"   ",sim2_name," minus ",sim1_name)) +
    theme(legend.position="bottom") +
    theme(legend.title=element_blank()) +
    theme(plot.title=element_text(size=10,face="bold", color ="black")) +
    xlab("Hour") +
    scale_x_datetime(date_breaks = "4 hours", date_labels = "%H") +
    ylab("Flux (watts/square meter)") +
    ylim(round_any(min(df.i$value), 10, floor),round_any(max(df.i$value), 10, ceiling)) +
    facet_wrap(~sim+constituent, nrow=1)
  
  p.flux.change
  
  
  ggsave(file=paste0(out_dir,name,"_SOLAR_FLUX_compare_km_",plot.kms[i],".png"),
         plot=p.flux,
         height=5,
         width=8,
         units="in")
  
  ggsave(file=paste0(out_dir,name,"_SOLAR_FLUX_change_km_",plot.kms[i],".png"),
         plot=p.flux.change,
         height=5,
         width=8,
         units="in")
  
}

