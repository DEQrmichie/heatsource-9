
#---------------------------------------------------------------------------------------
# Import heat source hourly temperature outputs from two simulations, 
# calculate 7DADM, calculae the change, and plot longitudinal data
#---------------------------------------------------------------------------------------

library(ggplot2)
library(plyr)

fun_dir <- "E:/GitHub/Rscripts/heatsource"

#name <- "Example_River"
name <- "Yach"

sim1_name <- "hs8.0.8"
sim2_name <- "hs9.0.0"

plot.sims <- c(sim1_name,sim2_name)

#out_dir <-"C:/workspace/Heatsource_example/comparisions/"
#sim1_dir <- "C:/workspace/Heatsource_example/outputs_808/"
#sim2_dir <- "C:/workspace/Heatsource_example/hs9/outputs/"

out_dir <-"C:/workspace/Heatsource_example/comparisions/"
sim1_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_17_test/outputs/"
sim2_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_6/outputs/"

source(paste0(fun_dir,"/","R_functions_hs_read_and_format.R"))

setwd(out_dir)

#--  Read temps and calc 7dadm ------------------------------------

data1 <- calc.7dadm2(sim1_dir,sim1_name,hs8 = TRUE)
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

rm(data1,data2,datad)

#--  Calc Summary Stats ------------------------------------

# calc max, mean, min for all dates by sim, constituent, and stream_km
data0.summary <- calc.summary(data=data0)

#-- Plot longitudinal summary for each simulation and differnece between them  ---------------


# longitudinal plot of dT summary
p.dT<- ggplot(data=data0.summary[data0.summary$sim == "Change" & data0.summary$constituent== "7DADM Temperature",], aes(x=Stream_km)) +
  geom_ribbon(aes(ymax = max, ymin = min, fill="Min & Max Range"), alpha = 0.6) + 
  geom_line(aes(y=mean, linetype="Mean")) +
  geom_line(aes(y=0.3, linetype="HUA")) +
  scale_linetype_manual(values=c("Mean"="solid", "HUA"="dashed")) +
  scale_fill_manual(values="skyblue") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Change 7DADM (C)") +
  facet_wrap(~sim+constituent, nrow=2)
p.dT 

ggsave(file=paste0(out_dir,name,"_dT_7DADM.png"),
       plot=p.dT,
       height=3,
       width=6.75,
       units="in")

setwd(fun_dir)
