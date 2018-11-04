
#---------------------------------------------------------------------------------------
# Import heat source effective shade output from two simulations
# and plot the change and absolute results.
#---------------------------------------------------------------------------------------

library(reshape2)
library(ggplot2)
library(zoo)

name <- "Yach_s2_02_03"

plot.date <- c("08/15/2005")

sim1_name <- "Current Condition"
sim2_name <- "Restored Vegetation"

fun_dir <- "E:/GitHub/Rscripts/heatsource"
out_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_02_VEG/s2_03/outputs_solar/"
sim1_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_1/hs8/outputs_solar/"
sim2_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_02_VEG/s2_03/outputs_solar/"

source(paste0(fun_dir,"/","R_functions_hs_read_and_format.R"))

#setwd(out_dir)

#--  Read data ------------------------------------

sim1 <- read.hs.outputs(output_dir=sim1_dir, file_name="Shade", 
                            constituent_name="Effective Shade", 
                            sim_name=sim1_name, hs8=TRUE)

sim2 <- read.hs.outputs(output_dir=sim2_dir, file_name="Shade", 
                            constituent_name="Effective Shade", 
                            sim_name=sim2_name, hs8=TRUE)

sim1$value <- sim1$value*100
sim2$value <- sim2$value*100

# sort by Stream_km and Date
sim1 <- sim1[with(sim1, order(-Stream_km,Datetime)), ]
sim2 <- sim2[with(sim2, order(-Stream_km,Datetime)), ]

df <- merge(sim1, sim2, by=c("Datetime", "Date", "hour", "Stream_km", "constituent"),all=TRUE)

colnames(df) <- c("Datetime", "Date", "hour", "Stream_km", "constituent", "sim1", "es1", "sim2", "es2")

# Change in Effective Shade
df$change<- df$es2 - df$es1

# Just get predictions on the plot date
df.plot <- df[df$Date %in% plot.date,]

# longitudinal plot of change
p.change <- ggplot(data=df.plot, aes(x=Stream_km)) +
  geom_line(aes(y=change)) +
  scale_linetype_manual(values=c("solid")) +
  scale_fill_manual(values="skyblue") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Change in Effective Shade % points") +
  facet_wrap(~Date, nrow = length(plot.date))
p.change 

ggsave(file=paste0(out_dir,name,"_ES_CHANGE.png"),
       plot=p.change,
       height=3,
       width=6.75,
       units="in")

# longitudinal plot of ES
p.es <- ggplot(data=df.plot, aes(x=Stream_km)) +
  geom_ribbon(aes(ymax = es2, ymin = es1, fill="Change"), alpha = 0.6) + 
  geom_line(aes(y=es1, linetype="es1")) +
  geom_line(aes(y=es2, linetype="es2")) +
  guides(color=guide_legend(override.aes=list(linetype=c(0,1))))  +
  scale_linetype_manual(values=c("es1"="solid", "es2"="dotted"),
                        labels=c(sim1_name, sim2_name)) +
  scale_fill_manual(values="skyblue") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Effective Shade %") + 
  ylim(0,100) +
  facet_wrap(~Date, nrow = length(plot.date))
p.es

ggsave(file=paste0(out_dir,name,"_ES_COMPARE.png"),
       plot=p.es,
       height=3,
       width=6.75,
       units="in")

# Output csv for GIS
write.csv(df.plot[,c(4,2,5:10)], file=paste(out_dir,name,"_ES_COMPARE_",gsub("/","_", plot.date),".csv",sep=""),row.names = FALSE)

