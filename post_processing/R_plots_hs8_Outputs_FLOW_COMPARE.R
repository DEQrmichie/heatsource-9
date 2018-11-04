
#---------------------------------------------------------------------------------------
# Import heat source 8 stream flow rate output from two simulations and plot results
#---------------------------------------------------------------------------------------

library(reshape2)
library(ggplot2)
library(zoo)

name <- "Yach_s2_03_01"

plot.date <- c("08/15/2005")

sim1.lab <- "Current Flow"
sim2.lab <- "Restored Flow"

outpath <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_03_FLOW/s3_01/outputs/"
sim1_path <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_16/outputs/"
sim2_path <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_03_FLOW/s3_01/outputs/"

setwd(outpath)

sim1.raw <- read.table(paste0(sim1_path,"Hyd_Flow.txt"),
                       sep="",dec=".",skip=2,header=TRUE, 
                       stringsAsFactors = FALSE, na.strings = "NA")

sim2.raw <- read.table(paste0(sim2_path,"Hyd_Flow.txt"),
                       sep="",dec=".",skip=2,header=TRUE, 
                       stringsAsFactors = FALSE, na.strings = "NA")

# Convert data from wide to long
sim1.l <- melt(sim1.raw, id.vars =c("Datetime"),variable_name="Stream_km")
colnames(sim1.l) <- c("Datetime", "Stream_km","flow1")
sim1.l$Stream_km <- as.numeric(gsub(pattern="X", replacement="", sim1.l$Stream_km, ignore.case = FALSE,fixed = FALSE))

sim2.l <- melt(sim2.raw, id.vars =c("Datetime"),variable_name="Stream_km")
colnames(sim2.l) <- c("Datetime", "Stream_km","flow2")
sim2.l$Stream_km <- as.numeric(gsub(pattern="X", replacement="", sim2.l$Stream_km, ignore.case = FALSE,fixed = FALSE))

# Format Datetime
sim1.l$Datetime <- round(as.POSIXct('1899-12-30')+(floor(sim1.l$Datetime)
                                                    *60*60*24)-3600,"day")
sim2.l$Datetime <- round(as.POSIXct('1899-12-30')+(floor(sim2.l$Datetime)
                                                    *60*60*24)-3600,"day")

# Date formatting
sim1.l$Date <- format(sim1.l$Datetime, "%m/%d/%Y")
# Date formatting
sim2.l$Date <- format(sim2.l$Datetime, "%m/%d/%Y")


# Calculate the daily average flow rate
sim1 <- as.data.frame(aggregate(sim1.l$flow1,by=list(sim1.l$Date,sim1.l$Stream_km), 
                                FUN=mean,na.rm = TRUE))
colnames(sim1)[1] <- "Date"
colnames(sim1)[2] <- "Stream_km"
colnames(sim1)[3] <- "flow1"


sim2 <- as.data.frame(aggregate(sim2.l$flow2,by=list(sim2.l$Date,sim2.l$Stream_km), 
                                FUN=mean,na.rm = TRUE))
colnames(sim2)[1] <- "Date"
colnames(sim2)[2] <- "Stream_km"
colnames(sim2)[3] <- "flow2"

# sort by Stream_km and Date
sim1 <- sim1[with(sim1, order(-Stream_km,Date)), ]
sim2 <- sim2[with(sim2, order(-Stream_km,Date)), ]

df <- merge(sim1, sim2, by=c("Date","Stream_km"),all=TRUE)

colnames(df) <- c("Date","Stream_km", "flow1", "flow2")

# convert from cms to cfs
df$flow1 <- df$flow1*35.3147
df$flow2 <- df$flow2*35.3147

# Change in flow between two simulations
df$dflow <- df$flow2 - df$flow1

# Just get predictions on the plot date
df.plot <- df[df$Date %in% plot.date,]

# longitudinal plot of the change in flow
p.dflow <- ggplot(data=df.plot, aes(x=Stream_km)) +
  geom_line(aes(y=dflow)) +
  scale_linetype_manual(values=c("solid")) +
  scale_fill_manual(values="skyblue") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Change in flow (cfs)") +
  facet_wrap(~Date, nrow = length(plot.date))
p.dflow 

ggsave(file=paste0(outpath,name,"_FLOW_CHANGE.png"),
       plot=p.dflow,
       height=3,
       width=6.75,
       units="in")

# longitudinal plot of flow by date
p.flow <- ggplot(data=df.plot, aes(x=Stream_km)) +
  geom_ribbon(aes(ymax = flow2, ymin = flow1, fill="Change"), alpha = 0.6) + 
  geom_line(aes(y=flow1, linetype="flow1")) +
  geom_line(aes(y=flow2, linetype="flow2")) +
  guides(color=guide_legend(override.aes=list(linetype=c(0,1))))  +
  scale_linetype_manual(values=c("flow1"="solid", "flow2"="dotted"),
                        labels=c(sim1.lab, sim2.lab)) +
  scale_fill_manual(values="skyblue") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Flow Rate (cfs)") + 
  facet_wrap(~Date, nrow = length(plot.date))
p.flow

ggsave(file=paste0(outpath,name,"_FLOW_COMPARE.png"),
       plot=p.flow,
       height=3,
       width=6.75,
       units="in")

