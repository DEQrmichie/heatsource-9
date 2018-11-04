#---------------------------------------------------------------------------------------
# Plot heat source 9 accretion and trib flow and temp:
#---------------------------------------------------------------------------------------

library (reshape2)
library (ggplot2)

sim_name <- "Yach_s1_16"

out_dir <-"F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_16/inputs_as_csv/"
sim_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_16/inputs_as_csv/"
obs_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/"

flows.p <- read.table(paste0(sim_dir,"inflows_for_plotting.csv"),
                     sep=",",dec=".",skip=0,header=TRUE, 
                     stringsAsFactors = FALSE, na.strings = "NA")

inflow <- read.table(paste0(sim_dir,"inflow.csv"),
                    sep=",",dec=".",skip=0,header=TRUE, 
                    stringsAsFactors = FALSE, na.strings = "NA")

inflow.n <- read.table(paste0(obs_dir,"inflow_natural.csv"),
                     sep=",",dec=".",skip=0,header=TRUE, 
                     stringsAsFactors = FALSE, na.strings = "NA")

acc <- read.table(paste0(sim_dir,"accretion.csv"),
                     sep=",",dec=".",skip=0,header=TRUE, 
                     stringsAsFactors = FALSE, na.strings = "NA")

obs.raw <- read.table(paste0(obs_dir,"ObservedData_tribs.csv"),
                      sep=",",dec=".",skip=0,header=TRUE, 
                      stringsAsFactors = FALSE, na.strings = "NA")


flow.key <- read.table(paste0(obs_dir,"inflow_key.csv"),
                      sep=",",dec=".",skip=0,header=TRUE, 
                      stringsAsFactors = FALSE, na.strings = "NA")

# clean up accrection df
colnames(acc) <- c("stream_id","node_id","Stream_km","acc_flow","temperature","outflow")
acc[is.na(acc$acc_flow),c("acc_flow")] <- 0
acc[is.na(acc$outflow),c("outflow")] <- 0

flow.cols <- grep("flow", tolower(colnames(inflow)), value = FALSE)

#flows <- inflow[,c(1,flow.cols)]
flows.n <- inflow.n[,c(1,flow.cols)]

# wide to long input flows
flows.l <- melt(flows.p, id.vars =c("datetime"),variable_name="inflow_id")
colnames(flows.l) <- c("Datetime","inflow_id", "ccc_flow")
flows.l$inflow_id <- as.numeric(gsub(pattern="Flow", replacement="", flows.l$inflow_id, ignore.case = FALSE,fixed = FALSE))
flows.l <- merge(flows.l,flow.key,by="inflow_id", all.x=TRUE)

# wide to long natural flows
flows.n.l <- melt(flows.n, id.vars =c("datetime"),variable_name="inflow_id")
colnames(flows.n.l) <- c("Datetime","inflow_id", "nat_flow")
flows.n.l$inflow_id <- as.numeric(gsub(pattern="Flow", replacement="", flows.n.l$inflow_id, ignore.case = FALSE,fixed = FALSE))
flows.n.l <- merge(flows.n.l,flow.key,by="inflow_id", all.x=TRUE)


#flows.all <- merge(flows.l,flows.n.l[,c(2:4)],by=c("Stream_km","datetime"))

# sum the current flow in the downstream direction by date
flows.m <- t(data.matrix(flows.p[,c(2:32)], rownames.force = NA))
flows.sum <- t(apply(flows.m, 2, cumsum))

# sum the natural flow in the downstream direction by date
flows.m <- t(data.matrix(flows.n[,c(2:32)], rownames.force = NA))
flows.n.sum <- t(apply(flows.m, 2, cumsum))

# Format Datetime
inflow$datetime <- as.POSIXct(inflow$datetime,format="%m/%d/%Y %H:%M") #m/dd/yyyy hh:mm
obs.raw$Datetime <- as.POSIXct(obs.raw$Datetime,format="%m/%d/%Y %H:%M") #m/dd/yyyy hh:mm
flows.l$Datetime <- as.POSIXct(flows.l$Datetime,format="%m/%d/%Y %H:%M") #m/dd/yyyy hh:mm
flows.n.l$Datetime <- as.POSIXct(flows.n.l$Datetime,format="%m/%d/%Y %H:%M") #m/dd/yyyy hh:mm

flows.in <- flows.l[flows.l$variable=="inflow",]
flows.out <- flows.l[flows.l$variable=="withdrawal",]

# Get the observed temperature data
obs <- obs.raw[obs.raw$ParameterCode == "STEMPH",]

# Only use observed data for the same time period as the model
Tstart <-min(inflow$datetime)
Tend <-max(inflow$datetime)

obs <- obs[obs$Datetime >= Tstart & obs$Datetime <= Tend,]

# Drop these cols
obs$Notes <- NULL
obs$ParameterCode <- NULL
obs$ParameterName <- NULL

p.temps <- ggplot(data=obs, aes(x=Datetime)) +
  geom_line(aes(y=value)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Date") +
  ylab("Hourly Temperature (C)") +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y") +
  scale_y_continuous(limits = c(5, 20)) +
  facet_wrap(~SiteName, ncol=2)
p.temps

p.flow.in <- ggplot(data=flows.in, aes(x=Datetime)) +
  geom_line(aes(y=ccc_flow*35.3147)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Date") +
  ylab("Flow Rate (cfs)") +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y") +
  facet_wrap(~sim_name, ncol=3)
p.flow.in

ggsave(file=paste0(out_dir,sim_name,"_TRIB_FLOWS.png"),
       plot=p.flow.in,
       height=9.5,
       width=6.75,
       units="in")

p.flow.out <- ggplot(data=flows.out, aes(x=Datetime)) +
  geom_line(aes(y=ccc_flow*35.3147)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Date") +
  ylab("Withdrawal Rate (cfs)") +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%m/%Y") +
  facet_wrap(~sim_name, ncol=3)
p.flow.out

ggsave(file=paste0(out_dir,sim_name,"_WITHDRAWALS.png"),
       plot=p.flow.out,
       height=9.5,
       width=6.75,
       units="in")

p.acc <- ggplot(data=acc, aes(x=Stream_km)) +
  geom_line(aes(y=acc_flow*35.3147)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Accretion Flow (cfs)")
p.acc

ggsave(file=paste0(out_dir,sim_name,"_ACCRECTION.png"),
       plot=p.acc,
       height=3,
       width=6.75,
       units="in")
