#---------------------------------------------------------------------------------------
# Plot heat source 9 land cover model input parameters including:

# Topographic Angle (meters)
# Lancover height (m)
#---------------------------------------------------------------------------------------

library (reshape2)
library (ggplot2)
library(plyr)

sim_name <- "Yach_s1_16"

out_dir <-"F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_16/inputs_as_csv/"
sim_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_16/inputs_as_csv/"

lcdata.raw <- read.table(paste0(sim_dir,"lcdata.csv"),
                    sep=",",dec=".",skip=0,header=TRUE, 
                    stringsAsFactors = FALSE, na.strings = "NA")

lccodes <- read.table(paste0(sim_dir,"lccodes.csv"),
                     sep=",",dec=".",skip=0,header=TRUE, 
                     stringsAsFactors = FALSE, na.strings = "NA")

#---------------------------------------------------------------------------------------
# Height
lcdata <- lcdata.raw[,c(3, c(grep("LC_", colnames(lcdata.raw))))]

# Convert data from wide to long
lcdata.l <- melt(lcdata, id.vars =c("STREAM_KM"))

# merge landcover height data
lcdata.l <- merge(lcdata.l,lccodes, by.x="value",by.y="CODE")

# Get the transect number from the landcover key
lcdata.l$transect = as.numeric(gsub("_","",substring(lcdata.l$variable,first=5,last=6)))

# mean value by transect
lcdata.ht  <- as.data.frame(aggregate(lcdata.l$HEIGHT,
                                      by=list(lcdata.l$transect,
                                              lcdata.l$STREAM_KM), 
                                      FUN=mean,na.rm = TRUE))

colnames(lcdata.ht)[1] <- "transect"
colnames(lcdata.ht)[2] <- "Stream_km"
colnames(lcdata.ht)[3] <- "value"

# remove zero since we are using zones
lcdata.ht <- lcdata.ht[lcdata.ht$transect > 0,]

lcdata.ht$transect <- factor(lcdata.ht$transect)

# This works for 8 transects only
lcdata.ht$transect <- revalue(lcdata.ht$transect, 
                              c("1"="Northeast",
                                "2"="East",
                                "3"="Southeast",
                                "4"="South",
                                "5"="Southwest",
                                "6"="West",
                                "7"="Northwest",
                                "8"="North"))

p.1 <- ggplot(data=lcdata.ht, aes(x=Stream_km)) +
  geom_line(aes(y=value)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Height (meters)") +
  facet_wrap(~transect, ncol = 1)
p.1

ggsave(file=paste0(out_dir,sim_name,"_HEIGHT.png"),
       plot=p.1,
       height=9.5,
       width=6.75,
       units="in")

#---------------------------------------------------------------------------------------
# Topo

topo <- lcdata.raw[,c("STREAM_KM","TOPO_W", "TOPO_S", "TOPO_E")]

# change col names
colnames(topo) <- c("Stream_km","West", "South", "East")

# Convert data from wide to long
topo.l <- melt(topo, id.vars =c("Stream_km"))


p.2 <- ggplot(data=topo.l, aes(x=Stream_km)) +
  geom_line(aes(y=value)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Topographic Shade Angle (Degrees)") +
  facet_wrap(~variable, ncol = 1)
p.2

ggsave(file=paste0(out_dir,sim_name,"_TOPO.png"),
       plot=p.2,
       height=5,
       width=6.75,
       units="in")
