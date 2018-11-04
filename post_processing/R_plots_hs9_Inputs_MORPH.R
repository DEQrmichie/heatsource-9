#---------------------------------------------------------------------------------------
# Plot heat source 9 morphology model input parameters including:

# Channel Elevation (meters)
# Channel Gradient
# Channel Bottom Width (meters)
# Channel Angle Ration (z)
# Manning's n
# Sediment Thermal conducutity
# Sediment Thickness
# Percent Hyporheic
# Porosity
#---------------------------------------------------------------------------------------

library (reshape2)
library (ggplot2)

sim_name <- "Yach_s1_16"

out_dir <-"F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_16/inputs_as_csv/"
sim_dir <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/s1_16/inputs_as_csv/"
#obs_path <- "F:/WorkSpace/Mid_Coast/Heat_Source/Yachats/sim_01_CCC/"


morph <- read.table(paste0(sim_dir,"morph.csv"),
                    sep=",",dec=".",skip=0,header=TRUE, 
                    stringsAsFactors = FALSE, na.strings = "NA")

colnames(morph) <- c("stream_id","node_id","stream_km","elevation","gradient",
                     "bottom_width", "channel_angle_z", "mannings_n", 
                     "sed_thermal_conductivity", "sed_thermal_diffusivity", 
                     "sed_hyporheic_thicknesss", "hyporheic_percent", "porosity")


#---------------------------------------------------------------------------------------
# Channel Elevation (meters)

p.1 <- ggplot(data=morph, aes(x=stream_km)) +
  geom_line(aes(y=elevation)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Elevation (meters)")

ggsave(file=paste0(out_dir,sim_name,"_ELEVATION.png"),
       plot=p.1,
       height=3,
       width=6.75,
       units="in")

#---------------------------------------------------------------------------------------
# Channel Gradient

p.2 <- ggplot(data=morph, aes(x=stream_km)) +
  geom_line(aes(y=gradient)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Gradient")

ggsave(file=paste0(out_dir,sim_name,"_GRADIENT.png"),
       plot=p.2,
       height=3,
       width=6.75,
       units="in")

#---------------------------------------------------------------------------------------
# Channel Bottom Width (meters)

p.3 <- ggplot(data=morph, aes(x=stream_km)) +
  geom_line(aes(y=bottom_width)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Bottom Widith (meters)")

ggsave(file=paste0(out_dir,sim_name,"_WIDTH_BOTTOM.png"),
       plot=p.3,
       height=3,
       width=6.75,
       units="in")
#---------------------------------------------------------------------------------------
# Channel Angle Ratio (z)

p.4 <- ggplot(data=morph, aes(x=stream_km)) +
  geom_line(aes(y=channel_angle_z)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Channel Angle Z")

ggsave(file=paste0(out_dir,sim_name,"_CHANNEL_ANGLE_Z.png"),
       plot=p.4,
       height=3,
       width=6.75,
       units="in")


#---------------------------------------------------------------------------------------
# Manning's n

p.5 <- ggplot(data=morph, aes(x=stream_km)) +
  geom_line(aes(y=mannings_n)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Manning's n")

ggsave(file=paste0(out_dir,sim_name,"_MANNINGS.png"),
       plot=p.5,
       height=3,
       width=6.75,
       units="in")

#---------------------------------------------------------------------------------------
# Sediment Thermal conducutity

p.6 <- ggplot(data=morph, aes(x=stream_km)) +
  geom_line(aes(y=sed_thermal_conductivity)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Sed. Thermal Conductivity (W/m/*C)") +
  ylim(0,5)

ggsave(file=paste0(out_dir,sim_name,"_SED_THERMAL_COND.png"),
       plot=p.6,
       height=3,
       width=6.75,
       units="in")
#---------------------------------------------------------------------------------------
# Sediment Thickness

p.7 <- ggplot(data=morph, aes(x=stream_km)) +
  geom_line(aes(y=sed_hyporheic_thicknesss)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Thickness (meters)") +
  ylim(0,2)

ggsave(file=paste0(out_dir,sim_name,"_SED_THICKNESS.png"),
       plot=p.7,
       height=3,
       width=6.75,
       units="in")

#---------------------------------------------------------------------------------------
# Percent Hyporheic

p.8 <- ggplot(data=morph, aes(x=stream_km)) +
  geom_line(aes(y=hyporheic_percent)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Hyporheic fraction of total flow") +
  ylim(0,0.02)

ggsave(file=paste0(out_dir,sim_name,"_HYPORHEIC.png"),
       plot=p.8,
       height=3,
       width=6.75,
       units="in")
#---------------------------------------------------------------------------------------
# Porosity

p.9 <- ggplot(data=morph, aes(x=stream_km)) +
  geom_line(aes(y=porosity)) +
  guides(color=guide_legend(override.aes=list(shape=c(16),linetype=c(1))))  +
  scale_colour_manual(values=c("Black")) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  xlab("Model Stream Kilometer") +
  ylab("Porosity") +
  ylim(0,1)

ggsave(file=paste0(out_dir,sim_name,"_POROSITY.png"),
       plot=p.9,
       height=3,
       width=6.75,
       units="in")
