#---------------------------------------------------------------------------------------
# Plot heat source 7 land cover model input parameters including:

# Channel Width (m)
# Topographic Angles
# Mean Landcover height (m)
# Mean Landcover density (%)
#---------------------------------------------------------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)

name <- "Antelope Creek"

# plot height and width
h <- 4
w <- 6.75

sim1_name <- "Current Condition"
sim1_dir <- "T:/TMDL_ER/Klamath_and_Lost_Rivers/Heatsource/Lost_River_Tribs/Antelope/"
sim1_file <- "HS7.Antelope.Current.Veg.xlsm"

sheet_name = "TTools Data"

fun_dir <- "E:/GitHub/heatsource-9/post_processing"
out_dir <- "T:/TMDL_ER/Klamath_and_Lost_Rivers/Heatsource/Lost_River_Tribs/Antelope/R/"

source(paste0(fun_dir,"/","R_functions_hs_read_and_format.R"))


#--  Read data ------------------------------------

lcdata.raw <- read.hs7.landcover(output_dir=sim1_dir, file_name=sim1_file, sim_name=sim1_name, sheet_name=sheet_name)


#---------------------------------------------------------------------------------------
# Channel Width

lcdata.width <- lcdata.raw %>%
  dplyr::select(Stream_km, width) %>%
  tidyr::gather(key="legend",value="value", -Stream_km) %>%
  dplyr::mutate(Stream_km=as.numeric(Stream_km),
                value=as.numeric(value))

p.wid <- lcdata.width %>%
  ggplot(aes(x=Stream_km, y=value)) +
  geom_line(size=1) +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        panel.background = element_rect(fill="white", colour = "black"),
        strip.background =element_rect(fill="white", colour = "black"),
        panel.grid.major =element_blank(),
        plot.title = element_text(size=12, hjust = 0.5)) +
  xlab("Distance from Willow Valley Reservoir (Kilometers)") +
  ylab("Active Channel Width (m)") + 
  ylim(0,NA)
p.wid

ggsave(file=paste0(out_dir,name,"_LC_ChanWidth_",sim1_name,".png"),
       plot=p.wid,
       height=h,
       width=w,
       units="in")

#---------------------------------------------------------------------------------------
# Topo

lcdata.topo <- lcdata.raw %>%
  dplyr::select(Stream_km, topo_w, topo_s, topo_e) %>%
  tidyr::gather(key="legend",value="value", -Stream_km) %>%
  dplyr::mutate(Stream_km=as.numeric(Stream_km),
                value=as.numeric(value))

p.topo <- lcdata.topo %>%
  ggplot(aes(x=Stream_km, y=value)) +
  geom_line(size=1) +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        panel.background = element_rect(fill="white", colour = "black"),
        strip.background =element_rect(fill="white", colour = "black"),
        panel.grid.major =element_blank(),
        plot.title = element_text(size=12, hjust = 0.5)) +
  xlab("Distance from Willow Valley Reservoir (Kilometers)") +
  ylab("Topographic Shade Angle (degrees)") + 
  facet_wrap(~legend, ncol = 1, labeller = as_labeller(c("topo_w"="West", "topo_s"="South", "topo_e"="East")))
p.topo

ggsave(file=paste0(out_dir,name,"_Topo_",sim1_name,".png"),
       plot=p.topo,
       height=h,
       width=w,
       units="in")

#---------------------------------------------------------------------------------------
# Height

lcdata.ht <- lcdata.raw %>%
  dplyr::select(Stream_km,height_l, height_r) %>%
  tidyr::gather(key="legend",value="value", -Stream_km) %>%
  dplyr::mutate(Stream_km=as.numeric(Stream_km),
                value=as.numeric(value))
    
  # longitudinal plot of ES for Heat Source
p.ht <- lcdata.ht %>%
  ggplot(aes(x=Stream_km, y=value)) +
  geom_line(size=1) +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        panel.background = element_rect(fill="white", colour = "black"),
        strip.background =element_rect(fill="white", colour = "black"),
        panel.grid.major =element_blank(),
        plot.title = element_text(size=12, hjust = 0.5)) +
  xlab("Distance from Willow Valley Reservoir (Kilometers)") +
  ylab("Mean Land Cover Height (m)") + 
  ylim(0,35) + 
  facet_wrap(~legend, ncol = 1, labeller = as_labeller(c("height_l"="Left Bank", "height_r"="Right Bank")))
p.ht

ggsave(file=paste0(out_dir,name,"_LC_Height_",sim1_name,".png"),
       plot=p.ht,
       height=h,
       width=w,
       units="in")

#---------------------------------------------------------------------------------------
# Canopy

lcdata.ca <- lcdata.raw %>%
  dplyr::select(Stream_km,density_l, density_r) %>%
  tidyr::gather(key="legend",value="value", -Stream_km) %>%
  dplyr::mutate(Stream_km=as.numeric(Stream_km),
                value=as.numeric(value))

p.ca <- lcdata.ca %>%
  ggplot(aes(x=Stream_km, y=value)) +
  geom_line(size=1) +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        panel.background = element_rect(fill="white", colour = "black"),
        strip.background =element_rect(fill="white", colour = "black"),
        panel.grid.major =element_blank(),
        plot.title = element_text(size=12, hjust = 0.5)) +
  scale_y_continuous(limits=c(0,1), labels = scales::percent_format(accuracy = 1)) +
  xlab("Distance from Willow Valley Reservoir (Kilometers)") +
  ylab("Mean Land Cover Density (%)") +
  facet_wrap(~legend, ncol = 1, labeller = as_labeller(c("density_l"="Left Bank", "density_r"="Right Bank")))
p.ca

ggsave(file=paste0(out_dir,name,"_LC_Density_",sim1_name,".png"),
       plot=p.ca,
       height=h,
       width=w,
       units="in")
