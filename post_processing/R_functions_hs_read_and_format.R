#---------------------------------------------------------------------------------------
# Functions to read and format heat source outputs
#---------------------------------------------------------------------------------------

read.hs7.outputs <- function(output_dir, file_name, sheet_name) {
  # Function to read any hourly output from heat source 7. Returns the data
  # as a dataframe. Excel workbook needs to be saved as .xlsx. .xls do not seem to work.
  
  library(readxl)
  
  excel.data <- read_excel(path=paste0(output_dir,"/",file_name), sheet=sheet_name, skip=14, na = c("","N/A", " "))
  
  new_cols <- excel.data[[4]]
  hs7.data <- data.frame(t(excel.data[,c(5:ncol(excel.data))]))
  colnames(hs7.data) <- new_cols
  
  hs7.data$Datetime <- as.numeric(rownames(hs7.data))
  rownames(hs7.data) <-NULL
  
  return(hs7.data)
}

read.hs8.outputs <- function(output_dir, file_name, sheet_name=NULL) {
  # Function to read any hourly output from heat source 8. Return the data
  # as a dataframe.
  
  # remove the extension in case it was added
  base_name <- gsub("\\..*","",file_name)
  
  hs8.data <- read.table(paste0(output_dir, base_name,".txt"),
                         sep="",dec=".",skip=2,header=TRUE, 
                         stringsAsFactors = FALSE, na.strings = "NA")
  return(hs8.data)
  
}

read.hs9.outputs <- function(output_dir, file_name, sheet_name=NULL) {
  # Function to reads any hourly output from heat source 9. Returns the data
  # as a dataframe.
  
  # remove the extension in case it was added
  base_name <- gsub("\\..*","", file_name)
  
  hs9.data <- read.table(paste0(output_dir, base_name,".csv"),
                         sep=",",dec=".",skip=6,header=TRUE, 
                         stringsAsFactors = FALSE, na.strings = "NA")
  return(hs9.data)
}

read.flux.outputs <-function(output_dir, sim_name, hs_ver=8) {
  # Reads all the hourly flux output from heat source, does some formatting, calculates total flux by
  # summing all the fluxes, and returns the data as a dataframe in long format. 
  # Simulation name and flux constituents are added as an ID variable.
  
  library(reshape2)
  library(lubridate)
  
  # Assign the correct read function based on model version
  if (as.integer(hs_ver) == 7) {
    stop("Reading flux outputs for heat source version 7 has not been implemented yet. Sorry.")
  }
    
  if (as.integer(hs_ver) == 8) {
    read.hs.outputs = read.hs8.outputs
    long <- "Heat_TR"
  }
  
  if (as.integer(hs_ver) == 9) {
    read.hs.outputs = read.hs9.outputs
    long <- "Heat_Long"
  }
  
  # read the data
  flux.cond <- read.hs.outputs(output_dir,"Heat_Cond")
  flux.conv <- read.hs.outputs(output_dir,"Heat_Conv")
  flux.evap <- read.hs.outputs(output_dir,"Heat_Evap")
  flux.long <- read.hs.outputs(output_dir,long)
  flux.sr4  <- read.hs.outputs(output_dir,"Heat_SR4")
  
  flux.cond$constituent <- "Bed Conduction"
  flux.conv$constituent <- "Air Convection"
  flux.evap$constituent <- "Evaporation"
  flux.long$constituent <- "Longwave"
  flux.sr4$constituent <- "Solar Radiation"
  
  flux.raw <- rbind(flux.cond,flux.conv,flux.evap,flux.long,flux.sr4)
  
  flux.raw$sim <- sim_name
  
  # Convert data from wide to long
  flux.l <- melt(flux.raw, id.vars =c("Datetime","sim", "constituent"),variable.name=c("Stream_km"))
  colnames(flux.l) <- c("Datetime","sim", "constituent","Stream_km","value")
  
  flux.l$Stream_km <- as.numeric(gsub(pattern="X", replacement="", 
                                      flux.l$Stream_km, ignore.case = FALSE,fixed = FALSE))
  
  # Convert back to wide to sum total flux
  flux.w  <- reshape(flux.l, timevar="constituent", idvar = c("Stream_km","sim", "Datetime"), 
                     direction="wide")
  
  colnames(flux.w) <- c("Datetime","sim", "Stream_km",
                        "Bed Conduction","Air Convection",
                        "Evaporation","Longwave","Solar Radiation")
  
  flux.w$Total <- rowSums(flux.w[,4:8])
  
  # now back to long format
  flux.l <- melt(flux.w, id.vars =c("Datetime","sim", "Stream_km"))
  colnames(flux.l) <- c("Datetime", "sim", "Stream_km","constituent","value")
  
  flux.l$Datetime <-round_date(as.POSIXct(( flux.l$Datetime*60*60*24), origin="1899-12-30", tz="GMT"), unit = "minute")
  flux.l$Date <- format(flux.l$Datetime,"%m/%d/%Y")
  flux.l$hour <-as.integer(format(flux.l$Datetime, "%H"))
  
  return(flux.l)
  
}

read.hs7.shade <- function(output_dir, file_name, sim_name, constituent_name="Effective Shade", 
                           statistic_name="Percent", sheet_name="Chart-Shade") {
  # Function to read effective shade output from heat source 7. Returns the data
  # as a dataframe. Excel workbook needs to be saved as .xlsx. .xls do not seem to work.
  
  library(readxl)
  library(lubridate)
  
  excel.data <- read_excel(path=paste0(output_dir,"/",file_name), sheet=sheet_name, skip=12, na = c("","N/A", " "),
                           col_names=c("Stream_km", "Datetime", "value"),
                           col_types =c("numeric","numeric","numeric"))
  
  excel.data$constituent <- constituent_name
  excel.data$statistic <- statistic_name
  excel.data$sim <- sim_name
  
  excel.data$Datetime <- as.numeric(excel.data$Datetime)
  
  excel.data$Datetime <-round_date(as.POSIXct((excel.data$Datetime*60*60*24), origin="1899-12-30", tz="GMT"), unit = "minute")
  
  excel.data$Date <- format(excel.data$Datetime,"%m/%d/%Y")
  
  excel.data$hour <-as.integer(format(excel.data$Datetime, "%H"))
  
  return(excel.data)
  
}

read.solar.flux.outputs <-function(output_dir, sim_name, hs_ver=8) {
  # Reads all the hourly solar flux outputs from heat source, does some formatting, and returns the data as 
  # a dataframe in long format. 
  # Simulation name and flux constituents are added as an ID variable.
  
  library(reshape2)
  library(lubridate)
  
  # Assign the correct read function based on model version and   
  # read the data
  
  if (as.integer(hs_ver) == 7) {
    stop("Reading solar flux outputs for heat source version 7 has not been implemented yet. Sorry.")
  }
  
  
  if (as.integer(hs_ver) == 8) {
    read.hs.outputs = read.hs8.outputs
    
    # read the data
    flux.sr1 <- read.hs.outputs(output_dir,"Heat_SR1")
    flux.sr4 <- read.hs.outputs(output_dir,"Heat_SR4")
    flux.sr6 <- read.hs.outputs(output_dir,"Heat_SR6")
    
    flux.sr1$constituent <- "SR1: Solar Radiation above Topo"
    flux.sr4$constituent <- "SR4: Solar Radiation above Stream"
    flux.sr6$constituent <- "SR6: Solar Radiation received by Stream"
    
    flux.raw <- rbind(flux.sr1, flux.sr4, flux.sr6)
    
  }
  
  if (as.integer(hs_ver) == 9) {
    read.hs.outputs = read.hs9.outputs
    
    # read the data
    flux.sr1 <- read.hs.outputs(output_dir,"Heat_SR1")
    flux.sr4 <- read.hs.outputs(output_dir,"Heat_SR4")
    flux.sr6 <- read.hs.outputs(output_dir,"Heat_SR6")
    
    flux.sr1$constituent <- "SR1: Solar Radiation above Topo"
    flux.sr4$constituent <- "SR4: Solar Radiation above Stream"
    flux.sr6$constituent <- "SR6: Solar Radiation received by Stream"
    
    # These are availibele in hs9 but not in hs8 so they are left out
    #flux.sr2 <- read.hs.outputs(output_dir,"Heat_SR2")
    #flux.sr3 <- read.hs.outputs(output_dir,"Heat_SR3")
    #flux.sr5 <- read.hs.outputs(output_dir,"Heat_SR5")
    
    #flux.sr2$constituent <- "SR2: Solar Radiation below Topo"
    #flux.sr3$constituent <- "SR3: Solar Radiation below Land Cover"
    #flux.sr5$constituent <- "SR5: Solar Radiation entering Stream"
    
    flux.raw <- rbind(flux.sr1, flux.sr4, flux.sr6)
    #flux.raw <- rbind(flux.sr1, flux.sr2, flux.sr3, flux.sr4, flux.sr5, flux.sr6)
    
  }
  
  flux.raw$sim <- sim_name
  
  # Convert data from wide to long
  flux.l <- melt(flux.raw, id.vars =c("Datetime","sim", "constituent"),variable.name=c("Stream_km"))
  colnames(flux.l) <- c("Datetime","sim", "constituent","Stream_km","value")
  
  flux.l$Stream_km <- as.numeric(gsub(pattern="X", replacement="", 
                                      flux.l$Stream_km, ignore.case = FALSE,fixed = FALSE))
  
  # Convert back to wide to sum total flux
  #flux.w  <- reshape(flux.l, timevar="constituent", idvar = c("Stream_km","sim", "Datetime"), 
  #                   direction="wide")
  
  #colnames(flux.w) <- c("Datetime","sim", "Stream_km",
  #                      "Solar Radiation above Topo","Solar Radiation above Topo",
  #                      "Solar Radiation received by Stream")
  
  # now back to long format
  #flux.l <- melt(flux.w, id.vars =c("Datetime","sim", "Stream_km"))
  #colnames(flux.l) <- c("Datetime", "sim", "Stream_km","constituent","value")
  
  flux.l$Datetime <-round_date(as.POSIXct((flux.l$Datetime*60*60*24), origin="1899-12-30", tz="GMT"), unit = "minute")
  flux.l$Date <- format(flux.l$Datetime,"%m/%d/%Y")
  flux.l$hour <-as.integer(format(flux.l$Datetime, "%H"))
  
  return(flux.l)
  
}

read.hs.outputs <- function(output_dir, file_name, constituent_name, statistic_name, sim_name, hs_ver=9, sheet_name=NULL) {
  # Reads any output from heat source version 7-9. Does some formatting, and returns the data
  # as a dataframe in long format. 
  # Simulation name, constituent, statistic are strings
  # hours are added as ID variables.
  # sheet_name is only used for heat source 7
  
  library(reshape2)
  library(lubridate)
  
  # Assign the correct read function based on model version
  if (as.integer(hs_ver) == 7) {
    read.hs.outputs = read.hs7.outputs
    }
  
  if (as.integer(hs_ver) == 8) {
    read.hs.outputs = read.hs8.outputs
  }
  
  if (as.integer(hs_ver) == 9) {
    read.hs.outputs = read.hs9.outputs
  }
  
  # read the data
  data.raw <- read.hs.outputs(output_dir, file_name, sheet_name)
  
  data.raw$constituent <- constituent_name
  data.raw$statistic <- statistic_name
  data.raw$sim <- sim_name
  
  # Convert data from wide to long
  data.l <- melt(data.raw, id.vars =c("Datetime","sim", "constituent","statistic"),variable.name=c("Stream_km"))
  colnames(data.l) <- c("Datetime","sim", "constituent","statistic","Stream_km","value")
  
  data.l$Stream_km <- as.numeric(gsub(pattern="X", replacement="", 
                                      data.l$Stream_km, ignore.case = FALSE,fixed = FALSE))
  
 
  data.l$Datetime <-round_date(as.POSIXct((data.l$Datetime*60*60*24), origin="1899-12-30", tz="GMT"), unit = "minute")
  
  data.l$Date <- format(data.l$Datetime,"%m/%d/%Y")
  
  data.l$hour <-as.integer(format(data.l$Datetime, "%H"))
  
  return(data.l)
  
}

read.hs.inputs <- function (input_dir, file_name) {
  
  data.raw <- read.table(paste0(input_dir,file_name),
                         sep=",",dec=".",skip=0,header=TRUE, 
                         stringsAsFactors = FALSE, na.strings = "NA")
  
  return(data.raw)
}

read.obs <- function(obs_dir, file_name) {
  # Read observation data
  
  obs.raw <- read.table(paste0(obs_dir,file_name),
                        sep=",", dec=".",skip=0, header=TRUE, 
                        stringsAsFactors = FALSE, na.strings = "NA")#,
                        #quote="")
  
  # Format Datetime
  obs.raw$Datetime <- as.POSIXct(obs.raw$Datetime,format="%m/%d/%Y %H:%M", tz="GMT") #6/17/2003 0:00
  
  # Add character Date
  obs.raw$Date <- format(obs.raw$Datetime, "%m/%d/%Y")
  
  return(obs.raw)
  
}

obs2simkm <- function(obs_dir, file_name, constituentCode, simkm) {
  # Returns a dataframe matching the observation site to the closest model 
  # simulation kilometer. 
  # simkm is a vector of unique model simulation kilometers
  
  # Read obs data
  obs.raw <- read.obs(obs_dir=obs_dir, file_name=file_name)
  
  # get the constituent
  obs <- obs.raw[obs.raw$constituentCode == constituentCode,]
  
  # Drop these cols
  obs$Notes <- NULL
  obs$constituentCode <- NULL
  obs$constituent <- NULL
  obs$statistic <- NULL
  
  # just get the sites and stream_km
  obs.wide  <- reshape(obs, timevar= "Datetime", idvar = c("Stream_km","SiteName"), 
                       direction="wide", drop = c("Date","value","group"))
  
  # get unique obs km
  obs.km <- as.numeric(unique(obs$Stream_km))
  
  sim.km <- numeric(0)
  
  for (i in 1:length(obs.km)) {
    # this makes a vector of the model km that is closest to the obs km
    sim.km[i] <-simkm[which(abs(simkm-obs.km[i])==min(abs(simkm-obs.km[i])))[1]]
  }
  
  km <- cbind(data.frame(obs.km),data.frame(sim.km))
  km <- merge(km,obs.wide,by.x="obs.km",by.y="Stream_km")
  
  return(km)
  
}

calc.7dadm <- function(df) {
  # Returns a dataframe of the daily maximum and 7DADM temperatures 
  # by date in long format. Constituent and statistic names are added as an ID variable.
  
  # df input is a dataframe. df must have the following columns names:
  # Datetime (with datetime as POSIXlt)
  # sim (with simulation name in character format)
  # Stream_km (as numeric)
  # value (numeric hourly stream temperature)       
  
  
  library(zoo)
  library(reshape2)
  
  # Add character Date
  df$Date <- format(df$Datetime, "%m/%d/%Y")
  
  data.max <- as.data.frame(aggregate(df$value,by=list(df$Date,df$Stream_km, df$sim), 
                                      FUN=max, na.rm = TRUE))
  
  colnames(data.max) <- c("Date","Stream_km", "sim", "max")
  
  # sort by Stream_km and Date
  data.max <- data.max[with(data.max, order(sim,-Stream_km,Date)), ]
  
  # change infinite values to NA
  data.max$max <- ifelse(is.infinite(data.max$max),NA,data.max$max)
  
  # Calculate the 7 day max rolling average
  data.max$sdadm <- ave(data.max$max, data.max$Stream_km, FUN = 
                          function(x) rollapply(zoo(x), 7, mean, fill = NA, align = "right")) 
  
  # convert to long format
  data.l <- melt(data.max, id.vars =c("Date","sim","Stream_km"))
  
  colnames(data.l) <- c("Date","sim","Stream_km","statistic","value")
  
  # make sure it is sorted
  data.l <- data.l[with(data.l, order(sim,statistic,-Stream_km,Date)), ]
  
  data.l$constituent <- unique(df$constituent)
  
  # rename the statistics
  data.l$statistic <- gsub(pattern="max", replacement="Daily Maximum Temperature", 
                             data.l$statistic, ignore.case = FALSE, fixed = FALSE)
  
  data.l$statistic <- gsub(pattern="sdadm", replacement="7DADM Temperature", 
                             data.l$statistic, ignore.case = FALSE,fixed = FALSE)
  
  
  
  return(data.l)
}

calc.7dadm2 <- function(output_dir, sim_name, hs_ver=9, file_name="Temp_H2O", sheet_name=NULL) {
  # Same function as calc.7DADM but instad of passing a dataframe as the input this function 
  # reads the the model hourly temperature output using read.hs.outputs().
  # Function does some formatting, and returns a dataframe of the hourly, daily maximum, and 7DADM temperatures 
  # by date in long format. Simulation name and statistics are added as an ID variable.
  
  library(zoo)
  library(reshape2)
  
  data.raw <- read.hs.outputs(output_dir=output_dir, 
                              file_name=file_name, 
                              constituent_name="Temperature",
                              statistic_name = "hourly",
                              sim_name=sim_name, 
                              hs_ver=hs_ver,
                              sheet_name=sheet_name)
  
  # Add character Date
  data.raw$Date <- format(data.raw$Datetime, "%m/%d/%Y")
  
  data.max <- as.data.frame(aggregate(data.raw$value,by=list(data.raw$Date,data.raw$Stream_km), 
                                      FUN=max, na.rm = TRUE))
  
  colnames(data.max) <- c("Date","Stream_km", "max")
  
  data.max$sim <- sim_name
  
  # sort by Stream_km and Date
  data.max <- data.max[with(data.max, order(-Stream_km,Date)), ]
  
  # Calculate the 7 day max rolling average
  data.max$sdadm <- ave(data.max$max, data.max$Stream_km, FUN = 
                          function(x) rollapply(zoo(x), 7, mean, fill = NA, align = "right")) 
  
  # convert to long format
  data.l <- melt(data.max, id.vars =c("Date","sim","Stream_km"))
  
  colnames(data.l) <- c("Date","sim","Stream_km","statistic","value")
  
  data.l$constituent <- constituent_name
  
  # rename the stats
  data.l$statistic <- gsub(pattern="max", replacement="Daily Maximum Temperature", 
                             data.l$statistic, ignore.case = FALSE, fixed = FALSE)
  
  data.l$statistic <- gsub(pattern="sdadm", replacement="7DADM Temperature", 
                             data.l$statistic, ignore.case = FALSE,fixed = FALSE)
  
  return(data.l)
}

calc.summary <- function(df) {
  # Function to calculate maximum, median, and minimum values for all dates by sim, constituent, statistic, and stream_km,
  # df input must be a dataframe in long format with the following column names:
  # Date (as character), 
  # Stream_km (as numeric), 
  # sim (as string), 
  # constituent (as string),
  # statistic (as string)
  # value (as numeric)
  # Returns a dataframe
  
  data.summary <- as.data.frame(as.list(aggregate(df$value, by=c(list(df$Stream_km), 
                                                                 list(df$sim),
                                                                 list(df$constituent),
                                                                 list(df$statistic)), 
                                                  FUN = function(x) c(mn = min(x,na.rm=TRUE),
                                                                      med = median(x,na.rm=TRUE), 
                                                                      mx = max(x,na.rm=TRUE)))))
  # fix colnames
  colnames(data.summary) <- c("Stream_km", "sim", "constituent", "statistic", "min", "median","max")
  
  # undo factors
  data.summary$sim <- as.character(data.summary$sim)
  data.summary$constituent <- as.character(data.summary$constituent)
  data.summary$statistic <- as.character(data.summary$statistic)
  
  return(data.summary)
  
}