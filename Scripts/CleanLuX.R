

light_cleanup<-function (data.path, light.serial, output.path, tf_write = FALSE,  recursive_tf = FALSE) {
  file.names.Cal<-here(data.path,paste0(light.serial,".csv"))
  
  LightLog <- read_csv(file.names.Cal, skip = 1, col_names = FALSE) %>%
    rename("Date" = "X2",
           "TempInSitu" = "X3",
           "Lux" = "X4")%>%
    select(!X1)%>%
    mutate(Date = mdy_hms(Date))
  
   
  return(LightLog)
}


############
#devtools::install_github("dbarnas/mooreasgd") # if package has updated since last run
library(tidyverse)
library(lubridate)
library(gsw)
library(here)
library(gridExtra)
library(mooreasgd)

###################################
### File Paths
###################################

### Input
# Path to folder storing logger .csv files
path.log<-here("Data","LoggerData","Raw") # Logger in situ file path (CT and Water Level files)
file.date <- "2023-03-20" # logger date used in file name(s)

### Output
# Path to store logger files
path.output<-here("Data","LoggerData","Processed") # Output file path


###################################
### Logger Serial Numbers
###################################

Light_Serial <- "12"

###################################
### Logger Launch and Retrieval dates
###################################

# Log dates
start.date <- ymd('2023-10-24')
end.date <- ymd('2023-11-03')

###################################
### Import calibration and launch records
###################################

# Read in files that are updated with calibration and launch information
#calibration.log<-read_csv(here("Data","Tris_Calibration.csv")) # Calibration time logs
launch.log<-read_csv(here("Data","Launch_Log.csv")) %>%  # Launch time logs
  filter(Log_Type == "pendent")


############################################################
### Read in Logger Files
############################################################

# cleanup function pulled from 'mooreasgd' package
# Reads in raw csv and returns tidied csv for the probe with the specified serial number

# In Situ pH file
Light.data <- light_cleanup(data.path = path.log, light.serial = Light_Serial) 

## parse the date ####
launch.log <- launch.log %>%
  mutate(time_start = mdy_hm(time_start), # convert to time
         time_end = mdy_hm(time_end),
         start  = date(time_start), # extract the date
         end = date(time_end)) %>%
  filter(CowTag == Light_Serial, # pull out the right serial number
         start == ymd(start.date),
         end == ymd(end.date))
############################################################
### In Situ Logger Data
############################################################
LightLog<-Light.data %>% # extract the data you need
  filter(between(Date,launch.log$time_start, launch.log$time_end))%>%
  select(Date, TempInSitu, Lux)



LightLog %>% 
    ggplot(aes(x = Date, y = Lux, color = TempInSitu)) + 
    geom_line() + 
    theme_bw() +
    labs(x = "Date", color = "Temperature (C)") +
    ggtitle(paste("Lux",Light_Serial))
  
  
# write out the clean data
write_csv(LightLog, paste0(path.output,"/Calibrated_",Light_Serial,".csv"))


#### Temperatures #####

temp_cleanup<-function (data.path, light.serial, output.path) {
  file.names.Cal<-here(data.path,paste0(light.serial,".csv"))
  
  TempLog <- read_csv(file.names.Cal, skip = 2, col_names = FALSE) %>%
    rename("Date" = "X2",
           "TempInSitu" = "X3")%>%
    select(!X1)%>%
    mutate(Date = mdy_hms(Date))
  
  return(TempLog)
}  

Light_Serial <- "47"

###################################
### Logger Launch and Retrieval dates
###################################

# Log dates
start.date <- ymd('2023-10-25')
end.date <- ymd('2023-11-02')

###################################
### Import calibration and launch records
###################################

# Read in files that are updated with calibration and launch information
#calibration.log<-read_csv(here("Data","Tris_Calibration.csv")) # Calibration time logs
launch.log<-read_csv(here("Data","Launch_Log.csv")) %>%  # Launch time logs
  filter(Log_Type == "temp")


############################################################
### Read in Logger Files
############################################################

# cleanup function pulled from 'mooreasgd' package
# Reads in raw csv and returns tidied csv for the probe with the specified serial number

# In Situ pH file
Temp.data <- temp_cleanup(data.path = path.log, light.serial = Light_Serial) %>%
  select(Date, TempInSitu)

## parse the date ####
launch.log <- launch.log %>%
  mutate(time_start = mdy_hm(time_start), # convert to time
         time_end = mdy_hm(time_end),
         start  = date(time_start), # extract the date
         end = date(time_end)) %>%
  filter(CowTag == Light_Serial, # pull out the right serial number
         start == ymd(start.date),
         end == ymd(end.date))
############################################################
### In Situ Logger Data
############################################################
TempLog<-Temp.data %>% # extract the data you need
  filter(between(Date,launch.log$time_start, launch.log$time_end))%>%
  select(Date, TempInSitu)



TempLog %>% 
  ggplot(aes(x = Date, y = TempInSitu)) + 
  geom_line() + 
  theme_bw() +
  labs(x = "Date", y = "Temperature (C)") 


# write out the clean data
write_csv(TempLog, paste0(path.output,"/Calibrated_",Light_Serial,".csv"))

