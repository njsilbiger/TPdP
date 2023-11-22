
#Load Libraries
library(tidyverse)
library(lubridate)
library(here)
library(patchwork)


## read in the metadata
metadata<-read_csv(here("Data","CowTagMetaData.csv"))%>%
  mutate(CowTagID = as.character(CowTagID))

### read in all the data ####
path<-here("Data","LoggerData","Processed")
files <- dir(path = path,pattern = ".csv", full.names = TRUE)

# create a function to read in and clean the files
read_fun<-function(name){
  read_csv({{name}})
   }

AllData<-files %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_fun,.id = "filename") %>%
  mutate(filename = str_split_i(filename, "/",10), # clean the ID
         CowTagID = str_split_i(filename, "_",-1),
         CowTagID = str_split_i(CowTagID,".csv",1))%>%
  left_join(metadata)%>%  # join with the metadata
  select(Site,CowTagID, Date, TempInSitu,Salinity_psu, Lux, Depth_logger, Time_of_Depth, Lat, Long)


## Read in the waterlevel data
path_WL<-here("Data","LoggerData","WL_Processed")
filesWL <- dir(path = path_WL,pattern = ".csv", full.names = TRUE)

WL<-filesWL %>%
  set_names()%>% # set's the id of each list to the file name
  map_df(read_fun,.id = "filename") %>%
  mutate(filename = str_split_i(filename, "/",10), # clean the ID
         Site = str_split_i(filename, "_",-2))%>%
  mutate(Date = mdy_hm(Date)) %>%
  select(Site, Date, Depth_m) %>%
  mutate(Site = ifelse(Site == "LaSource","La Source",Site)) # La Source does not have a space in the name

#### Add the WL data for each site
AllData<-AllData %>%
  left_join(WL, by = c("Site", "Date")) %>%
  filter(between(Date, mdy_hm("10-25-2023 17:00"), mdy_hm("11-02-2023 05:00"))) %>% # make everything the same date for all files
  mutate(CowTagID = factor(as.numeric(CowTagID)))


## Temperature plots
AllData %>%
  # filter(CowTagID != 5 ,# the seep loggers
  #        CowTagID != 41) %>%
  ggplot(aes(x = Date, y = TempInSitu, color = Depth_logger))+
  geom_line()+
  labs(x = "",
       y = expression(paste("Temperature",~degree,"C")),
       color = "Logger Depth (m)")+
  scale_color_continuous(trans = "log")+
  facet_wrap(~Site)+
  theme_bw()

ggsave(here("Output","TempDataAll.png"), width = 8, height = 4)

# Salinity plots
AllData %>%
  # filter(CowTagID != 5 ,# the seep loggers
  #        CowTagID != 41) %>%
  ggplot(aes(x = Date, y = Salinity_psu, color = Depth_logger))+
  geom_line()+
  labs(x = "",
       y = "Salinity (psu)",
       color = "Logger Depth (m)")+
  scale_color_continuous(trans = "log")+
  facet_wrap(~Site)+
  theme_bw()

ggsave(here("Output","SalinityAll.png"), width = 8, height = 4)

### plot depth vs salinity and temperature for logger 41, the lagoon seep

P_lagoon<-AllData %>%
  filter(CowTagID == 41)%>%
  ggplot(aes(x = Depth_m+.1, y = Salinity_psu))+
  geom_point(aes(color = TempInSitu), alpha = 0.5)+
  coord_trans(x = "log", xlim = c(0.01,0.9))+
  geom_smooth(method = lm, formula = y ~ splines::ns(log(x), 4), color = "black")+
  scale_colour_viridis_c(limits = c(24.5,29.6), trans = "log")+
  ylim(18,37)+
  labs(x = "Water Height (m)",
       y = "Salinity (psu)",
       title = "Lagoon",
       color = expression(paste("Temperature",~degree,"C")))+
  theme_bw()+
  theme(legend.position = "bottom",
        title = element_text(hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)
        )

# temp on x
P_Lagoon_Temp<-AllData %>%
  filter(CowTagID == 41)%>%
  ggplot(aes(x = TempInSitu, y = Salinity_psu))+
  geom_point(aes(color = TempInSitu), alpha = 0.5)+
 # coord_trans(x = "log", xlim = c(0.01,0.9))+
  geom_smooth(method = lm, formula = y ~ splines::ns(log(x), 4), color = "black")+
  scale_colour_viridis_c(limits = c(24.5,29.6), trans = "log")+
  ylim(18,37)+
  xlim(24.5,29.6)+
  labs(y = "Salinity (psu)",
       title = " ",
       x = expression(paste("Temperature",~degree,"C")),
       color = expression(paste("Temperature",~degree,"C")))+
  theme_bw()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)
  )

# La Source
P_Source<-AllData %>%
  filter(CowTagID == 5)%>%
  ggplot(aes(x = Depth_m+.1, y = Salinity_psu))+
  geom_point(aes(color = TempInSitu), alpha = 0.5)+
  coord_trans(x = "log")+
  scale_colour_viridis_c(limits = c(24.5,29.6), trans = "log")+
  ylim(18,37)+
  labs(x = "Water Height (m)",
       y = "Salinity (psu)",
       title = "La Source",
       color = expression(paste("Temperature",~degree,"C")))+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14)
          )

P_Source_Temp<-AllData %>%
  filter(CowTagID == 5)%>%
  ggplot(aes(x = TempInSitu, y = Salinity_psu))+
  geom_point(aes(color = TempInSitu), alpha = 0.5)+
 # coord_trans(x = "log")+
  scale_colour_viridis_c(limits = c(24.5,29.6), trans = "log")+
  geom_smooth(method = lm, 
              formula = y ~ splines::ns(log(x), 2), 
              color = "black")+
  ylim(18,37)+
  xlim(24.5,29.6)+
  labs(x = expression(paste("Temperature",~degree,"C")),
       y = "Salinity (psu)",
       title = " ",
       color = expression(paste("Temperature",~degree,"C")))+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5)
  )

(P_lagoon+P_Source)/(P_Lagoon_Temp+P_Source_Temp)+
  plot_layout(guides = "collect")&theme(legend.position = "none", title = element_text(size = 16), plot.title = element_text(hjust = 0.5))
  
ggsave(here("Output","Seeps.png"), width = 8, height = 8)


##### Extract the in situ Temps for pH ###
pHData<-read_csv(here("Data","pHProbe_Data.csv")) %>%
  mutate(Date = mdy_hms(paste(Date,as.character(SamplingTime))))

TempOnly<-AllData %>%
  mutate(CowTagID = as.numeric(as.character(CowTagID)))%>%
  select(Site,CowTagID,Date,TempInSitu) %>%
  right_join(pHData, by = c("Date","Site","CowTagID"))%>%
  arrange(Date,Site,CowTagID)

write_csv(TempOnly,here("Data","pHProbe_temp.csv"))       
