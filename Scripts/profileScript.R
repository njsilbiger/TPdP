### Profile Script###
### By Nyssa Silbiger ###
### Created on 2023-11-04 ####


### Load libraries ######
library(here)
library(tidyverse)
library(lubridate)
library(patchwork)

### Read in Data #####
# The Start and stop times
ProfileTimes<-read_csv(here("Data","ProfileTimes.csv"))%>%
  mutate(ProfileStart = mdy_hm(ProfileStart),
         ProfileStop = mdy_hm(ProfileStop)
  )

# The CT data
CT<-read_csv(here("Data","LoggerData","Processed","Calibrated_ProfileBothSites.csv")) %>%
  select(Date, TempInSitu, Salinity_psu)

# The DepthData
WL<-read_csv(here("Data","LoggerData","Processed","Calibrated_ProfilerBothSites_WL.csv"))%>%
  mutate(Date = mdy_hms(Date))%>%
  select(Date, Depth_m)

# Process data ###########
# Join the CT and Depth dat 
AllProfiles<-left_join(CT,WL) %>%
  mutate(cast = NA)

 # add cast number by all start and stop times   
 for(j in 1:nrow(ProfileTimes)){
    AllProfiles$cast<-ifelse(AllProfiles$Date >= ProfileTimes$ProfileStart[j] & AllProfiles$Date <= ProfileTimes$ProfileStop[j] , ProfileTimes$Cast_Num[j],AllProfiles$cast)
  }
  
# add in the site names and drop the extra data
AllProfiles<- AllProfiles %>%
  mutate(Site = case_when(cast %in% c(1:6)~ "La Source",
                          cast %in% c(7:12)~ "Nordhoff")) %>%
  drop_na(cast) %>%
  filter(Salinity_psu>20)# drop the data when sensor was out of the water

# plot the raw data.
AllProfiles %>%
  ggplot(aes(x = Salinity_psu, y = -Depth_m, color = factor(cast)))+
  geom_point()+
 # geom_line()
  facet_wrap(~Site)

# Average by depth bin by rounding to the nearest 10cm
AvgProfile<-AllProfiles %>%
  mutate(Depth_round = round(Depth_m,1))%>%
  group_by(Site, Depth_round)%>%
  summarise(Temp = mean(TempInSitu, na.rm = TRUE),
            Sal = mean(Salinity_psu, na.rm = TRUE))%>%
  ungroup()
# %>%
#   mutate(cast = factor(cast))

# plot the average profiles
#Salinity
SalinityProf<-AvgProfile %>%
  ggplot(aes(x = Sal, y = -Depth_round, color = Site))+
  geom_point(alpha = 0.2)+
  geom_smooth(orientation = "y", se = FALSE, method = lm, formula = y ~ splines::bs(x, 2))+
  scale_color_manual(values = c("#122A64","#01c3e6"))+
  #facet_wrap(~Site)+
  labs(y = "Depth (m)",
       x = "Salinity (psu)",
       color = "")+
  theme_bw()+
  theme(
    legend.position = c(0.45, 0.95),
    legend.justification = c("right", "top"),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )
   
TempProf<-AvgProfile %>%
  ggplot(aes(x = Temp, y = -Depth_round, color = Site))+
  geom_point(alpha = 0.2)+
  geom_smooth(orientation = "y", se = FALSE, method = lm, formula = y ~ splines::bs(x, 2))+
  scale_color_manual(values = c("#122A64","#01c3e6"))+
  #facet_wrap(~Site)+
  labs(y = "",
       x = expression(paste("Temperature",~degree,"C")),
       color = "")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_blank() 
  )

SalinityProf+TempProf
ggsave(here("Output","Profile.png"), width = 8, height = 4)
