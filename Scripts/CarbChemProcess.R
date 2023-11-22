## process pH for Orion
# Created by Nyssa Silbiger
# Edited on 11/06/2023

library(tidyverse)
library(seacarb)
library(broom)
library(here)
library(lubridate)
library(calecopal)
library(ggridges)


## bring in pH calibration files and raw data files
pHcalib<-read_csv(here("Data","Tris_Calibration_Log.csv")) %>%
  mutate(TrisCalDate = ymd(TrisCalDate))

pHData<-read_csv(here("Data","pHProbe_temp_edit.csv"))%>%
  mutate(TrisCalDate = ymd(TrisCalDate))

## take the mV calibration files by each date and use them to calculate pH
pHSlope<-pHcalib %>%
  filter(HOBO_Orion =="Orion") %>% # extract only the orion data
  nest_by(TrisCalDate)%>%
  mutate(fitpH = list(lm(mVTris~TTris, data = data))) %>% # linear regression of mV and temp of the tris
  reframe(broom::tidy(fitpH)) %>% # make the output tidy
  select(TrisCalDate, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) # put slope and intercept in their own column

pHSlope<-pHData %>%
  right_join(pHSlope, by = "TrisCalDate") %>% # join with the pH sample data
  mutate(mVTris = TempInLab*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  drop_na(TempInSitu)%>%
  drop_na(mV) %>%
  mutate(pH = pH(Ex=mV,Etris=mVTris,S=Salinity_In_Field,T=TempInLab))  # calculate pH of the samples using the pH seacarb function


# The TA data is missing from some bottles because of lack of water, but it really doesnt affect the pH calc.
# I am replacing the missing values with 2300 for the pH calculation then converting it back to NA

NoTA<-which(is.na(pHSlope$TA))
pHSlope$TA[NoTA]<-2300

#Now calculate pH
pHSlope <-pHSlope%>%
  mutate(pH_insitu = pHinsi(pH = pH, ALK = TA, Tinsi = TempInSitu, Tlab = TempInLab, 
                            S = Salinity_In_Field,Pt = 0, k1k2 = "m10", kf = "dg")) %>%
  # mutate(pH_insitu = pHinsi(pH = pH, ALK = TA, Tinsi = TempInSitu, Tlab = TempInLab, 
  #                           S = Salinity,Pt = Phosphate_umolL, k1k2 = "m10", kf = "dg")) %>%
  select(!pH) %>% # I only need the in situ pH calculation so remove this
  rename(pH = pH_insitu) %>% # rename it 
  ungroup() 

pHSlope$TA[NoTA]<-NA # make TA na again for the missing values

### Make some plots
pHSlope %>%
  filter(Seep_Reef == "Seep") %>%
  ggplot(aes(x = Salinity_In_Field, y = pH, color = Site))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "Salinity (psu)",
       y = "pH (total)",
       color = "")+
  scale_color_manual(values = c("#122A64","#01c3e6"))+
  theme_bw()

ggsave(here("Output","pH_Seep.png"), width = 4, height = 3)


pHSlope %>%
  filter(Seep_Reef == "Reef",
         Site == "Lagoon") %>%
  ggplot(aes(x = Salinity_In_Field-0.8, y = pH, color = Tide))+
  geom_point(aes(shape = Day_Night))+
  geom_smooth(method = "lm")+
  labs(x = "Salinity (psu)",
       y = "pH (total)",
       color = "Tide",
       shape = "")+
 # facet_wrap(~Site, scales = "free")+
  theme_bw()

ggsave(here("Output","pH_Lagoon_Reef.png"), width = 4, height = 3)
