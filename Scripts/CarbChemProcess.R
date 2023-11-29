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

# pHData<-read_csv(here("Data","pHProbe_temp_edit.csv"))%>%
#   mutate(TrisCalDate = ymd(TrisCalDate))

pHData<-read_csv(here("Data","pHProbe_Data.csv"))%>%
  mutate(TrisCalDate = ymd(TrisCalDate))

# Bring in metadata
meta<-read_csv(here("Data","CowTagMetaData.csv"))

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
  mutate(TA = as.numeric(TA),
    pH_insitu = pHinsi(pH = pH, ALK = TA, Tinsi = TempInSitu, Tlab = TempInLab, 
                            S = Salinity_In_Field,Pt = 0, k1k2 = "m10", kf = "dg")) %>%
  # mutate(pH_insitu = pHinsi(pH = pH, ALK = TA, Tinsi = TempInSitu, Tlab = TempInLab, 
  #                           S = Salinity,Pt = Phosphate_umolL, k1k2 = "m10", kf = "dg")) %>%
  select(!pH) %>% # I only need the in situ pH calculation so remove this
  rename(pH = pH_insitu) %>% # rename it 
  ungroup() 

pHSlope$TA[NoTA]<-NA # make TA na again for the missing values

#### Calculate the CO2 data #####
pHSlope_filtered <- pHSlope %>%
  drop_na(TA)

AllCO2<-carb(8, pHSlope_filtered$pH, pHSlope_filtered$TA/1000000, S=pHSlope_filtered$Salinity_In_Field, T=pHSlope_filtered$TempInSitu, Patm=1, P=0, Pt=0, Sit=0,
             k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential", 
             warn="y", eos="eos80")

AllCO2 <- AllCO2 %>%
  mutate(ALK = ALK*1000000,
         CO2 = CO2*1000000,
         CO3 = CO3*1000000,
         DIC = DIC*1000000,
         HCO3 = HCO3*1000000) %>% # convert everything back to umol %>%
  select(DIC, pCO2 = pCO2insitu, CO2, CO3, HCO3, OmegaCalcite, OmegaAragonite) 

AllCO2 <- pHSlope_filtered %>%
  select(Site, CowTagID, SeepCode, Date, SamplingTime)%>%
  bind_cols(AllCO2) %>%
  right_join(pHSlope) %>%
  select(Site, CowTagID, SeepCode, Date, Day_Night, SamplingTime, Tide, Seep_Reef,DIC, pCO2,OmegaAragonite,TA, pH, Salinity_In_Field, TempInSitu)%>%
  mutate(Day_Night = factor(Day_Night, levels = c("Dawn","Noon","Dusk")))

### Make some plots
AllCO2 %>%
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


AllCO2 %>%
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

## plot TA vs DIC

AllCO2 %>%
  filter(Site!= "Lagoon",
         #DIC < 2200
              )%>%
  ggplot(aes(x = DIC*Salinity_In_Field/37, y = TA*Salinity_In_Field/37, color = Site))+
  geom_point()+
  geom_smooth(method = "lm")+
#  geom_label(aes(label = CowTagID))+
  facet_wrap(~Seep_Reef, scales = "free")

AllCO2 %>%
   filter(Site == "Lagoon")%>%
  ggplot(aes(x = DIC*Salinity_In_Field/37, y = TA*Salinity_In_Field/37, color = Site))+
  geom_point()+
  geom_smooth(method = "lm")+
  #  geom_label(aes(label = CowTagID))+
  facet_wrap(~Seep_Reef, scales = "free")

AllCO2 %>%
  filter(CowTagID == 5)%>%
  ggplot(aes(x = Salinity_In_Field, y = pH))+
  geom_point()


AllCO2 %>%
  ggplot(aes(x = Site, color = Day_Night, y = TA))+
  geom_boxplot()+
  facet_wrap(~Seep_Reef, scales = "free")

AllCO2 %>%
  ggplot(aes(x = Site, color = Day_Night, y = TempInSitu))+
  geom_boxplot()+
  facet_wrap(~Seep_Reef)

AllCO2 %>%
  left_join(meta)%>%
  filter(Site != "Lagoon",
         !CowTagID %in%c(1,12,5))%>%
  ggplot(aes(x = log(Depth_logger), color = Site, y = TempInSitu))+
  geom_point()+
  geom_smooth(method = "lm")+
#  geom_label(aes(label = CowTagID))+
  facet_wrap(~Day_Night)

### check temps for 1 and 12... they are too high and probably wrong... look at calibration

## Plot TA vs Salinity
AllCO2 %>%
  filter(Site != "Lagoon")%>%
  droplevels()%>%
  ggplot(aes(color = Site, x = Day_Night , y = pH))+
  geom_boxplot()+
  #geom_label(aes(label = CowTagID))+
  facet_wrap(~Seep_Reef)

