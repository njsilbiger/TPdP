## process pH for Orion
# Created by Nyssa Silbiger
# Edited on 4/06/2024

library(tidyverse)
library(seacarb)
library(broom)
library(here)
library(lubridate)
library(calecopal)
library(ggridges)
library(patchwork)
library(lme4)
library(lmerTest)
library(interactions)


## bring in pH calibration files and raw data files
pHcalib<-read_csv(here("Data","Tris_Calibration_Log.csv")) %>%
  mutate(TrisCalDate = ymd(TrisCalDate))

# pHData<-read_csv(here("Data","pHProbe_temp_edit.csv"))%>%
#   mutate(TrisCalDate = ymd(TrisCalDate))

pHData<-read_csv(here("Data","pHProbe_Data.csv"))%>%
  mutate(TrisCalDate = ymd(TrisCalDate))

# Bring in metadata
meta<-read_csv(here("Data","CowTagMetaData.csv"))

# fdom data
fdom<-read_csv(here("Data","fDOM_Tahiti.csv"))

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
  mutate(pH = pH(Ex=mV,Etris=mVTris,S=Salinity_In_Lab,T=TempInLab))  # calculate pH of the samples using the pH seacarb function


# The TA data is missing from some bottles because of lack of water, but it really doesnt affect the pH calc.
# I am replacing the missing values with 2300 for the pH calculation then converting it back to NA


#Now calculate pH
pHSlope <-pHSlope%>%
  mutate(TA = as.numeric(TA),
    pH_insitu = pHinsi(pH = pH, ALK = TA, Tinsi = TempInSitu, Tlab = TempInLab, 
                            S = Salinity_In_Lab,Pt = 0, k1k2 = "m06", kf = "dg", ks = "d")) %>%
  # mutate(pH_insitu = pHinsi(pH = pH, ALK = TA, Tinsi = TempInSitu, Tlab = TempInLab, 
  #                           S = Salinity,Pt = Phosphate_umolL, k1k2 = "m10", kf = "dg")) %>%
  select(!pH) %>% # I only need the in situ pH calculation so remove this
  rename(pH = pH_insitu) %>% # rename it 
  ungroup() 


#### Calculate the CO2 data #####

AllCO2<-carb(8, pHSlope$pH, pHSlope$TA/1000000, S=pHSlope$Salinity_In_Lab, T=pHSlope$TempInSitu, Patm=1, P=0, Pt=0, Sit=0,
             k1k2="m06", kf="dg", ks="d", pHscale="T", b="u74", gas="potential", 
             warn="y", eos="eos80")

AllCO2 <- AllCO2 %>%
  mutate(ALK = ALK*1000000,
         CO2 = CO2*1000000,
         CO3 = CO3*1000000,
         DIC = DIC*1000000,
         HCO3 = HCO3*1000000) %>% # convert everything back to umol %>%
  select(DIC, pCO2 = pCO2insitu, CO2, CO3, HCO3, OmegaCalcite, OmegaAragonite) 

AllCO2 <- pHSlope %>%
  select(Site, CowTagID, SeepCode, Date, SamplingTime)%>%
  bind_cols(AllCO2) %>%
  right_join(pHSlope) %>%
  select(Site, CowTagID, SeepCode, Date, Day_Night, SamplingTime, Tide, Seep_Reef,DIC, pCO2,OmegaAragonite,TA, pH, Salinity_In_Lab, TempInSitu,NN, NO2, PO4, SiO2, NH3 )%>%
  mutate(Day_Night = factor(Day_Night, levels = c("Dawn","Noon","Dusk")))

### Make some plots
AllCO2 %>%
  filter(Seep_Reef == "Seep") %>%
  ggplot(aes(x = Salinity_In_Lab, y = pH, color = Site))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "Salinity (psu)",
       y = "pH (total)",
       color = "")+
  scale_color_manual(values = c("#122A64","#01c3e6"))+
  theme_bw()

ggsave(here("Output","pH_Seep.png"), width = 4, height = 3)


AllCO2 %>%
  mutate(Day_Night = factor(Day_Night, levels = c("Dawn", "Noon","Dusk")))%>%
  filter(Seep_Reef == "Reef",
         Site == "Lagoon"
         ) %>%
  ggplot(aes(x = Salinity_In_Lab, y = pH, color = Day_Night))+
  geom_point(aes())+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c("#fbc540","#ec5c04","#ad3304"))+
  labs(x = "Salinity (psu)",
       y = "pH (total)",
       color = "",
       shape = "")+
 # facet_wrap(~Site, scales = "free")+
  theme_bw()

ggsave(here("Output","pH_Lagoon_Reef.png"), width = 4, height = 3)

## plot TA vs DIC

TADIC<-AllCO2 %>%
  filter(#Site!= "Lagoon",
         Seep_Reef == "Reef"
         #DIC < 2200
              )%>%
  ggplot(aes(x = DIC*Salinity_In_Lab/36, y = TA*Salinity_In_Lab/36, color = Site))+
  geom_point()+
  geom_smooth(method = "lm")+
 # geom_label(aes(label = CowTagID))+
  scale_color_manual(values = c("#122A64","#01c3e6","#5F9EA0"))+
  labs(x = expression(paste("Salinity-normalized DIC (",mu,"mol kg"^-1,")")),
       y = expression(paste("Salinity-normalized TA (",mu,"mol kg"^-1,")")),
  )+
  facet_wrap(~Day_Night, ncol = 1)+
  theme_bw()+
  theme(legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.background =element_blank(),
    strip.text = element_text(size = 16, face = "bold"))
      

AllCO2 <-AllCO2 %>%
  mutate(TA_salnorm = TA*Salinity_In_Lab/36, # salinity normalized
         DIC_salnorm = DIC*Salinity_In_Lab/36)

#### bring in the fdom data
AllCO2<-AllCO2 %>% 
  left_join(fdom) 

modTADIC<-lm(TA_salnorm~DIC_salnorm*Day_Night*Site , data =AllCO2 %>%
                 filter(#Site!= "Lagoon",
                   Seep_Reef == "Reef") )

modTADIC<-lmer(TA_salnorm~DIC_salnorm*Day_Night*Site +(1|CowTagID), data =AllCO2 %>%
                 filter(#Site!= "Lagoon",
                   Seep_Reef == "Reef") )

anova(modTADIC)
summary(modTADIC)

ss <- sim_slopes(modTADIC, pred = DIC_salnorm, modx = Site, mod2 = Day_Night, johnson_neyman = FALSE)
plot(ss)

# extract the individual slopes
#### need to work on this
test<-bind_rows(data.frame(ss$slopes[[2]]),data.frame(ss$slopes[[3]])) %>%
  mutate_at(.vars = c("Est.","S.E.","X2.5.","X97.5.", "t.val.","p"), as.numeric)%>%
  mutate(Day_Night = c("Noon","Noon","Noon", "Dusk","Dusk","Dusk"))

a<-data.frame(ss$slopes[[1]]) %>%
  mutate_at(.vars = c("Est.","S.E.","X2.5.","X97.5.", "t.val.","p"), as.numeric)%>%
  mutate(Day_Night = c("Dawn","Dawn","Dawn"))

slopes<-bind_rows(a,test) %>%
  mutate(Day_Night = factor(Day_Night, levels = c("Dusk","Noon","Dawn")))%>%
  rename(Site = Value.of.Site)

testint<-bind_rows(data.frame(ss$ints[[2]]),data.frame(ss$ints[[3]])) %>%
  mutate_at(.vars = c("Est.","S.E.","X2.5.","X97.5.", "t.val.","p"), as.numeric)%>%
  mutate(Day_Night = c("Noon","Noon","Noon", "Dusk","Dusk","Dusk"))

aint<-data.frame(ss$ints[[1]]) %>%
  mutate_at(.vars = c("Est.","S.E.","X2.5.","X97.5.", "t.val.","p"), as.numeric)%>%
  mutate(Day_Night = c("Dawn","Dawn","Dawn"))

ints <-bind_rows(aint,testint) %>%
  mutate(Day_Night = factor(Day_Night, levels = c("Dusk","Noon","Dawn")))%>%
  rename(Site = Value.of.Site)

P_estimate<-slopes %>%
  ggplot(aes(color = Site, x = Est., y = Day_Night))+
  geom_point(size = 4, position = position_dodge(0.2))+
  geom_errorbarh(aes(xmin = Est. - S.E.,xmax = Est. + S.E.  ), height = 0.01, position = position_dodge(0.2))+
  scale_color_manual(values = c("#122A64","#01c3e6","#5F9EA0"))+
  labs(x = "TA/DIC Slopes",
       y = "",
       color = "")+
  theme_bw()+
  theme(legend.position = c(0.25, 0.9),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.y = element_blank()
        )

P_ints<-ints %>%
  ggplot(aes(color = Site, x = Est., y = Day_Night))+
  geom_point(size = 4)+
  geom_errorbarh(aes(xmin = Est. - S.E.,xmax = Est. + S.E.  ), height = 0.01)+
  scale_color_manual(values = c("#122A64","#01c3e6","#5F9EA0"))+
  labs(x = "Intercept",
       y = "",
       color = "")+
  theme_bw()+
  theme(legend.position = c(0.15, 0.1),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.y = element_blank()
  )


TADIC|P_estimate

ggsave(here("Output","TADIC.png"), width = 9, height = 8)

AllCO2 %>%
  filter(Site== "Lagoon",
         Seep_Reef == "Reef"
         #DIC < 2200
  )%>%
  ggplot(aes(x = DIC*Salinity_In_Lab/36, y = TA*Salinity_In_Lab/36, color = Day_Night))+
  geom_point()+
  geom_smooth(method = "lm")+
 # geom_label(aes(label = CowTagID))
  facet_wrap(~Day_Night, scales = "free")

AllCO2 %>%
  filter(#Site!= "Lagoon",
         Seep_Reef == "Reef"
         #DIC < 2200
  )%>%
  group_by(Site, Day_Night)%>%
  summarise(salmean = mean(Salinity_In_Lab, na.rm = TRUE),
            salse = sd(Salinity_In_Lab, na.rm = TRUE)/sqrt(n()))%>%
  ggplot(aes(x = Site, y = salmean, color = Day_Night))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = salmean - salse, ymax = salmean+salse), width = 0.1)

sal1<-AllCO2 %>%
  filter(#Site!= "Lagoon",
    Seep_Reef == "Reef"
    #DIC < 2200
  )%>%
  group_by(Site,Seep_Reef, Day_Night)%>%
  summarise(salmean = mean(Salinity_In_Lab, na.rm = TRUE),
            salse = sd(Salinity_In_Lab, na.rm = TRUE)/sqrt(n()))%>%
  ggplot(aes(x = Site, y = salmean, color = Day_Night))+
  geom_point(size = 3, position = position_dodge(.2))+
  geom_errorbar(aes(ymin = salmean - salse, ymax = salmean+salse), width = 0.1, position = position_dodge(.2))+
  scale_color_manual(values = c("#fbc540","#ec5c04","#ad3304"))+
  labs(y = "Mean salinity (psu)",
         x = "",
         color = "")+
 # facet_wrap(~Seep_Reef, scale = "free")+
  theme_bw()+
  theme(legend.position = "none",
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
       # axis.text.y = element_blank()
  )

sal2<-AllCO2 %>%
  filter(#Site== "Lagoon",
         Seep_Reef == "Seep"
         #DIC < 2200
  )%>%
  group_by(Site,Seep_Reef, Day_Night)%>%
  summarise(salmean = mean(Salinity_In_Lab, na.rm = TRUE),
            salse = sd(Salinity_In_Lab, na.rm = TRUE)/sqrt(n()))%>%
  ggplot(aes(x = Site, y = salmean, color = Day_Night))+
  geom_point(size = 3, position = position_dodge(.2))+
  geom_errorbar(aes(ymin = salmean - salse, ymax = salmean+salse), width = 0.1, position = position_dodge(.2))+
  scale_color_manual(values = c("#fbc540","#ec5c04","#ad3304"))+
  labs(y = "Mean salinity (psu)",
       x = "",
       color = "")+
  facet_wrap(~Site, scale = "free")+
  theme_bw()+
  theme(#legend.position = c(0.15, 0.1),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
  #  axis.text.y = element_blank(),
  strip.text = element_blank()
  )


sal1|sal2
ggsave(here("Output","SalMean.png"), width = 10, height = 4)

AllCO2 %>%
   filter(Site == "Lagoon",
          Seep_Reef == "Reef")%>%
  ggplot(aes(x = DIC*Salinity_In_Lab/36, y = TA*Salinity_In_Lab/36, color = Day_Night))+
  geom_point()+
  geom_smooth(method = "lm")
#  geom_label(aes(label = CowTagID))+
 # facet_wrap(~Day_Night)

  
AllCO2 %>%
  #filter(CowTagID == 5)%>%
  # something off with 10 noon
  filter(!CowTagID %in% c(5,41))%>%
  ggplot(aes(x = DIC, y = pH))+
  geom_point()+
#  geom_label(aes(label = CowTagID))+
  geom_smooth(method = "lm")+
  facet_wrap(~Site)

AllCO2 %>%
  filter(Seep_Reef == "Seep")%>%
  ggplot(aes(x = Salinity_In_Lab, y = TA, color = Site))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c("#122A64","#01c3e6"))+
  facet_wrap(~Site, scale = "free", ncol = 1)+
  labs(y = expression(paste("TA (",mu,"mol kg"^-1,")")),
       x = "Salinity (psu)")+
  theme_bw()+
  theme(legend.position = "none",
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        strip.background = element_blank())

ggsave(here("Output","TA_sal.png"), width = 4, height = 8)

AllCO2 %>%
  filter(Seep_Reef == "Seep")%>%
  ggplot(aes(x = Salinity_In_Lab, y = TempInSitu))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Site, scale = "free")+
  theme_bw()

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
#  filter(Site == "Lagoon")%>%
  droplevels()%>%
  ggplot(aes(color = Site, x = Day_Night , y = TA))+
  geom_boxplot()+
  #geom_label(aes(label = CowTagID))+
  facet_wrap(~Seep_Reef, scale = "free")


AllCO2 %>%
  left_join(meta)%>%
  filter(Site != "Lagoon",
         !CowTagID %in%c(1,12,5))%>%
  ggplot(aes(x = Salinity_In_Lab, color = Site, y = TempInSitu))+
  geom_point()+
  geom_smooth(method = "lm")+
  #  geom_label(aes(label = CowTagID))+
  facet_wrap(~Day_Night, scale = "free")

#### Turbinaria data #####
NN<-meta %>%
  ggplot(aes(x = del15N, y = N_percent))+
  geom_point()+
  geom_smooth(data =  meta %>% filter(Site != "Nordhoff"),method = "lm")+
  geom_point(data = meta %>% filter(CowTagID %in% c(41,5)), color = "red")+
  facet_wrap(~Site, scales = "free")+
  theme_bw()

ggsave(here("Output","NvsN.png"), width = 8, height = 4)

# mean and SE for the different values

TurbMeans<-meta %>%
  filter(CowTagID != 41)%>% # remove the lagoon seep
  select(Site, del15N:C_percent)%>%
  pivot_longer(cols = del15N:C_percent)%>%
  group_by(Site, name)%>%
  summarise(MeanValue = mean(value,na.rm = TRUE),
            MeanSE = sd(value, na.rm = TRUE)/sqrt(n()))

## pull out the seep samples only
SeepN<-meta %>% 
  filter(CowTagID %in% c(5,41))%>%
  select(Site, del15N:C_percent)%>%
  pivot_longer(cols = del15N:C_percent) %>%
  rename(MeanValue = value)
  

TMeans<-TurbMeans %>%
  mutate(Site = factor(Site, levels = c("Lagoon","La Source","Nordhoff")))%>%
  filter(name %in% c("del15N","N_percent"))%>%
  ggplot(aes(x = Site, y = MeanValue))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = MeanValue-MeanSE, ymax = MeanValue+MeanSE), width = 0.1)+
#  geom_point(data = SeepN, aes(x = Site, y = MeanValue), color = "red")+
  facet_wrap(~name, scales = "free")+
  theme_bw()

#TMeans|NN
ggsave(here("Output","TurbData.png"), width = 8, height = 4)

SeepN %>%
  mutate(Site = factor(Site, levels = c("Lagoon","La Source","Nordhoff")))%>%
  filter(name %in% c("del15N","N_percent"))%>%
  ggplot(aes(x = Site, y = MeanValue))+
  geom_point(size = 3)+
  #  geom_point(data = SeepN, aes(x = Site, y = MeanValue), color = "red")+
  facet_wrap(~name, scales = "free")+
  theme_bw()

#TMeans|NN
ggsave(here("Output","TurbData_Seep.png"), width = 8, height = 4)

#### water column nutrient plots

AllCO2 %>%
  ggplot(aes(x = SiO2, y = NN, color = Site))+
  geom_point()

AllCO2 %>%
#  filter(Seep_Reef == "Seep")%>%
#  filter(Site != "Lagoon")%>%
  ggplot(aes(x = SiO2, y = NN, color = Site))+
  geom_point()+
  facet_wrap(Seep_Reef~Site)


AllCO2 %>%
  filter(Seep_Reef == "Seep")%>%
#  filter(Site != "Lagoon")%>%
  ggplot(aes(x = SiO2, y = Salinity_In_Lab, color = Site))+
  geom_point()

AllCO2 %>%
 # filter(Site != "Lagoon")%>%
  drop_na(NN,SiO2)%>%
  ggplot(aes(x = SiO2, y = NN, color = Site))+
  geom_point()+
 # coord_trans(x = "log", y = "log")+
  geom_smooth(method = "lm")+
  facet_wrap(~Site)

AllCO2 %>%
  # filter(Site != "Lagoon")%>%
  filter(Seep_Reef == "Seep")%>%
  drop_na(NN,SiO2)%>%
  ggplot(aes(x = pH, y = NN, color = Site))+
  geom_point()+
  # coord_trans(x = "log", y = "log")+
  geom_smooth(method = "lm")

AllCO2 %>%
# filter(!CowTagID %in% c(13,16))%>%
  filter(Seep_Reef == "Reef")%>%
  ggplot(aes(x = log(SiO2), y = NN, color = Site))+
  geom_point()+
#  geom_label(aes(label = CowTagID))+
#  coord_trans(x = "log", y = "log")+
  #geom_smooth(method = "lm")
  facet_wrap(~Day_Night)


AllCO2 %>%
  left_join(meta) %>%
  filter(Seep_Reef == "Reef")%>%
 # filter(!CowTagID %in% c(13,16))%>%
  ggplot(aes(x = Depth_logger, y = NN, color = Site))+
  geom_point()+
  facet_wrap(Site~Day_Night, scale = "free")


AllCO2 %>%
  left_join(meta) %>%
  group_by(Site, CowTagID)%>%
  summarise_at(vars(DIC:DIC_salnorm, del15N:C_percent), .funs = mean)%>%
  drop_na()%>%
  ggplot(aes(x = NN, y = N_percent, color = Site))+
  geom_point()+
#  geom_label(aes(label = CowTagID))+
  facet_wrap(~Site, scale = "free")


AllCO2 %>%
#  left_join(meta) %>%
  drop_na(NN)%>%
  filter(Seep_Reef == "Reef")%>%
  ggplot(aes(x = Site, y = NN, color = Day_Night))+
  geom_boxplot()+
  geom_jitter(position = position_dodge(width = 0.8))
#  geom_label(aes(label = CowTagID),position = position_dodge(width = 0.8))


AllCO2 %>%
  #  left_join(meta) %>%
  drop_na(NN)%>%
  ggplot(aes(x = PO4, y = NN))+
  geom_point(aes(color = Day_Night))+
  facet_wrap(~CowTagID, scale = "free")

AllCO2 %>%
 # filter(!CowTagID %in% c(13,16))%>%
  #filter(Seep_Reef == "Reef")%>%
  group_by(Site, Seep_Reef)%>%
  summarise(NN_mean = mean(NN, na.rm = TRUE),
            NN_se = sd(NN,na.rm = TRUE)/sqrt(n()))%>%
  mutate(Site = factor(Site, levels = c("Lagoon", "La Source","Nordhoff")))%>%
  ggplot(aes(x = Site))+
  geom_point(aes(y = NN_mean), size = 3)+
  geom_errorbar(aes(ymin = NN_mean - NN_se, ymax = NN_mean+NN_se), width = .1)+
  labs(y = expression(paste("Mean nitrate + nitrite (",mu,"mol L"^-1,")")))+
  facet_wrap(~Seep_Reef, scales = "free")+
  theme_bw()
ggsave(here("Output","NN_mean.png"), width = 8, height = 4)

modN<-lmer(NN~Site +(1|Day_Night), data = AllCO2 %>%
             filter(!CowTagID %in% c(13,16)) %>%
             filter(Seep_Reef == "Reef"))
anova(modN)
modN

## plot fdom

AllCO2 %>%
#  filter(Seep_Reef == "Reef")%>%
  ggplot(aes(x = Salinity_In_Lab, y = MarineHumic+UVHumic+VisibleHumic, color = Site))+
  geom_point()+
  #  geom_label(aes(label = CowTagID))+
  #  coord_trans(x = "log", y = "log")+
  geom_smooth(method = "lm")+
  labs(x = "Salinity (psu)", y = "Humic-like fDom (raman units)")+
  facet_wrap(Seep_Reef~Site, scales = "free")

ggsave(here("Output","humics.png"), width = 8, height = 4)


AllCO2 %>%
  #  filter(Seep_Reef == "Reef")%>%
  ggplot(aes(x = Salinity_In_Lab, y = Tyrosine+Tryptophan+Phenylalanine))+
  geom_point()+
  #  geom_label(aes(label = CowTagID))+
  #  coord_trans(x = "log", y = "log")+
  geom_smooth(method = "lm")+
  labs(x = "Salinity (psu)", y = "Proteinaceous-like fDom (raman units)")+
  facet_wrap(Seep_Reef~Site, scales = "free")

ggsave(here("Output","prot.png"), width = 8, height = 4)

AllCO2 %>%
  #  filter(Seep_Reef == "Reef")%>%
  ggplot(aes(x = Salinity_In_Lab, y = MC))+
  geom_point()+
  #  geom_label(aes(label = CowTagID))+
  #  coord_trans(x = "log", y = "log")+
  geom_smooth(method = "lm")+
  labs(x = "Salinity (psu)", y = "M to C")+
  facet_wrap(Seep_Reef~Site, scales = "free")
ggsave(here("Output","MC.png"), width = 8, height = 4)

### pca for fdom ###

## make a pca of all the fdom data for the reef
pcadata<-AllCO2  %>%
  filter(Seep_Reef == "Reef") %>%
  ungroup()%>%
  dplyr::select(UVHumic:MC) %>%
  drop_na()

# run a pca
pca<-prcomp(pcadata, scale. = TRUE, center = TRUE)

# calculate percent explained by each PC
perc.explained<-round(100*pca$sdev/sum(pca$sdev),1)

# Extract the scores and loadings
PC_scores <-as_tibble(pca$x[,1:2])

PC_loadings<-as_tibble(pca$rotation)%>%
  bind_cols(labels = rownames(pca$rotation))


data_pca<-AllCO2  %>%
  filter(Seep_Reef == "Reef") %>%
  drop_na(UVHumic:MC)%>%
  bind_cols(PC_scores)%>%
  mutate(Site = factor(Site, levels = c("Lagoon","La Source","Nordhoff")))


p1<-data_pca %>%
  #mutate(time_point = factor(time_point))%>%
  ggplot(aes(x = PC1, y = PC2, color = Site, shape = Day_Night))+
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  # scale_shape_manual(values = c(1, 22,15,16))+
  scale_colour_manual(values = c("#01c3e6","#122A64","#5F9EA0"))+
  scale_fill_manual(values = c("#01c3e6","#122A64","#5F9EA0"))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ggforce::geom_mark_ellipse(
    aes(fill = Site, label = paste(Site, Day_Night), color =Site),
    alpha = .35, show.legend = FALSE,  label.buffer = unit(1, "mm"), con.cap=0, tol = 0.05)+
  geom_point(size = 2) +
  labs(
    x = paste0("PC1 ","(",perc.explained[1],"%)"),
    y = paste0("PC2 ","(",perc.explained[2],"%)"))+
  theme_bw()+
  theme(legend.position = "none",
     #   axis.text.x = element_blank(),
    #    axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        strip.background = element_blank(),
        #      strip.text = element_blank()
  )+
  facet_wrap(~Site)
  
  

p2<-PC_loadings %>%
  ggplot(aes(x=PC1, y=PC2, label=labels))+
  geom_richtext(aes(x = PC1*8+0.1, y = PC2*8+.1 ), show.legend = FALSE, size = 5, fill=NA, label.colour = NA) +
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1*8,yend=PC2*8),size = 1.2,
               arrow=arrow(length=unit(0.1,"cm")))+
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  labs(
    y = "",
    x = "")+
  #y = paste0("PC2 ","(",perc.explained_both[2],"%)"))+
  #  scale_color_manual(values = wes_palette("Darjeeling1"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = c(0.75, 0.75),
        legend.position = "none",
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))


p1/plot_spacer()/p2+plot_layout(heights = c(2,-0.1, 1))

ggsave(here("Output","fdom_pca_reef.png"), width = 10, height = 6)

##### Do a unique pca by site to see how the correlates differ.  Do for both seep and reef

### Now the same with the alll biogeochem data 
## make a pca of all the fdom data for the reef
pcadata<-AllCO2  %>%
  filter(Seep_Reef == "Reef") %>%
  ungroup()%>%
  dplyr::select(Salinity_In_Lab,pH, TA, NN, PO4, SiO2, UVHumic:MC) %>%
  drop_na()

# run a pca
pca<-prcomp(pcadata, scale. = TRUE, center = TRUE)

# calculate percent explained by each PC
perc.explained<-round(100*pca$sdev/sum(pca$sdev),1)

# Extract the scores and loadings
PC_scores <-as_tibble(pca$x[,1:2])

PC_loadings<-as_tibble(pca$rotation)%>%
  bind_cols(labels = rownames(pca$rotation))


data_pca<-AllCO2  %>%
  filter(Seep_Reef == "Reef") %>%
  drop_na(Salinity_In_Lab,pH, TA, NN, PO4, SiO2, UVHumic:MC)%>%
  bind_cols(PC_scores)%>%
  mutate(Site = factor(Site, levels = c("Lagoon","La Source","Nordhoff")))


p1<-data_pca %>%
  #mutate(time_point = factor(time_point))%>%
  ggplot(aes(x = PC1, y = PC2, color = Site, shape = Day_Night))+
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  # scale_shape_manual(values = c(1, 22,15,16))+
  scale_colour_manual(values = c("#01c3e6","#122A64","#5F9EA0"))+
  scale_fill_manual(values = c("#01c3e6","#122A64","#5F9EA0"))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ggforce::geom_mark_ellipse(
    aes(fill = Site, label = paste(Site, Day_Night), color =Site),
    alpha = .35, show.legend = FALSE,  label.buffer = unit(1, "mm"), con.cap=0, tol = 0.05)+
  geom_point(size = 2) +
  labs(
    x = paste0("PC1 ","(",perc.explained[1],"%)"),
    y = paste0("PC2 ","(",perc.explained[2],"%)"))+
  theme_bw()+
  theme(legend.position = "none",
        #   axis.text.x = element_blank(),
        #    axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        strip.background = element_blank(),
        #      strip.text = element_blank()
  )+
  facet_wrap(~Site)



p2<-PC_loadings %>%
  ggplot(aes(x=PC1, y=PC2, label=labels))+
  geom_richtext(aes(x = PC1*8+0.1, y = PC2*8+.1 ), show.legend = FALSE, size = 5, fill=NA, label.colour = NA) +
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_segment(data = PC_loadings, aes(x=0,y=0,xend=PC1*8,yend=PC2*8),size = 1.2,
               arrow=arrow(length=unit(0.1,"cm")))+
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  labs(
    y = "",
    x = "")+
  #y = paste0("PC2 ","(",perc.explained_both[2],"%)"))+
  #  scale_color_manual(values = wes_palette("Darjeeling1"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = c(0.75, 0.75),
        legend.position = "none",
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))


p1/plot_spacer()/p2+plot_layout(heights = c(2,-0.1, 1))

ggsave(here("Output","fdom_pca_reef.png"), width = 10, height = 6)
