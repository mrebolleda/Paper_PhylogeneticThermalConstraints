#
#
# Here we obtain climate series data for the 15 sites sampled
# We calculate 
# 1) temperature and precipitation soil bioclimatic variables
# 2) annual temperature distributions per site and climate profile
# 3) plant growth season 
# 4) temperature distributions per site and climate profile during the growth season
# 5) growth season temperature distributions per site and climate profile
#
#
#### Setting the directory and packages ####
setwd("PATH_TO_YOUR_DIRECTORY")

library(data.table)
library(mgcv)
library(tidyverse)
library(lubridate)
library(rTPC)
library(nls.multstart)
library(broom)
library(dplyr)
library(multcomp)
library(ggrepel)
library(ggforce)
library(MuMIn)
library(viridis)
library(viridisLite)
library(tidyr)
library(purrr)
library(data.table)
library(ggplot2)
library(ggfortify)
library(ggtree)
library(extrafont)
library(Rttf2pt1)
font_import("Futura")
loadfonts("win")
fonts()
library(stringr)
library(readxl)


#### Loading all temperature data from dendra database ####

# We downloaded available temperature values from September 2010 through December 2022

dendra_data <- fread("temperature-data-ucnrs-240502.csv")

#### Processing temperature data for 13 reserves ####

# We first extract Year, Month, and Day from the Time column
dendra_data <- dendra_data %>%
  mutate(Year = year(Time),
         Month = month(Time),
         Day = day(Time))

# We then reshape the dataframe to long format and extract Site names (Runtime ~ 1 hr)
dendra_data_long <- dendra_data %>%
  pivot_longer(cols = -c(Time, Year, Month, Day), names_to = "Site_Metric", values_to = "Value") %>%
  separate(Site_Metric, into = c("Site", "Metric"), sep = "(?=Soil)", extra = "merge")%>%
  pivot_wider(names_from = Metric, values_from = Value)

# We remove outliers due to sensor malfunction or fire events that would distort the average values
dendra_data_long <- dendra_data_long %>% filter(!(Site == "Bodega " & Year < 2019 )) # remove years where Bodega thermometer was very erratic
dendra_data_long <- dendra_data_long %>% filter (!`Soil Temp 50 mm Max degC`>= 65)
dendra_data_long <- dendra_data_long %>% filter (!`Soil Temp 50 mm Min degC`<= -23)

# We then average the daily extremes to estimate mean values. We do this to minimize differences between reserves that have different time resolution. 
# While some reserves provide measurements every 10 min, others only provide daily minima and maxima
dendra_data_long$`Soil Temp 50 mm Avg` <- (dendra_data_long$`Soil Temp 50 mm Max degC` + dendra_data_long$`Soil Temp 50 mm Min degC`)/2

#dendra_data_long <- dendra_data_long %>% dplyr::select (!`Soil Temp 50 mm Avg degC`)
dendra_data_long <- dendra_data_long %>% filter(!is.na(`Soil Temp 50 mm Avg`))

write.csv(dendra_data_long, "dendra_data_long.csv",row.names = FALSE)


dendra_data_long <- fread("dendra_data_long.csv")

dendra_data_long_by_day <- dendra_data_long %>%
  group_by(Site, Year, Month, Day) %>%
  summarise(
    `Soil Temp 50 mm Max degC` = max(`Soil Temp 50 mm Max degC`, na.rm = TRUE), # max day temperature
    `Soil Temp 50 mm Min degC` = min(`Soil Temp 50 mm Min degC`, na.rm = TRUE), # mix day temperature
    `Soil Temp 50 mm Avg` = mean(`Soil Temp 50 mm Avg`, na.rm = TRUE), # mean temperature of a day
    `Soil Temp 50 mm Range` =  (`Soil Temp 50 mm Max degC` - `Soil Temp 50 mm Min degC`) # temperature range of a day
  )


#### Summarizing the data by complete day, month, year, and fully ####



# Summarize by year on a monthly scale
dendra_data_summarized_by_month <- dendra_data_long_by_day %>%
  group_by(Site, Year, Month) %>%
  summarise(
    MT = mean(`Soil Temp 50 mm Avg`, na.rm = TRUE), # mean temperature of a given day in a month and year
    MTR = max(`Soil Temp 50 mm Max degC`, na.rm = TRUE) - min(`Soil Temp 50 mm Min degC`, na.rm = TRUE), # temperature range of a given month and year
    MDTR = mean(`Soil Temp 50 mm Range`, na.rm = TRUE), # mean daily temperature range for each day of a given month 
    MTmax = max(`Soil Temp 50 mm Max degC`, na.rm = TRUE), # maximum temperature of a given month and year
    MMTmax = mean(`Soil Temp 50 mm Max degC`, na.rm = TRUE), # mean maximum temperature of a given month and year
    MTmin = min(`Soil Temp 50 mm Min degC`, na.rm = TRUE), # minimum temperature of a given month and year
    MMTmin = mean(`Soil Temp 50 mm Min degC`, na.rm = TRUE),
    MT_sd = sd(`Soil Temp 50 mm Avg`, na.rm = TRUE),
    MTR_sd = sd(`Soil Temp 50 mm Max degC` - `Soil Temp 50 mm Min degC`,na.rm = TRUE), 
    MTmax_sd = sd(`Soil Temp 50 mm Max degC`, na.rm = TRUE), 
    MTmin_sd = sd(`Soil Temp 50 mm Max degC`, na.rm = TRUE)
  )

# We filter out NA values, group by Site and Year, and keep full to not bias calculations towards a specific season
dendra_data_by_month_filtered <- dendra_data_summarized_by_month %>%
  filter(!is.na(MT_sd) & !is.na(MTR_sd) & !is.na(MTmax_sd) & !is.na(MTmin_sd))%>%
  group_by(Site, Year) %>%
  filter(n() == 12) %>%
  ungroup()

# Summarize by year on a yearly scale
dendra_data_summarized_by_year <- dendra_data_by_month_filtered %>%
  group_by(Site, Year) %>%
  summarise(
    AT = mean(MT, na.rm = TRUE), # mean temperature of a given year
    ADTR = mean(MDTR, na.rm = TRUE), # mean daily temperature range for each day of a given year
    AMTmax = max(MTmax, na.rm = TRUE), # maximum temperature of the warmest month of a given year
    AMTmin = min(MTmin, na.rm = TRUE), # minimum temperature of the coldest month of a given year
    ATR = (AMTmax - AMTmin), # temperature range of a given year
    AMTR = mean(MTR, na.rm = TRUE), # minimum temperature of a given month and year
    Seasonality = (sd(MT, na.rm = TRUE)) # SD of the mean annual temperature => seasonality 
  )

dendra_data_summarized_by_year %>%
  ggplot(aes(x = AT, y = ATR, color = Site))+
  geom_point(size = 4)+
  geom_point(data = dendra_data_summarized_by_year,aes(x = mean(AT),y = mean(ATR),group = Site, color = Site),size = 6)+
  theme_classic()


# We obtain bioclimatic variables by summarizing the data by Site across years
dendra_data_summarized_by_site <- dendra_data_summarized_by_year %>%
  group_by(Site) %>%
  summarise(
    Seasonality = mean(Seasonality,na.rm = TRUE), # SD of the mean annual temperature => seasonality 
    MAT = mean(AT, na.rm = TRUE), # mean temperature of a given year (SBIO1)
    MATR = mean(ATR, na.rm = TRUE), # mean temperature annual range (SBIO7) 
    MMTR = mean(AMTR, na.rm = TRUE), # mean temperature range per month of a given year (SBIO2)
    MDTR = mean(ADTR, na.rm = TRUE), # mean temperature range per day of a given year (SBIO2)
    MATmax = max(AMTmax), # maximum temperature of the hottest month (SBIO5)
    MMTmax = mean(AMTmax, na.rm = TRUE), # mean maximum temperature per month
    MATmin = min(AMTmin), # minimum temperature of the coldest month (SBIO6)
    MMTmin = mean(AMTmin, na.rm = TRUE), # mean minimum temperature per month
    Isothermality = (MMTR/MATR)*100 # (SBIO3)
  )


dendra_data_summarized_by_site %>% ggplot(aes(x = MAT, y = MATR,color = Site))+
  geom_point(size = 7)+
  scale_color_viridis_d()+
  theme_classic()

dendra_data_summarized_by_site %>% ggplot(aes(x = MMTR, y = MATR,color = Site))+
  geom_point(size = 7)+
  scale_color_viridis_d()+
  theme_classic()

dendra_data_summarized_by_site %>% ggplot(aes(x = MDTR, y = MATR,color = Site))+
  geom_point(size = 7)+
  scale_color_viridis_d()+
  theme_classic()

dendra_data_summarized_by_site %>% ggplot(aes(x = MDTR, y = MMTR,color = Site))+
  geom_point(size = 7)+
  scale_color_viridis_d()+
  theme_classic()

dendra_data_summarized_by_site %>% ggplot(aes(x = Site, y = MDTR/MMTR,color = Site))+
  geom_point(size = 7)+
  scale_color_viridis_d()+
  theme_classic()

dendra_data_summarized_by_site %>% ggplot(aes(x = Site, y = MMTR/MATR,color = Site))+
  geom_point(size = 7)+
  scale_color_viridis_d()+
  theme_classic()


#### Loading Merced Vernal Pools Nature Reserve (MVPNR) temperature data ####

# The Merced Vernal Pools Nature Reserve does not have data uploaded to the Dendra database as of 2023, 
# so we obtained data from the San Joaquin River weather station (https://water.ca.gov/). Data was available
# only at a depth of 100 mm, which can underestimate the thermal variability experienced on site by the sampled bacteria

merced_data <- fread("UCM_temperature_data.csv")

# We first extract Year, Month, and Day from the Time column
merced_data <- merced_data %>%
  mutate(`DATE TIME` = mdy_hm(`DATE TIME`),
         Year = year(`DATE TIME`),
         Month = month(`DATE TIME`),
         Day = day(`DATE TIME`))

# Group the data by STATION_ID, Year, Month, and Day, then summarize
merced_data_summarized <- merced_data %>%
  group_by(STATION_ID, Year, Month, Day) %>%
  summarise(
    `Soil Temp 100 mm Avg` = mean(`Soil 100mm Temperature degC`, na.rm = TRUE),
    `Soil Temp 100 mm Min degC` = min(`Soil 100mm Temperature degC`, na.rm = TRUE),
    `Soil Temp 100 mm Max degC` = max(`Soil 100mm Temperature degC`, na.rm = TRUE),
    `Soil Temp 100 mm Range` =  (`Soil Temp 100 mm Max degC` - `Soil Temp 100 mm Min degC`) # temperature range of a day
    
  )

#### Processing temperature data for MVPNR ####

# Summarize on a monthly scale

merced_data_summarized_month <- merced_data_summarized %>%
  group_by(Month) %>%
  summarise(
    MT = mean(`Soil Temp 100 mm Avg`, na.rm = TRUE), # mean temperature of a given month
    MTR = max(`Soil Temp 100 mm Max degC`, na.rm = TRUE) - min(`Soil Temp 100 mm Min degC`, na.rm = TRUE), # temperature range of a given month
    MDTR = mean(`Soil Temp 100 mm Range`, na.rm = TRUE), # mean daily temperature range for each day of a given month
    MTmax = max(`Soil Temp 100 mm Max degC`, na.rm = TRUE), # maximum temperature of a given month
    MMTmax = mean(`Soil Temp 100 mm Max degC`, na.rm = TRUE), # mean maximum temperature of a given month
    MTmin = min(`Soil Temp 100 mm Min degC`, na.rm = TRUE), # minimum temperature of a given month
    MMTmin = mean(`Soil Temp 100 mm Min degC`, na.rm = TRUE), # mean minimum temperature of a given month
    MT_sd = sd(`Soil Temp 100 mm Avg`, na.rm = TRUE), # SD of the mean temperature of a given month
    MTR_sd = sd(`Soil Temp 100 mm Max degC` - `Soil Temp 100 mm Min degC`, na.rm = TRUE), # SD of the mean temperature range of a given month
    MTmax_sd = sd(`Soil Temp 100 mm Max degC`, na.rm = TRUE), # SD of the maximum temperature of a given month
    MTmin_sd = sd(`Soil Temp 100 mm Min degC`, na.rm = TRUE) # SD of the minimum temperature of a given month
  )

# We obtain bioclimatic variables by summarizing the data by Site across years

merced_data_filtered_full <- merced_data_summarized_month %>%
  summarise(
    Seasonality = (sd(MT, na.rm = TRUE)), # SD of the mean annual temperature => seasonality 
    MDTR_sd = sd(MDTR, na.rm = TRUE), # SD of the mean daily temperature range per month of a given year
    MAT = mean(MT, na.rm = TRUE), # mean temperature of a given year (SBIO1)
    MATR = max(MTmax, na.rm = TRUE) - min(MTmin, na.rm = TRUE), # mean temperature annual range (SBIO7) 
    MMTR = mean(MTR, na.rm = TRUE), # mean temperature range per month of a given year (SBIO2)
    MDTR = mean(MDTR, na.rm = TRUE), # mean temperature range per day of a given year (SBIO2)
    MATmax = max(MTmax), # maximum temperature of the hottest month (SBIO5)
    MMTmax = mean(MTmax, na.rm = TRUE), # mean maximum temperature per month of a given year
    MATmin = min(MTmin), # minimum temperature of the coldest month (SBIO6)
    MMTmin = mean(MTmin, na.rm = TRUE), # mean minimum temperature per month of a given year
    Isothermality = (MMTR/MATR)*100 # (SBIO3)
  )

merced_data_filtered_full$Site <- "Merced"
  


#### Calculating offset of air vs soil temperatures in Boyd Deep Canyon (DC) ####
# The pinyon flats site was located near the Boyd Deep Canyon reserve, and we only had access to air temperature data from its coordinates.
# We decided to estimate the soil temperatures by first calculating the offset of air and soil temperatures in Boyd Deep Canyon from Dendra data

# We first read the BDC offset data
offset_data <- fread("air-soil-offset-bdc.csv")

# We then extract Year, Month, and Day from the Time column
offset_data <- offset_data %>%
  mutate(Time = mdy_hm(Time),
         Year = year(Time),
         Month = month(Time),
         Day = day(Time))

# And summarize by daily variation
offset_data_summarized <- offset_data %>%
  group_by(Year, Month, Day) %>%
  summarise(
    `Deep Canyon Soil Temp 50 mm Avg` = mean(`Deep Canyon Soil Temp 50 mm Avg`, na.rm = TRUE),
    `Deep Canyon Soil Temp 50 mm Min degC` = min(`Deep Canyon Soil Temp 50 mm Min degC`, na.rm = TRUE),
    `Deep Canyon Soil Temp 50 mm Max degC` = max(`Deep Canyon Soil Temp 50 mm Max degC`, na.rm = TRUE),
    `Deep Canyon Air Temp Avg degC` = mean(`Deep Canyon Air Temp Avg degC`, na.rm = TRUE),
    `Deep Canyon Air Temp Min degC` = min(`Deep Canyon Air Temp Min degC`, na.rm = TRUE),
    `Deep Canyon Air Temp Max degC` = max(`Deep Canyon Air Temp Max degC`, na.rm = TRUE)
  )

offset_data_summarized <- offset_data_summarized %>% filter(!is.infinite(`Deep Canyon Soil Temp 50 mm Max degC`))

# We then reshape the dataframe to long format, extract Site names, and remove wrong values
offset_data_long <- offset_data_summarized %>%
  pivot_longer(cols = -c(Year, Month, Day), names_to = "Site_Metric", values_to = "Value") %>%
  separate(Site_Metric, into = c("Site", "Metric"), sep = "(?=Soil|Air)", extra = "merge") 


offset_data_summarized_month <- offset_data_summarized %>%
  group_by(Month) %>%
  summarise(
    MT_Soil = mean(`Deep Canyon Soil Temp 50 mm Avg`, na.rm = TRUE),
    MTR_Soil = max(`Deep Canyon Soil Temp 50 mm Max degC`, na.rm = TRUE) - min(`Deep Canyon Soil Temp 50 mm Min degC`, na.rm = TRUE),
    MDTR_Soil = mean(`Deep Canyon Soil Temp 50 mm Max degC` - `Deep Canyon Soil Temp 50 mm Min degC`, na.rm = TRUE),
    MTmax_Soil = max(`Deep Canyon Soil Temp 50 mm Max degC`, na.rm = TRUE),
    MMTmax_Soil = mean(`Deep Canyon Soil Temp 50 mm Max degC`, na.rm = TRUE),
    MTmin_Soil = min(`Deep Canyon Soil Temp 50 mm Min degC`, na.rm = TRUE),
    MMTmin_Soil = mean(`Deep Canyon Soil Temp 50 mm Min degC`, na.rm = TRUE),
    MT_sd_Soil = sd(`Deep Canyon Soil Temp 50 mm Avg`, na.rm = TRUE),
    MTR_sd_Soil = sd(`Deep Canyon Soil Temp 50 mm Max degC` - `Deep Canyon Soil Temp 50 mm Min degC`, na.rm = TRUE),
    MTmax_sd_Soil = sd(`Deep Canyon Soil Temp 50 mm Max degC`, na.rm = TRUE),
    MTmin_sd_Soil = sd(`Deep Canyon Soil Temp 50 mm Min degC`, na.rm = TRUE),
    MT_Air = mean(`Deep Canyon Air Temp Avg degC`, na.rm = TRUE),
    MTR_Air = max(`Deep Canyon Air Temp Max degC`, na.rm = TRUE) - min(`Deep Canyon Air Temp Min degC`, na.rm = TRUE),
    MDTR_Air = mean(`Deep Canyon Air Temp Max degC` - `Deep Canyon Air Temp Min degC`, na.rm = TRUE),
    MTmax_Air = max(`Deep Canyon Air Temp Max degC`, na.rm = TRUE),
    MMTmax_Air = mean(`Deep Canyon Air Temp Max degC`, na.rm = TRUE),
    MTmin_Air = min(`Deep Canyon Air Temp Min degC`, na.rm = TRUE),
    MMTmin_Air = mean(`Deep Canyon Air Temp Min degC`, na.rm = TRUE),
    MT_sd_Air = sd(`Deep Canyon Air Temp Avg degC`, na.rm = TRUE),
    MTR_sd_Air = sd(`Deep Canyon Air Temp Max degC` - `Deep Canyon Air Temp Min degC`, na.rm = TRUE),
    MTmax_sd_Air = sd(`Deep Canyon Air Temp Max degC`, na.rm = TRUE),
    MTmin_sd_Air = sd(`Deep Canyon Air Temp Min degC`, na.rm = TRUE)
  )


# We then calculate the mean offset values between soil and air temperatures for BDC
offset_data_summarized_month$offset_means <- offset_data_summarized_month$MT_Soil - offset_data_summarized_month$MT_Air
offset_data_summarized_month$offset_max <- offset_data_summarized_month$MTmax_Soil - offset_data_summarized_month$MTmax_Air
offset_data_summarized_month$offset_min <- offset_data_summarized_month$MTmin_Soil - offset_data_summarized_month$MTmin_Air
offset_data_summarized_month$Site <- "DC"

# Plot the offset mean by month
offset_data_summarized_month %>% ggplot (aes(x = as.factor(Month),y = offset_means))+
  geom_point(size=6)+
  geom_smooth(aes(group = Site))+
  theme_classic()

# We finally extract the offset means by month to apply to the Pinyon crest data
offset_soil_air_dc <- offset_data_summarized_month %>% dplyr::select(c(Month,offset_means))


#### Loading Pinyon Flats (PF) temperature data ####

# We could only access monthly values of air temperatures for Pinyon Crest from the 
# Boyd Deep Canyon weather station at Pinyon Crest (https://deepcanyon.ucnrs.org/weather-data/

pinyon_crest_data <- fread("Pinyon_crest_temperature_data.csv")

pinyon_crest_data <- pinyon_crest_data %>% dplyr::select(c(Year,Month,`Soil Temp 50 mm Min degC`,`Soil Temp 50 mm Max degC`,`Soil Temp 50 mm Avg`))

pinyon_crest_data <- left_join(pinyon_crest_data,offset_soil_air_dc) # We add information about the temperature offset of soil vs air from a nearby site

# We estimate soil temperatures by adding the offset values 
pinyon_crest_data <- pinyon_crest_data %>%
  mutate(
    `Soil Temp 50 mm Min degC` = `Soil Temp 50 mm Min degC` + offset_means,
    `Soil Temp 50 mm Max degC` = `Soil Temp 50 mm Max degC` + offset_means,
    `Soil Temp 50 mm Avg` = `Soil Temp 50 mm Avg` + offset_means
  ) %>% dplyr::select(!offset_means) 


# We summarize data by month
pinyon_crest_data_summarized_month <- pinyon_crest_data %>%
  group_by(Month) %>%
  summarise(
    MT = mean(`Soil Temp 50 mm Avg`, na.rm = TRUE), # mean temperature of a given month
    MTR = max(`Soil Temp 50 mm Max degC`, na.rm = TRUE) - min(`Soil Temp 50 mm Min degC`, na.rm = TRUE), # temperature range of a given month
    MDTR = mean(`Soil Temp 50 mm Max degC` - `Soil Temp 50 mm Min degC`, na.rm = TRUE), # mean temperature range for each day of a given month
    MTmax = max(`Soil Temp 50 mm Max degC`, na.rm = TRUE), # maximum temperature of a given month
    MMTmax = mean(`Soil Temp 50 mm Max degC`, na.rm = TRUE), # mean maximum temperature of a given month
    MTmin = min(`Soil Temp 50 mm Min degC`, na.rm = TRUE), # minimum temperature of a given month
    MMTmin = mean(`Soil Temp 50 mm Min degC`, na.rm = TRUE), # mean minimum temperature of a given month
    MT_sd = sd(`Soil Temp 50 mm Avg`, na.rm = TRUE), # SD of the mean temperature of a given month
    MTR_sd = sd(`Soil Temp 50 mm Max degC` - `Soil Temp 50 mm Min degC`, na.rm = TRUE), # SD of the mean temperature range of a given month
    MTmax_sd = sd(`Soil Temp 50 mm Max degC`, na.rm = TRUE), # SD of the maximum temperature of a given month
    MTmin_sd = sd(`Soil Temp 50 mm Min degC`, na.rm = TRUE) # SD of the minimum temperature of a given month
  )


pinyon_crest_data_filtered_full <- pinyon_crest_data_summarized_month %>%
  summarise(
    Seasonality = sd(MT, na.rm = TRUE), # SD of the mean annual temperature => seasonality (SBIO4)
    MDTR_sd = sd(MDTR, na.rm = TRUE), # SD of the mean daily temperature range per month of a given year
    MAT = mean(MT, na.rm = TRUE), # mean temperature of a given year (SBIO1)
    MATR = max(MTmax, na.rm = TRUE) - min(MTmin, na.rm = TRUE), # mean temperature annual range (SBIO7) 
    MMTR = mean(MTR, na.rm = TRUE), # mean temperature range per month of a given year (SBIO2)
    MDTR = mean(MDTR, na.rm = TRUE), # mean temperature range per day of a given year (SBIO2)
    MATmax = max(MTmax), # maximum temperature of the hottest month (SBIO5)
    MMTmax = mean(MTmax, na.rm = TRUE), # mean maximum temperature per month of a given year
    MATmin = min(MTmin), # minimum temperature of the coldest month (SBIO6)
    MMTmin = mean(MTmin, na.rm = TRUE), # mean minimum temperature per month of a given year
    MATmax_sd = sd(MTmax, na.rm = TRUE), # SD of the mean maximum temperature per month of a given year
    MATmin_sd = sd(MTmin, na.rm = TRUE), # SD of the mean minimum temperature per month of a given year
    Isothermality = (MMTR/MATR)*100 # (SBIO3)
  )

pinyon_crest_data_filtered_full$Site <- "Pinyon crest"


#### Making a condensed temperature database ####

# Bioclimatic Variables
#SBIO1 annual mean temperature (MAT)
#SBIO2 mean diurnal range (mean of monthly (max temp -  min temp)) (MMTR)  
#SBIO3 isothermality (SBIO2/SBIO7) (×100) (MMTR/MATR) (Isothermality)
#SBIO4 temperature seasonality (standard deviation ×100) (Seasonality)
#SBIO5 max temperature of warmest month (MATmax)
#BIO6 min temperature of coldest month (MATmin)
#BIO7 temperature annual range (SBIO5-SBIO6) (MATR)

temperature_database <- bind_rows(merced_data_filtered_full,dendra_data_summarized_by_site,pinyon_crest_data_filtered_full)%>% dplyr::select(!c(MDTR_sd,MATmax_sd,MATmin_sd))

temperature_database$Site <- factor(temperature_database$Site)
temperature_database$Site_initials <- factor(c("ME","AB","B","DC","EC","SJ","J","M","P","QR","S","V","W","Y","PF"))
temperature_database$`Climate profile` <- c("Hot stable","Hot variable","Cold stable","Hot variable","Hot variable","Cold variable",
                                            "Hot stable","Hot stable","Cold stable","Hot variable","Cold variable","Cold stable",
                                            "Cold variable","Cold stable","Hot variable")
temperature_database$Site <- factor(temperature_database$Site)
levels(temperature_database$Site) <- c("Anza Borrego HQ","Bodega","Boyd Deep Canyon","Elliot Chaparral","James San Jacinto","Jepson Prairie","McLaughlin",
                                       "Merced","Pinyon crest","Point Reyes","Quail Ridge","SNARL","Valentine","White Mountain","Yosemite Mariposa Grove")


write.csv(temperature_database,"temperature_database_tcs_240505.csv",row.names = FALSE)

temperature_database <- fread("temperature_database_tcs_240505.csv")


#### Obtaining monthly temperature averages by site and climate profile for plotting purposes ####

raw_dendra_data <- dendra_data_long %>% dplyr::select(!c(`Soil Temp 50 mm Max degC`,`Soil Temp 50 mm Min degC`,Time,Day)) %>%
  rename("Temperature" = "Soil Temp 50 mm Avg")

raw_merced_data <- merced_data %>% 
  dplyr::select(c(Year,Month,`Soil 100mm Temperature degC`)) %>%
  mutate(Site = "Merced") %>%
  rename("Temperature" = "Soil 100mm Temperature degC")

raw_pinyon_crest_data <- pinyon_crest_data %>% 
  mutate(Site = "Pinyon Crest") %>%
  rename("Temperature" = "Soil Temp 50 mm Avg") %>% 
  dplyr::select(!c(`Soil Temp 50 mm Max degC`,`Soil Temp 50 mm Min degC`))

raw_climate_data <- bind_rows(raw_dendra_data,raw_merced_data,raw_pinyon_crest_data)

summarized_climate_data <- raw_climate_data %>%
  group_by(Year,Month,Site) %>%
  summarize(Temperature = mean(Temperature,na.rm = TRUE)
)

summarized_climate_data <- summarized_climate_data %>% mutate("Site_initials"  = factor(case_when(Site == "Anza Borrego HQ " ~ "AB",
                                                                                                  Site == "Bodega " ~ "B",
                                                                                                  Site == "Deep Canyon " ~ "DC",
                                                                                                  Site == "Elliott " ~ "EC",
                                                                                                  Site == "James " ~ "SJ",
                                                                                                  Site == "Jepson " ~ "J",
                                                                                                  Site == "McLaughlin "  ~ "M",
                                                                                                  Site == "Merced" ~ "ME",
                                                                                                  Site == "Pinyon Crest" ~ "PF",
                                                                                                  Site == "Point Reyes " ~ "P",
                                                                                                  Site == "Quail Ridge " ~ "QR",
                                                                                                  Site == "SNARL " ~ "S",
                                                                                                  Site == "Valentine " ~ "V",
                                                                                                  Site == "WhiteMt Barcroft " ~ "W",
                                                                                                  Site == "Yosemite Mariposa Grove " ~ "Y")))

summarized_climate_data  <- summarized_climate_data %>% mutate("Climate profile"  = factor(case_when(Site_initials %in% c("AB","DC","PF","QR","EC") ~ "Hot variable",
                                                                                              Site_initials %in% c("B", "P","V","Y") ~ "Cold stable",
                                                                                              Site_initials %in% c("SJ","S","W") ~ "Cold variable",
                                                                                              Site_initials %in% c("ME","J","M") ~ "Hot stable")))



summarized_climate_data <- summarized_climate_data %>% filter(!is.na(Temperature))

summarized_climate_data$`Climate profile` <- factor(summarized_climate_data$`Climate profile`, 
                                                    levels = c("Cold stable", "Cold variable",
                                                               "Hot stable", "Hot variable"))



write.csv(summarized_climate_data, "summarized_climate_data.csv", row.names = FALSE)


summarized_climate_data %>% ggplot(aes(x = interaction(Year,Month), y = Temperature, color = Site_initials)) +
  geom_point(size = 3,alpha = 0.6) +
  geom_line(aes(group = interaction(Site,Year)),linewidth = 0.8)+
  scale_color_viridis_d(option = "viridis")+
  theme_classic(base_size = 18,base_family = "Futura Bk BT")+
  theme(
    legend.position = "left",
    panel.background = element_rect(color = "black"), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'))+
  #facet_wrap(~Site_initials,ncol = 5)+
  facet_wrap(~`Climate profile`,ncol = 1)


#### Calculating temperature probability distributions by climate profile ####

# First, we create temperature bins to get a rough estimate of the fraction of time a given Climate Profile falls in each bin
# We first define temperature bins
temperature_bins <- seq(-10, 65, by = 10)
temperature_labels <- cut(temperature_bins[-1], breaks = temperature_bins, include.lowest = TRUE)

all_combinations <- expand_grid(
  `Climate profile` = unique(summarized_climate_data$`Climate profile`),
  Temperature_bin = levels(cut(temperature_bins, breaks = temperature_bins, include.lowest = TRUE))
)

# Then we summarize and fill missing bins with proportion = 0
summarized_temperature_bins <- summarized_climate_data %>%
  mutate(Temperature_bin = cut(Temperature, breaks = temperature_bins, include.lowest = TRUE)) %>%
  group_by(`Climate profile`, Temperature_bin) %>%
  summarise(bin_count = n(), .groups = "drop") %>%
  # We join with all combinations to include missing bins
  right_join(all_combinations, by = c("Climate profile", "Temperature_bin")) %>%
  mutate(bin_count = replace_na(bin_count, 0)) %>%
  # And calculate total counts and proportions by profile
  group_by(`Climate profile`) %>%
  mutate(
    total_count = sum(bin_count, na.rm = TRUE),
    proportion = bin_count / total_count
  ) %>%
  ungroup()

summarized_temperature_bins$Temperature_bin <- factor(summarized_temperature_bins$Temperature_bin,
                                                     levels = c( "[-10,0]", "(0,10]", "(10,20]", "(20,30]", "(30,40]", "(40,50]", "(50,60]" ))


# Now we try to expand the predictions to a continuous distribution 
# We create a function to fit loess and predict proportions
fit_loess_and_predict <- function(data, temp_range) {
  # Fit loess model
  loess_model <- loess(proportion ~ Temperature_mid, data = data, span = 0.5)
  # Predict proportions for the specified temperature range
  predictions <- predict(loess_model, newdata = temp_range)
  # Return results as a tibble
  tibble(
    Temp_continuous = temp_range,
    Proportion_continuous = predictions
  )
}

# We calculate the midpoint of each bin using regular expressions
summarized_temperature_bins <- summarized_temperature_bins %>%
  mutate(
    lower = as.numeric(str_extract(Temperature_bin, "(?<=\\().+?(?=,)")),
    upper = as.numeric(str_extract(Temperature_bin, "(?<=,).+?(?=\\])")),
    Temperature_mid = (lower + upper) / 2
  ) %>%
  dplyr::select(-lower, -upper)  # Remove intermediate columns if not needed

summarized_temperature_bins <- summarized_temperature_bins %>%
  mutate(
    Temperature_mid = if_else(is.na(Temperature_mid), -5, Temperature_mid) # The negative bin needs an additional step
  )

# We then generate continuous proportion data
continuous_proportion_data <- summarized_temperature_bins %>%
  group_by(`Climate profile`) %>%
  nest() %>%  # Nest data by Climate profile
  mutate(
    # We define the temperature range dynamically for each profile, using its minimum and maximum value
    temp_range = purrr::map(data, ~ seq(
      min(.x$Temperature_mid, na.rm = TRUE),
      max(.x$Temperature_mid, na.rm = TRUE),
      by = 0.1
    )),
    # We fit the loess model and generate predictions
    loess_results = map2(data, temp_range, fit_loess_and_predict)
  ) %>%
  unnest(loess_results) %>%
  ungroup() %>%
  filter(!is.na(Proportion_continuous)) %>% 
  dplyr::select(!c(data, temp_range))

# We then plot the smoothed continuous data to visualize the probability distributions
ggplot(continuous_proportion_data, aes(x = Temp_continuous, y = Proportion_continuous, color = `Climate profile`, group = `Climate profile`)) +
  geom_line(size = 0.8) +  # Smoothed line
  geom_ribbon(aes(ymin = 0, ymax = Proportion_continuous, fill = `Climate profile`), alpha = 0.5) +  # Filled area under curve
  labs(
    x = "Temperature (°C)",
    y = "Proportion",
    title = "Temperature Probability Distributions for Climate Profiles"
  ) +
  scale_color_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  theme_classic(base_size = 14, base_family = "Futura Bk BT") +
  theme(
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )+
  facet_wrap(~`Climate profile`,ncol = 1)


write.csv(continuous_proportion_data, "climate_profile_temperature_probability.csv", row.names = FALSE)





#### Loading Precipitation data ########

# We downloaded publicly available precipitation data from the PRISM database. 
# We combined the data from the 15 reserves into one csv file
prism_data <- fread("rainfall_data_PRISM_240725.csv")

#### Processing precipitation data ####

# We reshape the dataframe to long format and extract Site names, as we did with the temperature data
prism_data_long <- prism_data %>%
  pivot_longer(cols = -c(Year, Month), names_to = "Site", values_to = "Monthly Precipitation (mm)")

prism_data_long %>% ggplot(aes(x = as.factor(Month), y = `Monthly Precipitation (mm)`))+
  geom_point()+
  geom_boxplot(alpha=0.6)+
  facet_wrap(~Site)+
  theme_classic()

# We then summarize on a monthly scale
prism_data_summarized_month <- prism_data_long %>%
  group_by(Site, Month) %>%
  summarise(
    MP = mean(`Monthly Precipitation (mm)`, na.rm = TRUE), # mean precipitation of a given month
    MPR = max(`Monthly Precipitation (mm)`, na.rm = TRUE) - min(`Monthly Precipitation (mm)`, na.rm = TRUE), # precipitation range of a given month
    MPmax = max(`Monthly Precipitation (mm)`, na.rm = TRUE), # maximum precipitation of a given month
    MPmin = min(`Monthly Precipitation (mm)`, na.rm = TRUE), # minimum precipitation of a given month
    MP_sd = sd(`Monthly Precipitation (mm)`, na.rm = TRUE) # SD of the mean precipitation of a given month
  )

prism_data_summarized_month %>% ggplot(aes(x = Month, y = MP))+
  geom_point()+
  facet_wrap(~Site)


# We obtain the bioclimatic variables on a yearly basis

# SBIO12 Annual Precipitation (Annual_precipitation)
# SBIO13 precipitation of the wettest month (MAPmax)
# SBIO14 precipitation of the wettest month (MAPmin)
# BIO15 SD of the mean annual precipitation (Precipitation seasonality)

prism_data_filtered_full <- prism_data_summarized_month %>%
  group_by(Site) %>%
  summarise(
    Annual_precipitation = sum(MP), # (SBIO12)
    Precipitation_seasonality = (sd(MP, na.rm = TRUE)), # SD of the mean annual precipitation => (BIO15) seasonality 
    MAPmax = max(MPmax), # precipitation of the wettest month (SBIO13)
    MAPmin = min(MPmin), # precipitation of the driest month (SBIO14)
  )

# And create the final precipitation database

precipitation_database <- prism_data_filtered_full

precipitation_database$Site <- factor(precipitation_database$Site)
precipitation_database$Site_initials <- factor(c("AB","B","DC","EC","SJ","J","M","ME","PF","P","QR","S","V","W","Y"))

write.csv(precipitation_database,"precipitation_database_tcs_250512.csv",row.names = FALSE)

#### Merging the temperature and precipitation databases ####

environment_database <-  left_join(temperature_database,precipitation_database)

write.csv(environment_database,"environment_database_tcs_250512.csv",row.names = FALSE)






#### Calculating PCA of bioclimatic variables ####

environment_database <- fread("environment_database_tcs_250512.csv") %>%
  dplyr::select(!c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10))


library(PCAtools)
library(grid)
library(ggalt)

View(environment_database)


#### S7 - PCA for environmental variables including precipitation data ####


raw_pca_data_env_variables_full <- environment_database %>% 
  #filter(temp %in% c(15,23,28,33,38,43)) %>% # keep only these temps
  dplyr::select(!c(Site,MAPmin)) %>% # select the columns we'll need
  na.omit() %>% # remove rows with NA values 
  as.data.frame() 

raw_pca_data_env_variables_full$Precipitation_seasonality <- (raw_pca_data_env_variables_full$Precipitation_seasonality)*0.01

pca_data_env_variables_full <- raw_pca_data_env_variables_full %>%
  dplyr::select(-Site_initials,-`Climate profile`) %>% # Remove non-numeric columns
  t() # transpose

rownames(raw_pca_data_env_variables_full) <- raw_pca_data_env_variables_full$Site_initials # add row names
colnames(pca_data_env_variables_full) <- raw_pca_data_env_variables_full$Site_initials # add column names
raw_pca_data_env_variables_full$`Climate profile` <- factor(raw_pca_data_env_variables_full$`Climate profile`)

# We perform PCA using PCAtools
pca_result_env_variables_full <- pca(pca_data_env_variables_full, metadata = raw_pca_data_env_variables_full,
                                scale = TRUE,center = TRUE) 


palette_colors_Spectral_Climate <-c("#584B9F","#22C4B3","#A71B4B","#F39B29")

# We construct a biplot coloring data by Order and shaping by climate profile
S7 <- PCAtools::biplot(pca_result_env_variables_full, 
                                        legendPosition = 'none',
                                        #boxedLabels = TRUE,
                                        fillBoxedLoadings = alpha("white", 0.5),
                                        pointSize = 7.5,
                                        alpha = 0.7,
                                        #encircle = TRUE,
                                        #encircleFill = FALSE,
                                        #encircleLineSize = 5,
                                        showLoadings = TRUE,
                                        ntopLoadings = 9,
                                        colLoadingsArrows = "transparent",
                                        drawConnectorsLoadings = TRUE,
                                        colConnectorsLoadings = "white",
                                        lengthLoadingsArrowsFactor = 1.7,
                                        colby = 'Climate profile',  # Color points by 'Order' column
                                        x = "PC1",  # Set x-axis as PC1
                                        y = "PC2",
                                        ylim = -6,
                                        lab = pca_result_env_variables_full$yvars,
                                        boxedLabels = FALSE,
                                        #lab = NULL,
                                        drawConnectors = FALSE,
                                        colConnectors = "transparent"
)+                                                                             
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  #scale_color_viridis_d(option = "inferno") +
  scale_color_manual(values = palette_colors_Spectral_Climate)+
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey60")+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey60")+
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
    panel.background = element_rect(color = "black",fill = "white"),
    plot.background = element_rect(fill = 'white', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )


PCs_env_full <- pca_result_env_variables_full$loadings

write.csv(PCs_env_full, "PCs_env_full.csv", row.names = FALSE)


pc_scores_full <- as.data.frame(pca_result_env_variables_full$rotated)
pc_scores_full$Site_initials <- rownames(pc_scores_full)



environment_database_full <- environment_database %>% left_join(pc_scores_full, by = "Site_initials")



write.csv(environment_database_full,"environment_database_tcs_full_250512.csv",row.names = FALSE)




environment_database_full_long <- environment_database_full %>%
  pivot_longer(cols = c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10),names_to = "PC",values_to = "PC_score")

environment_database_full_long$PC <- factor(environment_database_full_long$PC, levels = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))
environment_database_full_long$Site_initials <- factor(environment_database_full_long$Site_initials, levels = c("Y","B","V","P","W","S","SJ",
                                                                                                                    "M","J","ME","AB","DC","QR","PF","EC"))

f116 <- environment_database_full_long %>%
  ggplot(aes(x = Site_initials, y = PC, fill = PC_score, label = round(PC_score,digits = 2)))+
  #geom_tile(aes(alpha = (1/sd_delta_A_H)))+
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low="#DE5925",mid = "white", high="#0090B5")+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 
f116

environment_database_full_scaled <- environment_database_full %>%
  summarize(
    "Site" = Site,
    "Site_initials" = Site_initials,
    "Climate profile" = `Climate profile`,
    "Seasonality" = scale(Seasonality),              
    "MAT"      = scale(MAT),                
    "MATR"     = scale(MATR),                 
    "MMTR"                     = scale(MMTR),
    "MDTR"                      = scale(MDTR),
    "MATmax"                   = scale(MATmax),
    "MMTmax"                    = scale(MMTmax),
    "MATmin"                   = scale(MATmin),
    "MMTmin"                    = scale(MMTmin),
    "Isothermality"            = scale(Isothermality),
    "Annual_precipitation"     = scale(Annual_precipitation),
    "Precipitation_seasonality" = scale(Precipitation_seasonality),
    "MAPmax"                   = scale(MAPmax)
  )

environment_database_full_long_variables <- environment_database_full_scaled %>%
  pivot_longer(cols = c(Seasonality,MAT,MATR,MMTR,MDTR,MATmax,MMTmax,MATmin,MMTmin,Isothermality,Annual_precipitation,Precipitation_seasonality,MAPmax),names_to = "Bioclimatic Variable",values_to = "Value")

environment_database_full_long_variables$`Bioclimatic Variable` <- factor(environment_database_full_long_variables$`Bioclimatic Variable`, 
                                                                          levels = c(
                                                                            "Annual_precipitation",      
                                                                            "Precipitation_seasonality",
                                                                            "MAPmax",                   
                                                                            "MAT",                       
                                                                            "MMTmin",  
                                                                            "MMTmax",                   
                                                                            "MATmin",                   
                                                                            "MATmax",                    
                                                                            "MATR",                      
                                                                            "MMTR",                      
                                                                            "MDTR",                      
                                                                            "Isothermality",             
                                                                            "Seasonality"   
                                                                          ))
  

environment_database_full_long_variables$Site_initials <- factor(environment_database_full_long_variables$Site_initials, levels = c("Y","B","V","P","W","S","SJ","M",
                                                                                                                                        "J","ME","AB","DC","QR","PF","EC"))

f117 <- environment_database_full_long_variables %>%
  ggplot(aes(x = Site_initials, y = `Bioclimatic Variable`, fill = scale(Value), label = round(Value,digits = 2)))+
  #geom_tile(aes(alpha = (1/sd_delta_A_H)))+
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low="#DE5925",mid = "white", high="#0090B5")+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 
f117



#### PCA for environmental variables without precipitation ####


raw_pca_data_env_variables_noprec <- environment_database %>% 
  #filter(temp %in% c(15,23,28,33,38,43)) %>% # keep only these temps
  dplyr::select(!c(Site, MAPmin,Annual_precipitation,Precipitation_seasonality,MAPmax)) %>% # select the columns we'll need
  na.omit() %>% # remove rows with NA values 
  as.data.frame()


pca_data_env_variables_noprec <- raw_pca_data_env_variables_noprec %>%
  dplyr::select(-Site_initials,-`Climate profile`) %>% # Remove non-numeric columns
  t() # transpose

rownames(raw_pca_data_env_variables_noprec) <- raw_pca_data_env_variables_noprec$Site_initials # add row names
colnames(pca_data_env_variables_noprec) <- raw_pca_data_env_variables_noprec$Site_initials # add column names
raw_pca_data_env_variables_noprec$`Climate profile` <- factor(raw_pca_data_env_variables_noprec$`Climate profile`)

# We perform PCA using PCAtools
pca_result_env_variables_noprec <- pca(pca_data_env_variables_noprec, metadata = raw_pca_data_env_variables_noprec,
                                     scale = TRUE,center = TRUE) 


palette_colors_Spectral_Climate <-c("#584B9F","#22C4B3","#A71B4B","#F39B29")

# We construct a biplot coloring data by Order and shaping by climate profile
F1A.left <- PCAtools::biplot(pca_result_env_variables_noprec, 
                                             legendPosition = 'none',
                                             boxedLabels = FALSE,
                                             fillBoxedLoadings = alpha("white", 0.5),
                                             pointSize = 7.5,
                                             alpha = 0.7,
                                             #encircle = TRUE,
                                             #encircleFill = FALSE,
                                             #encircleLineSize = 5,
                                             showLoadings = TRUE,
                                             ntopLoadings = 9,
                                             colLoadingsArrows = "transparent",
                                             drawConnectorsLoadings = TRUE,
                                             colConnectorsLoadings = "white",
                                             lengthLoadingsArrowsFactor = 1.8,
                                             colby = 'Climate profile',  # Color points by 'Order' column
                                             x = "PC1",  # Set x-axis as PC1
                                             y = "PC2",
                                             ylim = -7,
                                             #xlim = (-7,7),
                                             lab = pca_result_env_variables_noprec$yvars,
                                             #lab = NULL,
                                             drawConnectors = FALSE,
                                             colConnectors = "transparent"
)+                                                                             
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  #scale_color_viridis_d(option = "inferno") +
  scale_color_manual(values = palette_colors_Spectral_Climate)+
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey60")+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey60")+
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
    panel.background = element_rect(color = "black",fill = "white"),
    plot.background = element_rect(fill = 'white', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )


PCs_env_noprec <- pca_result_env_variables_noprec$loadings

write.csv(PCs_env_noprec, "PCs_env_noprec.csv", row.names = FALSE)

View(pca_result_env_variables_noprec$loadings)



pc_scores_noprec <- as.data.frame(pca_result_env_variables_noprec$rotated)
pc_scores_noprec$Site_initials <- rownames(pc_scores_noprec)



environment_database_noprec <- environment_database %>% left_join(pc_scores_noprec, by = "Site_initials")



write.csv(environment_database_noprec,"environment_database_tcs_noprec_250512.csv",row.names = FALSE)

environment_database_noprec_long <- environment_database_noprec %>%
  pivot_longer(cols = c(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10),names_to = "PC",values_to = "PC_score")

environment_database_noprec_long$PC <- factor(environment_database_noprec_long$PC, levels = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))
environment_database_noprec_long$Site_initials <- factor(environment_database_noprec_long$Site_initials, levels = c("Y","B","V","P","W","S","SJ",
                                                                                                                    "M","J","ME","AB","DC","QR","PF","EC"))

f101 <- environment_database_noprec_long %>%
  ggplot(aes(x = Site_initials, y = PC, fill = PC_score, label = round(PC_score,digits = 2)))+
  #geom_tile(aes(alpha = (1/sd_delta_A_H)))+
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low="#DE5925",mid = "white", high="#0090B5")+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 
f101

environment_database_noprec_long_variables <- environment_database_noprec %>%
  pivot_longer(cols = c(Seasonality,MAT,MATR,MMTR,MDTR,MATmax,MMTmax,MATmin,MMTmin,Isothermality),names_to = "Bioclimatic Variable",values_to = "Value")

environment_database_noprec_long_variables$Site_initials <- factor(environment_database_noprec_long_variables$Site_initials, levels = c("Y","B","V","P","W","S","SJ","M",
                                                                                                                                        "J","ME","AB","DC","QR","PF","EC"))

f102 <- environment_database_noprec_long_variables %>%
  ggplot(aes(x = Site_initials, y = `Bioclimatic Variable`, fill = Value, label = round(Value,digits = 2)))+
  #geom_tile(aes(alpha = (1/sd_delta_A_H)))+
  geom_tile()+
  geom_text(size = 3)+
  scale_fill_gradient2(low="#DE5925",mid = "white", high="#0090B5")+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 
f102



#### Loading all chlorophyll data from AppEEARS output files ####
#
#
# We downloaded available temperature values from September 2010 through December 2022
# (Bodega and Point Reyes had water coordinates initially, so we had to add a download for the corrected coordinates) 
#
#

NDVI_500m_1 <- fread("TCS-reserves-MOD13A1-061-results.csv")
NDVI_500m_2 <- fread("TCS-reserves-2-MOD13A1-061-results.csv")

NDVI_500m_1 <- NDVI_500m_1 %>% filter(!ID %in% c("B","P"))

NDVI_500m <- bind_rows(NDVI_500m_1,NDVI_500m_2)

# Create a color palette for all the sites 

full_sites_palette <- c("darkgreen", "violet", "darkgoldenrod", "grey40",
                        "darkturquoise", "goldenrod2", "black",
                        "darkseagreen4", "turquoise2",  "plum2",
                        "darkseagreen2", "turquoise4", "goldenrod4", "grey80","plum4")


#### Processing NDVI data for 15 reserves ####


# We first extract Year, Month, and Day from the Time column
NDVI_500m_full <- NDVI_500m %>%
  mutate(Year = year(Date),
         Month = month(Date),
         Day = day(Date))

# We add a season colum 
NDVI_500m_full_season <- NDVI_500m_full %>%
  mutate(Season = case_when(
    Month %in% c(1,2,3) ~ "Winter",
    Month %in% c(4,5,6) ~ "Spring",
    Month %in% c(7,8,9) ~ "Summer",
    Month %in% c(10,11,12) ~ "Fall")
  )

# We add a season colum 
NDVI_500m_full_season <- NDVI_500m_full_season %>%
  mutate(Month_name = case_when(
    Month %in% c(1) ~ "Jan",
    Month %in% c(2) ~ "Feb",
    Month %in% c(3) ~ "Mar",
    Month %in% c(4) ~ "Apr",
    Month %in% c(5) ~ "May",
    Month %in% c(6) ~ "Jun",
    Month %in% c(7) ~ "Jul",
    Month %in% c(8) ~ "Aug",
    Month %in% c(9) ~ "Sep",
    Month %in% c(10) ~ "Oct",
    Month %in% c(11) ~ "Nov",
    Month %in% c(12) ~ "Dec")
  )

NDVI_500m_full_season$Season <- factor(NDVI_500m_full_season$Season, levels = c("Fall","Winter","Spring","Summer"))
NDVI_500m_full_season$Month_name <- factor(NDVI_500m_full_season$Month_name, levels = c("Jan","Feb","Mar",
                                                                                        "Apr","May","Jun",
                                                                                        "Jul","Aug","Sep",
                                                                                        "Oct","Nov","Dec"))

NDVI_500m_full %>% 
  filter(MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description == "VI produced with good quality") %>%
  ggplot(aes(x = Month,y = MOD13A1_061__500m_16_days_NDVI,color = MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description))+
  geom_point()+
  geom_smooth(method = "loess")+
  facet_wrap(~ID,ncol=5)+
  scale_color_viridis_d(option = "viridis")+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(color = "black",fill = "grey95"),
    plot.background = element_rect(fill = 'grey98', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )



NDVI_500m_full_season %>% 
  filter(MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description == "VI produced with good quality") %>%
  ggplot(aes(x = Season,y = MOD13A1_061__500m_16_days_NDVI,fill = ID))+
  geom_point()+
  geom_boxplot(alpha = 0.8)+
  facet_wrap(~ID,ncol=5)+
  scale_fill_manual(values = full_sites_palette)+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(color = "black",fill = "grey95"),
    plot.background = element_rect(fill = 'grey98', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )





# Add Climate profile data # 

environmental_database <- fread("environment_database_tcs_noprec_250512.csv")

NDVI_500m_full_season <- NDVI_500m_full_season %>% rename(
  "Site_initials" = "ID"
)

environmental_database_with_season <- left_join(NDVI_500m_full_season,environmental_database)


environmental_database_with_season %>% 
  filter(MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description == "VI produced with good quality") %>%
  ggplot(aes(x = Season,y = MOD13A1_061__500m_16_days_NDVI,fill = Site_initials))+
  geom_boxplot(alpha = 0.8)+
  facet_wrap(~`Climate profile`,ncol=4)+
  scale_fill_manual(values = full_sites_palette)+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(color = "black",fill = "grey95"),
    plot.background = element_rect(fill = 'grey98', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

environmental_database_with_season %>% 
  filter(MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description == "VI produced with good quality") %>%
  ggplot(aes(x = Month_name,y = MOD13A1_061__500m_16_days_NDVI,fill = Site_initials))+
  geom_boxplot(alpha = 0.8)+
  facet_wrap(~`Climate profile`,ncol=4)+
  scale_fill_manual(values = full_sites_palette)+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(color = "black",fill = "grey95"),
    plot.background = element_rect(fill = 'grey98', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

environmental_database_with_season %>% 
  filter(MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description == "VI produced with good quality") %>%
  ggplot(aes(x = Month,y = MOD13A1_061__500m_16_days_NDVI,color = Site_initials))+
  geom_point()+
  geom_smooth(method = "loess")+
  facet_wrap(~`Climate profile`,ncol=4)+
  #scale_color_viridis_d(option = "viridis")+
  scale_color_manual(values = full_sites_palette)+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(color = "black",fill = "grey95"),
    plot.background = element_rect(fill = 'grey98', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

environmental_database_with_season %>%
  filter(MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description == "VI produced with good quality") %>%
  mutate(Month = factor(Month, levels = c(7,8,9,10,11,12,1,2,3,4,5,6))) %>%
  ggplot(aes(x = Month, y = MOD13A1_061__500m_16_days_NDVI)) +
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf, fill = "lightsalmon", alpha = 0.3) +  # Summer (Jul–Sep)
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf, fill = "lightyellow", alpha = 0.3) +  # Fall (Oct–Dec)
  annotate("rect", xmin = 6.5, xmax = 9.5, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +    # Winter (Jan–Mar)
  annotate("rect", xmin = 9.5, xmax = 12.5, ymin = -Inf, ymax = Inf, fill = "lightgreen", alpha = 0.3) +  # Spring (Apr–Jun)
  geom_point(alpha = 0.9,aes(color = Site_initials)) +
  geom_smooth(method = "loess", aes(color = Site_initials,group = Site_initials)) +
  facet_wrap(~`Climate profile`, ncol = 4) +
  # scale_color_viridis_d(option = "viridis") +
  scale_color_manual(values = full_sites_palette) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(color = "black", fill = "grey95"),
    plot.background = element_rect(fill = 'grey98', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

environmental_database_with_season$Month


#### Summarize by month ####
# Prepare the data: mean NDVI by Site_initials and Month_name
ndvi_summary <- environmental_database_with_season %>%
  filter(MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description == "VI produced with good quality") %>%
  group_by(Site_initials, Month_name) %>%
  summarise(mean_NDVI = mean(MOD13A1_061__500m_16_days_NDVI, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Month_name, values_from = mean_NDVI) %>%
  replace(is.na(.), 0) %>%
  mutate(across(-Site_initials, rescale))  # Normalize values between 0 and 1 for radar plot

# Ensure consistent Month_nameal column order
Month_name_order <- c("Jan","Feb","Mar",
                      "Apr","May","Jun",
                      "Jul","Aug","Sep",
                      "Oct","Nov","Dec")
ndvi_summary <- ndvi_summary %>%
  dplyr::select(Site_initials, all_of(Month_name_order))

# Create a radar plot for each site
plots <- lapply(unique(ndvi_summary$Site_initials), function(site) {
  ggradar(ndvi_summary %>% filter(Site_initials == site),
          values.radar = element_blank(),
          grid.min = 0, grid.mid = 0.5, grid.max = 1,
          group.line.width = 1,
          group.point.size = 3,
          background.circle.colour = "white",
          gridline.min.linetype = "solid",
          gridline.mid.linetype = "solid",
          gridline.max.linetype = "solid",
          fill=TRUE,
          plot.title = paste(site))
})




radars <- plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],
                    plots[[6]],plots[[7]],plots[[8]],plots[[9]],plots[[10]],
                    plots[[11]],plots[[12]],plots[[13]],plots[[14]],plots[[15]],
                    ncol = 5,
                    align = "h",rel_widths = c(1, 1),axis = "l")
radars

#### Summarize by seasons ####
# Prepare the data: mean NDVI by Site_initials and Season
ndvi_summary <- environmental_database_with_season %>%
  filter(MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description == "VI produced with good quality") %>%
  group_by(Site_initials, Season) %>%
  summarise(mean_NDVI = mean(MOD13A1_061__500m_16_days_NDVI, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Season, values_from = mean_NDVI) %>%
  replace(is.na(.), 0) %>%
  mutate(across(-Site_initials, rescale))  # Normalize values between 0 and 1 for radar plot

# Ensure consistent Seasonal column order
Season_order <- c("Winter","Spring","Summer","Fall")
ndvi_summary <- ndvi_summary %>%
  dplyr::select(Site_initials, all_of(Season_order))

# Create a radar plot for each site
plots <- lapply(unique(ndvi_summary$Site_initials), function(site) {
  ggradar(ndvi_summary %>% filter(Site_initials == site),
          values.radar = element_blank(),
          grid.min = 0, grid.mid = 0.5, grid.max = 1,
          group.line.width = 1,
          group.point.size = 3,
          background.circle.colour = "white",
          gridline.min.linetype = "solid",
          gridline.mid.linetype = "solid",
          gridline.max.linetype = "solid",
          fill=TRUE,
          plot.title = paste(site))
})




radars <- plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],
                    plots[[6]],plots[[7]],plots[[8]],plots[[9]],plots[[10]],
                    plots[[11]],plots[[12]],plots[[13]],plots[[14]],plots[[15]],
                    ncol = 5,
                    align = "h",rel_widths = c(1, 1),axis = "l")
radars

#### Convert dates to julian format ####
colnames(environmental_database_with_season)


environmental_database_with_season$julian_date <- yday(environmental_database_with_season$Date)

#### Obtain derivative ####


ndvi_derivative_df <- environmental_database_with_season %>%
  filter(MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description == "VI produced with good quality") %>%
  group_by(Site_initials, julian_date,Site,`Climate profile`) %>%
  summarise(mean_ndvi = mean(MOD13A1_061__500m_16_days_NDVI, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site_initials) %>%
  arrange(Site_initials, julian_date) %>%
  mutate(
    idx = row_number(),
    n = n(),
    ndvi_derivative = map_dbl(idx, ~{
      # Get window indices, wrapping around if needed
      window_idx <- ((.x - 6 + 1):.x) %% n
      window_idx[window_idx == 0] <- n  # Fix index 0 to n
      dates <- julian_date[window_idx]
      ndvi_vals <- mean_ndvi[window_idx]
      dx <- diff(range(dates))
      dy <- ndvi_vals[6] - ndvi_vals[1]
      dy / dx
    })
  ) %>%
  dplyr::select(-idx, -n)


ndvi_derivative_df %>% 
  ggplot(aes(x = julian_date,y = ndvi_derivative,group = Site_initials,color = Site_initials))+
  annotate("rect", xmin = -Inf, xmax = 90, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +  # Summer (Jul–Sep)
  annotate("rect", xmin = 90, xmax = 181, ymin = -Inf, ymax = Inf, fill = "lightgreen", alpha = 0.3) +  # Fall (Oct–Dec)
  annotate("rect", xmin = 181, xmax = 272, ymin = -Inf, ymax = Inf, fill = "lightsalmon", alpha = 0.3) +    # Winter (Jan–Mar)
  annotate("rect", xmin = 272, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "lightyellow", alpha = 0.3) +  # Spring (Apr–Jun)
  geom_line(linewidth = 0.75, color = "grey60")+
  geom_point(size = 2.5, color = "gray8")+
  geom_point(size = 2)+
  #scale_color_viridis_d(optihttp://127.0.0.1:44673/graphics/ddc1cebe-0bce-4b60-a3d7-2b774c0987f4.pngon = "viridis")+
  facet_wrap(~`Climate profile`,ncol = 5)+
  scale_color_manual(values = full_sites_palette)+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(color = "black",fill = "grey95"),
    plot.background = element_rect(fill = 'grey98', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )


ndvi_derivative_df %>% 
  ggplot(aes(x = julian_date,y = ndvi_derivative,group = Site_initials,color = Site_initials))+
  annotate("rect", xmin = -Inf, xmax = 90, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +  # Summer (Jul–Sep)
  annotate("rect", xmin = 90, xmax = 181, ymin = -Inf, ymax = Inf, fill = "lightgreen", alpha = 0.3) +  # Fall (Oct–Dec)
  annotate("rect", xmin = 181, xmax = 272, ymin = -Inf, ymax = Inf, fill = "lightsalmon", alpha = 0.3) +    # Winter (Jan–Mar)
  annotate("rect", xmin = 272, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "lightyellow", alpha = 0.3) +  # Spring (Apr–Jun)
  geom_hline(yintercept = 0,color = "black",linetype = "dashed",linewidth = 0.5)+
  geom_line(linewidth = 0.75, color = "grey60")+
  geom_point(size = 2.5, color = "gray8")+
  geom_point(size = 2)+
  #scale_color_viridis_d(optihttp://127.0.0.1:44673/graphics/ddc1cebe-0bce-4b60-a3d7-2b774c0987f4.pngon = "viridis")+
  facet_wrap(~Site_initials,ncol = 5)+
  scale_color_manual(values = full_sites_palette)+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(color = "black",fill = "grey95"),
    plot.background = element_rect(fill = 'grey98', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )


ndvi_derivative_df %>% 
  ggplot(aes(x = julian_date,y = mean_ndvi,group = Site_initials,color = Site_initials))+
  annotate("rect", xmin = -Inf, xmax = 90, ymin = -Inf, ymax = Inf, fill = "lightblue", alpha = 0.3) +  # Summer (Jul–Sep)
  annotate("rect", xmin = 90, xmax = 181, ymin = -Inf, ymax = Inf, fill = "lightgreen", alpha = 0.3) +  # Fall (Oct–Dec)
  annotate("rect", xmin = 181, xmax = 272, ymin = -Inf, ymax = Inf, fill = "lightsalmon", alpha = 0.3) +    # Winter (Jan–Mar)
  annotate("rect", xmin = 272, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "lightyellow", alpha = 0.3) +  # Spring (Apr–Jun)
  geom_line(linewidth = 0.75, color = "grey60")+
  geom_point(size = 2.5, color = "gray8")+
  geom_point(size = 2)+
  #scale_color_viridis_d(option = "viridis")+
  facet_wrap(~`Climate profile`,ncol = 5)+
  scale_color_manual(values = full_sites_palette)+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(color = "black",fill = "grey95"),
    plot.background = element_rect(fill = 'grey98', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

ndvi_derivative_df

#### Find maximum derivative ####

ndvi_growth_season <- ndvi_derivative_df %>%
  group_by(Site_initials, Site, `Climate profile`) %>%
  group_modify(~ {
    df <- .
    
    # Find peak derivative and peak NDVI
    peak_index <- which.max(df$ndvi_derivative)
    growth_rate_peak <- max(df$ndvi_derivative, na.rm = TRUE)
    growth_rate_bottom = min(df$ndvi_derivative, na.rm = TRUE) 
    positive_growth_rate = 0 
    growth_rate_peak_julian_date <- df$julian_date[peak_index]
    growth_peak_julian_date <- df$julian_date[which.max(df$mean_ndvi)]
    growth_peak <- max(df$mean_ndvi, na.rm = TRUE)
    
    # Search before peak
    before_peak_indexes <- which(df$ndvi_derivative[1:peak_index] <= positive_growth_rate)
    # Search after peak
    after_peak_indexes <- which(df$ndvi_derivative[(peak_index+1):nrow(df)] <= positive_growth_rate) + peak_index
    
    # Determine season start and end
    if (length(before_peak_indexes) > 0 & length(after_peak_indexes) > 0) {
      growth_season_start <- df$julian_date[max(before_peak_indexes)]
      growth_season_end <- df$julian_date[min(after_peak_indexes)]
    } else if (length(before_peak_indexes) > 0) {
      growth_season_start <- df$julian_date[max(before_peak_indexes)]
      growth_season_end <- df$julian_date[min(before_peak_indexes)]
    } else if (length(after_peak_indexes) > 0) {
      growth_season_start <- df$julian_date[max(after_peak_indexes)]
      growth_season_end <- df$julian_date[min(after_peak_indexes)]
    } else {
      growth_season_start <- NA_real_
      growth_season_end <- NA_real_
    }
    
    tibble(
      growth_peak_julian_date,
      growth_rate_peak_julian_date,
      growth_peak,
      growth_rate_peak,
      growth_season_start,
      growth_season_end
    )
  }) %>%
  ungroup()

# Calculating dates of start and end, and length of the growth season

ndvi_growth_season <- ndvi_growth_season %>%
  mutate(
    growth_season_start_date = as.Date(growth_season_start - 1, origin = "2022-01-01"),
    growth_season_end_date   = as.Date(growth_season_end - 1, origin = "2022-01-01"),
    start_month_day = paste0(as.integer(format(growth_season_start_date, "%m")), "-", as.integer(format(growth_season_start_date, "%d"))),
    end_month_day   = paste0(as.integer(format(growth_season_end_date, "%m")), "-", as.integer(format(growth_season_end_date, "%d"))),
    start_month = as.integer(format(growth_season_start_date, "%m")),
    end_month   = as.integer(format(growth_season_end_date, "%m"))
  ) %>%
  mutate(growth_season_length = case_when(
    growth_season_start > growth_season_end ~ (365 - growth_season_start + growth_season_end),
    growth_season_start < growth_season_end ~ growth_season_end - growth_season_start))


#### Subset dates of temperature dataset ####

# We want to keep only the growth season, but handling the cases in which we need to wrap around the end of year
summarized_climate_data <- fread("summarized_climate_data.csv")

summarized_climate_data_growing_season <- summarized_climate_data %>%
  left_join(ndvi_growth_season %>% dplyr::select(Site_initials, start_month, end_month,growth_season_length),
            by = "Site_initials") %>%
  filter(
    # Handle non-wrapping seasons normally (e.g., April to September)
    (start_month <= end_month & Month >= start_month & Month <= end_month) |
      # Handle wrapping seasons (e.g., November to March)
      (start_month > end_month & (Month >= start_month | Month <= end_month))
  )

write.csv(summarized_climate_data_growing_season,"summarized_climate_data_growing_season.csv", row.names = FALSE )


summarized_climate_data_growing_season <- fread("summarized_climate_data_growing_season.csv")

#### Calculating temperature probability distributions by climate profile ####

# First, we create temperature bins to get a rough estimate of the fraction of time a given Climate Profile falls in each bin
# We first define temperature bins
temperature_bins <- seq(-10, 65, by = 5)
temperature_labels <- cut(temperature_bins[-1], breaks = temperature_bins, include.lowest = TRUE)

all_combinations <- expand_grid(
  `Climate profile` = unique(summarized_climate_data_growing_season$`Climate profile`),
  Temperature_bin = levels(cut(temperature_bins, breaks = temperature_bins, include.lowest = TRUE))
)

# Then we summarize and fill missing bins with proportion = 0
summarized_temperature_bins <- summarized_climate_data_growing_season %>%
  mutate(Temperature_bin = cut(Temperature, breaks = temperature_bins, include.lowest = TRUE)) %>%
  group_by(`Climate profile`, Temperature_bin) %>%
  summarise(bin_count = n(), .groups = "drop") %>%
  # We join with all combinations to include missing bins
  right_join(all_combinations, by = c("Climate profile", "Temperature_bin")) %>%
  mutate(bin_count = replace_na(bin_count, 0)) %>%
  # And calculate total counts and proportions by profile
  group_by(`Climate profile`) %>%
  mutate(
    total_count = sum(bin_count, na.rm = TRUE),
    proportion = bin_count / total_count
  ) %>%
  ungroup()

summarized_temperature_bins$Temperature_bin <- factor(summarized_temperature_bins$Temperature_bin,
                                                      levels = c( "[-10,-5]","(-5,0]", "(0,5]", "(10,15]" ,"(15,20]", "(20,25]", 
                                                                  "(25,30]", "(30,35]", "(35,40]", "(40,45]",
                                                                  "(45,50]", "(5,10]", "(50,55]", "(55,60]", "(60,65]" ))


# Now we try to expand the predictions to a continuous distribution 
# We create a function to fit loess and predict proportions
fit_loess_and_predict <- function(data, temp_range) {
  # Fit loess model
  loess_model <- loess(proportion ~ Temperature_mid, data = data, span = 0.5)
  # Predict proportions for the specified temperature range
  predictions <- predict(loess_model, newdata = temp_range)
  # Return results as a tibble
  tibble(
    Temp_continuous = temp_range,
    Proportion_continuous = predictions
  )
}

# We calculate the midpoint of each bin using regular expressions
summarized_temperature_bins <- summarized_temperature_bins %>%
  mutate(
    lower = as.numeric(str_extract(Temperature_bin, "(?<=\\().+?(?=,)")),
    upper = as.numeric(str_extract(Temperature_bin, "(?<=,).+?(?=\\])")),
    Temperature_mid = (lower + upper) / 2
  ) %>%
  dplyr::select(-lower, -upper)  # Remove intermediate columns if not needed

summarized_temperature_bins <- summarized_temperature_bins %>%
  mutate(
    Temperature_mid = if_else(is.na(Temperature_mid), -5, Temperature_mid) # The negative bin needs an additional step
  )

# We then generate continuous proportion data
continuous_proportion_data <- summarized_temperature_bins %>%
  group_by(`Climate profile`) %>%
  nest() %>%  # Nest data by Climate profile
  mutate(
    # We define the temperature range dynamically for each profile, using its minimum and maximum value
    temp_range = purrr::map(data, ~ seq(
      min(.x$Temperature_mid, na.rm = TRUE),
      max(.x$Temperature_mid, na.rm = TRUE),
      by = 0.1
    )),
    # We fit the loess model and generate predictions
    loess_results = map2(data, temp_range, fit_loess_and_predict)
  ) %>%
  unnest(loess_results) %>%
  ungroup() %>%
  filter(!is.na(Proportion_continuous)) %>% 
  dplyr::select(!c(data, temp_range))

# We remove negative probability values
continuous_proportion_data <- continuous_proportion_data %>%
  mutate(
    "Proportion_continuous" = ifelse(Proportion_continuous < 0, 0, Proportion_continuous)
  )


# We then plot the smoothed continuous data to visualize the probability distributions
ggplot(continuous_proportion_data, aes(x = Temp_continuous, y = Proportion_continuous, color = `Climate profile`, group = `Climate profile`)) +
  geom_line(size = 0.8) +  # Smoothed line
  geom_ribbon(aes(ymin = 0, ymax = Proportion_continuous, fill = `Climate profile`), alpha = 0.5) +  # Filled area under curve
  labs(
    x = "Temperature (°C)",
    y = "Proportion",
    title = "Temperature Probability Distributions for Climate Profiles"
  ) +
  scale_color_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  theme_classic(base_size = 14, base_family = "Futura Bk BT") +
  theme(
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )+
  facet_wrap(~`Climate profile`,ncol = 1)


write.csv(continuous_proportion_data, "climate_profile_temperature_probability_growing_season.csv", row.names = FALSE)


#### Calculating temperature probability distributions by site ####

# First, we create temperature bins to get a rough estimate of the fraction of time a given Climate Profile falls in each bin
# We first define temperature bins
temperature_bins <- seq(-10, 65, by = 5)
temperature_labels <- cut(temperature_bins[-1], breaks = temperature_bins, include.lowest = TRUE)

all_combinations <- expand_grid(
  Site_initials = unique(summarized_climate_data_growing_season$Site_initials),
  Temperature_bin = levels(cut(temperature_bins, breaks = temperature_bins, include.lowest = TRUE))
)

all_combinations <- all_combinations %>%
  mutate(`Climate profile` = 
           case_when(Site_initials %in% c("W","S","SJ") ~ "Cold variable",
                     Site_initials %in% c("Y","B","V","P") ~ "Cold stable",
                     Site_initials %in% c("AB","DC","QR","PF","EC") ~ "Hot variable",
                     Site_initials %in% c("M","J","ME") ~ "Hot stable"))


# Then we summarize and fill missing bins with proportion = 0
summarized_temperature_bins_site <- summarized_climate_data_growing_season %>%
  mutate(Temperature_bin = cut(Temperature, breaks = temperature_bins, include.lowest = TRUE)) %>%
  group_by(Site_initials, Temperature_bin, `Climate profile`) %>%
  summarise(bin_count = n(), .groups = "drop") %>%
  # We join with all combinations to include missing bins
  right_join(all_combinations, by = c("Site_initials","Climate profile", "Temperature_bin")) %>%
  mutate(bin_count = replace_na(bin_count, 0)) %>%
  # And calculate total counts and proportions by profile
  group_by(Site_initials) %>%
  mutate(
    total_count = sum(bin_count, na.rm = TRUE),
    proportion = bin_count / total_count
  ) %>%
  ungroup()

summarized_temperature_bins_site$Temperature_bin <- factor(summarized_temperature_bins_site$Temperature_bin,
                                                           levels = c( "[-10,-5]" ,"(-5,0]", "(0,5]", "(10,15]" ,"(15,20]", "(20,25]", 
                                                                       "(25,30]", "(30,35]", "(35,40]", "(40,45]",
                                                                       "(45,50]", "(5,10]", "(50,55]", "(55,60]", "(60,65]"))


# Now we try to expand the predictions to a continuous distribution 
# We create a function to fit loess and predict proportions
fit_loess_and_predict <- function(data, temp_range) {
  # Fit loess model
  loess_model <- loess(proportion ~ Temperature_mid, data = data, span = 0.5)
  # Predict proportions for the specified temperature range
  predictions <- predict(loess_model, newdata = temp_range)
  # Return results as a tibble
  tibble(
    Temp_continuous = temp_range,
    Proportion_continuous = predictions
  )
}

# We calculate the midpoint of each bin using regular expressions
summarized_temperature_bins_site <- summarized_temperature_bins_site %>%
  mutate(
    lower = as.numeric(str_extract(Temperature_bin, "(?<=\\().+?(?=,)")),
    upper = as.numeric(str_extract(Temperature_bin, "(?<=,).+?(?=\\])")),
    Temperature_mid = (lower + upper) / 2
  ) %>%
  dplyr::select(-lower, -upper)  # Remove intermediate columns if not needed

summarized_temperature_bins_site <- summarized_temperature_bins_site %>%
  mutate(
    Temperature_mid = if_else(is.na(Temperature_mid), -5, Temperature_mid) # The negative bin needs an additional step
  )

# We then generate continuous proportion data
continuous_proportion_data_site <- summarized_temperature_bins_site %>%
  group_by(Site_initials, `Climate profile`) %>%
  nest() %>%  # Nest data by Site_initials
  mutate(
    # We define the temperature range dynamically for each profile, using its minimum and maximum value
    temp_range = purrr::map(data, ~ seq(
      min(.x$Temperature_mid, na.rm = TRUE),
      max(.x$Temperature_mid, na.rm = TRUE),
      by = 0.1
    )),
    # We fit the loess model and generate predictions
    loess_results = map2(data, temp_range, fit_loess_and_predict)
  ) %>%
  unnest(loess_results) %>%
  ungroup() %>%
  filter(!is.na(Proportion_continuous)) %>% 
  dplyr::select(!c(data, temp_range))


continuous_proportion_data_site <- continuous_proportion_data_site %>%
  mutate(
    "Proportion_continuous" = ifelse(Proportion_continuous < 0, 0, Proportion_continuous)
  )


# We then plot the smoothed continuous data to visualize the probability distributions
ggplot(continuous_proportion_data_site, aes(x = Temp_continuous, y = Proportion_continuous, color = `Climate profile`, group = Site_initials)) +
  geom_line(size = 0.8) +  # Smoothed line
  geom_ribbon(aes(ymin = 0, ymax = Proportion_continuous, fill = `Climate profile`), alpha = 0.5) +  # Filled area under curve
  labs(
    x = "Temperature (°C)",
    y = "Proportion",
    title = "Temperature Probability Distributions for Site_initialss"
  ) +
  scale_color_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  theme_classic(base_size = 14, base_family = "Futura Bk BT") +
  theme(
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )+
  facet_wrap(~Site_initials,ncol = 3)


write.csv(continuous_proportion_data_site, "site_temperature_probability_growing_season.csv", row.names = FALSE)



#### Calculating more frequent temperatures and averages of the growth season by climate profile ####

continuous_proportion_data_growing_season <- fread("climate_profile_temperature_probability_growing_season.csv")

continuous_proportion_data_growing_season_summarized <- continuous_proportion_data_growing_season %>%
  group_by(`Climate profile`) %>%
  summarize(more_likely_temp = Temp_continuous[which.max(Proportion_continuous)],
            mean_temp = weighted.mean(Temp_continuous,Proportion_continuous))

continuous_proportion_data <- fread("climate_profile_temperature_probability.csv")

continuous_proportion_data_summarized <- continuous_proportion_data %>%
  group_by(`Climate profile`) %>%
  summarize(more_likely_temp = Temp_continuous[which.max(Proportion_continuous)],
            mean_temp = weighted.mean(Temp_continuous,Proportion_continuous))




