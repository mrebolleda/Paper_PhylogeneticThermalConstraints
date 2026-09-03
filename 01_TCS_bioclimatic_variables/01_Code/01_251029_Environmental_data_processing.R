######################################################
# Paper: Conserved upper thermal limits and small   #
# safety margins in soil copiotrophic bacteria.     #

#######################################################
####### Obtain and process environmental data ########
#######################################################
## Code: Ariel Favier (afavier@uci.edu)
## Last updated: Aug 19 2026
#######################################################


# Here we obtain climate series data for the 15 sites sampled
# We calculate 
# 1) temperature and precipitation soil bioclimatic variables
# 2) annual temperature distributions per site and climate profile
# 3) plant growth season 
# 4) temperature distributions per site and climate profile during the growth season
# 5) growth season temperature distributions per site and climate profile

# Load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(viridis)
library(viridisLite)
font_import("Futura")
loadfonts("win")
fonts()

# Setting the directory 
# parent <- ("YOUR PATH/Paper_PhylogeneticThermalConstraints/") set the path to the parent directory
paste(parent,"01_TCS_bioclimatic_variables/", sep="") %>% setwd

#### Loading all temperature data from dendra database ####

# We downloaded available temperature values from September 2010 through December 2022
dendra_data <- fread("02_Data/temperature-data-ucnrs-240502.csv.gz")

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

#fwrite(dendra_data_long, "03_Output/dendra_data_long.csv.gz")
#dendra_data_long <- fread("03_Output/dendra_data_long.csv.gz")

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

merced_data <- fread("02_Data/UCM_temperature_data.csv.gz")

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
offset_data <- read.csv("02_Data/air-soil-offset-bdc.csv.gz")

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
    `Deep Canyon Soil Temp 50 mm Avg` = mean(`Deep.Canyon.Soil.Temp.50.mm.Avg`, na.rm = TRUE),
    `Deep Canyon Soil Temp 50 mm Min degC` = min(`Deep.Canyon.Soil.Temp.50.mm.Min.degC`, na.rm = TRUE),
    `Deep Canyon Soil Temp 50 mm Max degC` = max(`Deep.Canyon.Soil.Temp.50.mm.Max.degC`, na.rm = TRUE),
    `Deep Canyon Air Temp Avg degC` = mean(`Deep.Canyon.Air.Temp.Avg.degC`, na.rm = TRUE),
    `Deep Canyon Air Temp Min degC` = min(`Deep.Canyon.Air.Temp.Min.degC`, na.rm = TRUE),
    `Deep Canyon Air Temp Max degC` = max(`Deep.Canyon.Air.Temp.Max.degC`, na.rm = TRUE)
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

pinyon_crest_data <- fread("02_Data/Pinyon_crest_temperature_data.csv.gz")

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


write.csv(temperature_database,"03_Output/temperature_database_tcs_240505.csv",row.names = FALSE)

#temperature_database <- fread("03_Output/temperature_database_tcs_240505.csv")


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



write.csv(summarized_climate_data, "03_Output/summarized_climate_data.csv", row.names = FALSE)


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



#### Loading Precipitation data ########

# We downloaded publicly available precipitation data from the PRISM database. 
# We combined the data from the 15 reserves into one csv file
prism_data <- fread("02_Data/rainfall_data_PRISM_240725.csv.gz")

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

fwrite(precipitation_database,"03_Output/precipitation_database_tcs_250512.csv.gz",row.names = FALSE)

#### Merging the temperature and precipitation databases ####

environment_database <-  left_join(temperature_database,precipitation_database)

fwrite(environment_database,"03_Output/environment_database_tcs_250512.csv.gz")
