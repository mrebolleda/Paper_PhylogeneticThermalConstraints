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


#### Loading all precipitation data from dendra database ####

# We downloaded available temperature values from September 2010 through December 2022

dendra_precipitation_data <- fread("rainfall_data_ucnrs_240722.csv")[,c(1,15:27)] 
prism_precipitation_data <- fread("rainfall_data_PRISM_240725.csv")

dendra_precipitation_data_long <- dendra_precipitation_data %>%
  pivot_longer(cols = c(!time),names_to = "Site",values_to = "precipitation")%>%
  filter(precipitation > -0.1 & precipitation < 3000)%>%
  mutate(Year = year(time),
         Month = month(time),
         Day = day(time))%>%
  mutate(part_of_year = case_when(
    Month < 10 ~ "beginning",
    Month >= 10 ~ "end"
  ))

# Getting annual precipitation values from dendra
dendra_precipitation_data_long_sum <- dendra_precipitation_data_long %>%
  group_by(Site,Year)%>%
  summarize(
    "Annual_precipitation" = sum(c(max(precipitation[part_of_year == "beginning"],na.rm = TRUE),
                                   max(precipitation[part_of_year == "end"],na.rm = TRUE)
  )
)
)

dendra_precipitation_data_long %>%
  ggplot(aes(x = time,y = precipitation,color = Site,group = Site))+
  geom_line()


dendra_precipitation_data_long_sum %>%
  ggplot(aes(x = factor(Year),y = Annual_precipitation,color = Site,group = Site))+
  geom_line()+
  facet_grid(~Site)


dendra_precipitation_data_long_sum <- dendra_precipitation_data_long_sum %>%
  mutate(Site = case_when(
    Site == "anza-borrego-rainfall-cumulative-mm" ~ "Anza Borrego HQ",
    Site == "bodega-rainfall-cumulative-mm" ~ "Bodega",
    Site == "deep-canyon-rainfall-cumulative-mm" ~ "Boyd Deep Canyon",
    Site == "elliott-rainfall-cumulative-mm" ~ "Elliot Chaparral",
    Site == "james-precipitation-cumulative-mm" ~ "James San Jacinto",
    Site == "jepson-rainfall-cumulative-mm" ~ "Jepson Prairie", 
    Site == "mclaughlin-rainfall-cumulative-mm" ~ "McLaughlin",
    Site == "point-reyes-rainfall-cumulative-mm" ~ "Point Reyes",
    Site == "quail-ridge-rainfall-cumulative-mm" ~ "Quail Ridge",
    Site == "snarl-rainfall-cumulative-mm" ~ "SNARL" ,
    Site == "valentine-precipitation-cumulative-mm" ~ "Valentine",
    Site == "whitemt-barcroft-rainfall-cumulative-mm" ~ "White Mountain",
    Site == "yosemite-mariposa-grove-rainfall-cumulative-mm" ~ "Yosemite Mariposa Grove"
  ))

dendra_precipitation_data_long_sum <- dendra_precipitation_data_long_sum %>%
  mutate("database" = "dendra")

#### Processing PRISM data ####

prism_precipitation_data_sum <- prism_precipitation_data %>%
  pivot_longer(cols = !c(Year,Month),names_to = "Site",values_to = "precipitation")%>%
  group_by(Site,Year)%>%
  summarize("Annual_precipitation" = sum(precipitation))%>%
  mutate("database" = "prism")

#### Merging dataframes ####

precipitation_database_comparison <- bind_rows(prism_precipitation_data_sum,dendra_precipitation_data_long_sum) %>%
  filter(Annual_precipitation > -0.1 & Annual_precipitation < 2000)

filtered_precipitation_database_comparison <- precipitation_database_comparison %>%
  group_by(Site, Year) %>%
  filter(sum(!is.na(Annual_precipitation)) == 2) %>%
  ungroup()%>%
  group_by(database,Site) %>%
  filter(sum(!is.na(Annual_precipitation)) > 5)%>%
  ungroup()

#### S1 - Plotting the comparison between precipitation measurements ####
S1 <- filtered_precipitation_database_comparison %>%
  ggplot(aes(x = factor(Year), y = Annual_precipitation, color = database,group = database))+
  geom_line()+
  geom_point(size = 2.5)+
  xlab("Year")+
  ylab("Annual Precipitation")+
  facet_grid(~Site,scales = "free_x")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent'),
    axis.text.x = element_text(angle = 60,hjust = 1)
  ) 
