######################################################
# Paper: Conserved upper thermal limits and small   #
# safety margins in soil copiotrophic bacteria.     #

###########################################################
####### Obtain density distributions of temperature ########
##########################################################
## Code: Ariel Favier (afavier@uci.edu)
## Last updated: Aug 19 2026
#######################################################

# Load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(viridis)
library(viridisLite)
font_import("Futura")
loadfonts("win")
fonts()
library(grid)
library(ggalt)

# # Setting the directory 
# parent <- ("YOUR PATH/Paper_PhylogeneticThermalConstraints/") set the path to the parent directory
paste(parent,"01_TCS_bioclimatic_variables/", sep="") %>% setwd

# Open data
summarized_climate_data <- fread("03_Output/summarized_climate_data.csv")



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


#write.csv(continuous_proportion_data, "climate_profile_temperature_probability.csv", row.names = FALSE)



#### Loading all chlorophyll data from AppEEARS output files ####
#
#
# We downloaded available temperature values from September 2010 through December 2022
# (Bodega and Point Reyes had water coordinates initially, so we had to add a download for the corrected coordinates) 
#
#

NDVI_500m_1 <- fread("02_Data/TCS-reserves-MOD13A1-061-results.csv")
NDVI_500m_2 <- fread("02_Data/TCS-reserves-2-MOD13A1-061-results.csv")

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

environmental_database <- fread("03_Output/environment_database_tcs_noprec_250512.csv")

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


write.csv(continuous_proportion_data, "03_Output/climate_profile_temperature_probability_growing_season.csv", row.names = FALSE)


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


write.csv(continuous_proportion_data_site, "03_Output/site_temperature_probability_growing_season.csv", row.names = FALSE)



#### Calculating more frequent temperatures and averages of the growth season by climate profile ####

continuous_proportion_data_growing_season <- fread("03_Output/climate_profile_temperature_probability_growing_season.csv")

continuous_proportion_data_growing_season_summarized <- continuous_proportion_data_growing_season %>%
  group_by(`Climate profile`) %>%
  summarize(more_likely_temp = Temp_continuous[which.max(Proportion_continuous)],
            mean_temp = weighted.mean(Temp_continuous,Proportion_continuous))

continuous_proportion_data <- fread("03_Output/climate_profile_temperature_probability.csv")

continuous_proportion_data_summarized <- continuous_proportion_data %>%
  group_by(`Climate profile`) %>%
  summarize(more_likely_temp = Temp_continuous[which.max(Proportion_continuous)],
            mean_temp = weighted.mean(Temp_continuous,Proportion_continuous))




