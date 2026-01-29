
###### Introduction ######
#
#
# This script is used to generate the main and supplementary figures of the paper
# All needed data is provided in the form of raw files and intermediate, pre-processed files, 
# for parts of the code that are computationally intensive.
# Additionally, output figures and datasets are also provided.
# To run any part of the script, the following three sections should be run first
# Main sections are indicated with six # symbols, as ###### Name of the section 
# Each main section has the necessary datafiles to be run independently from the rest of the script
# listed at the beginning
#
#

###### Setting the directory and packages ######
#
#
# These are all the required packages 
#
#
setwd("C:\\Users\\Ariel\\OneDrive - personalmicrosoftsoftware.uci.edu\\Desktop\\Growthcurves_4\\5) TCS_figures")

#setwd("PATH_TO_YOUR_DIRECTORY")
library(extrafont)
library(data.table)
library(purrr)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(MuMIn)
library(vegan)
library(ape)
library(cowplot)
library(lme4)
library(ggplot2)
library(sf)
library(scico)
library(ggtree)
library(phytools)
library(geiger)
library(PCAtools)
library(grid)
library(colorspace)
library(emmeans)
library(corrplot)
library(ggradar)
library(forcats)
library(ggside)
library(partR2)
library(rsq)
library(lmerTest)  
library(stringr)
###### Loading all data for thermal traits ######
#
#
# Listed below are all the files needed to run the script continuously from this point. 
# Necessary intermediate files will be generated in the following steps, 
# but are also provided at the start of each main section to save time
#
#
# TPC parameter output file of the 250324_TCS_TPC.R script
ave_params_gcplyr_combined <- fread("251029_tcs_params_gcplyr.csv")

# TPC fits output file of the 250324_TCS_TPC.R script (we replace negative values with 0)
ave_preds_gcplyr_combined <- fread("251029_tcs_fits_gcplyr.csv")


ave_results_gcplyr <- list(
  params = ave_params_gcplyr_combined,
  preds  = ave_preds_gcplyr_combined
)

site_levels <- c(
  "B","P","V","Y",
  "SJ","S","W",
  "ME","J","M",
  "AB","DC","QR","PF","EC"
)

climate_levels <- c(
  "Cold stable", "Cold variable",
  "Hot stable",  "Hot variable"
)

tax_cols <- c("Domain", "Class", "Order", "Family", "Genus", "Species")

# Here we convert the climate profile variable to a factor and sort the levels, 
# we do the same with the site variable, 
# and then we convert the taxonomy columns to factors as well

ave_results_gcplyr <- ave_results_gcplyr %>%
  purrr::map(~ .x %>%
               mutate(
                 Site_initials     = factor(Site_initials, levels = site_levels),
                 `Climate profile` = factor(`Climate profile`, levels = climate_levels)
               ) %>%
               mutate(across(any_of(tax_cols), as.factor))
  )

list2env(
  list(
    ave_params_gcplyr_combined = ave_results_gcplyr$params,
    ave_preds_gcplyr_combined  = ave_results_gcplyr$preds
  ),
  .GlobalEnv
)

# We clip the climate profile column with mapped seqIDs for future use
climate_profile_vector <- ave_params_gcplyr_combined %>% dplyr::select(c(seqID,`Climate profile`))

# We remove negative values of TPCs for plotting purposes
ave_preds_gcplyr_combined <- ave_preds_gcplyr_combined %>% 
  mutate(.fitted = ifelse(.fitted <= 0, 0, .fitted))

# For some analyses, it is convenient to have an indicator of the loss of fitness with every degree of temperature increase 
# above the thermal optimum, here defined as "superopt slope"
# along with its counterpart increase in fitness under the thermal optimum, "subopt_slope"
# It is worth noting that these are simplified versions of the Ea and EH parameters, energy of activation and deactivation
# Obtained after transforming the TPC into a triangle
# We first calculate the slope of performance decrease above Topt
ave_params_gcplyr_combined$superopt_slope <- -1*(ave_params_gcplyr_combined$rmax)/((ave_params_gcplyr_combined$ctmax)-(ave_params_gcplyr_combined$topt))
# And then we calculate the slope of performance increase below Topt
ave_params_gcplyr_combined$subopt_slope <- (ave_params_gcplyr_combined$rmax)/((ave_params_gcplyr_combined$topt)-(ave_params_gcplyr_combined$ctmin))



# We then load the taxonomic assignment database obtained as output of the "TCS_genome_quality_assessment.R" script
taxonomy <-  fread("Kbase_taxonomy_quality_filtered.csv")

taxonomy <- taxonomy %>% dplyr::select(c(seqID,Domain,Phylum,Class,Order,Family,Genus,Species))

ave_params_gcplyr_combined <- ave_params_gcplyr_combined %>% filter(seqID %in% taxonomy$seqID)
ave_preds_gcplyr_combined <- ave_preds_gcplyr_combined %>% filter(seqID %in% taxonomy$seqID)

taxonomy$seqID <- as.character(taxonomy$seqID)



# Our sampling was very enriched in some taxa, which can impact the significance of figures and, more importantly,
# statistical analyses. Here we take a preliminary step to only keep relatively abundant taxa, as those with at least
# >1% of the total collection. In this case, the Class Gammaproteobacteria concentrates the four more abundant Orders

ave_params_gcplyr_combined_common_taxa <- ave_params_gcplyr_combined %>% 
  filter(Class %in% c("Gammaproteobacteria"))%>%
  droplevels()

ave_preds_gcplyr_combined_common_taxa <- ave_preds_gcplyr_combined %>% 
  filter(Class %in% c("Gammaproteobacteria"))%>%
  droplevels()

# For practical purposes, we then split the original datasets into one per growth metric
ave_params_gcplyr_combined_common_taxa_auc <- ave_params_gcplyr_combined_common_taxa %>% 
  filter(metric == "auc")
ave_params_gcplyr_combined_common_taxa_rate <- ave_params_gcplyr_combined_common_taxa %>% 
  filter(metric == "rate")

ave_preds_gcplyr_combined_common_taxa_auc <- ave_preds_gcplyr_combined_common_taxa %>% 
  filter(metric == "auc")
ave_preds_gcplyr_combined_common_taxa_rate <- ave_preds_gcplyr_combined_common_taxa %>% 
  filter(metric == "rate")



# The filtering process outlined above does not completely remove rare taxa, and thus later
# we focus our analyses on the two majoritary taxa, Orders Enterobacterales and Pseudomonadales
# which comprise >90% of the collection

ave_params_gcplyr_combined_major_taxa <- ave_params_gcplyr_combined %>% 
  filter(Order %in% c("Enterobacterales","Pseudomonadales")) %>%
  droplevels()

ave_preds_gcplyr_combined_major_taxa <- ave_preds_gcplyr_combined %>% 
  filter(Order %in% c("Enterobacterales","Pseudomonadales")) %>%
  droplevels()


# For practical purposes, we also split the original datasets into one per growth metric
ave_params_gcplyr_combined_major_taxa_auc <- ave_params_gcplyr_combined_major_taxa %>% 
  filter(metric == "auc")
ave_params_gcplyr_combined_major_taxa_rate <- ave_params_gcplyr_combined_major_taxa %>% 
  filter(metric == "rate")

ave_preds_gcplyr_combined_major_taxa_auc <- ave_preds_gcplyr_combined_major_taxa %>% 
  filter(metric == "auc")
ave_preds_gcplyr_combined_major_taxa_rate <- ave_preds_gcplyr_combined_major_taxa %>% 
  filter(metric == "rate")


# We will visualize time series of temperature by site and climate profile in some figures
# We will also visualize the probability of different temperatures across climate profiles
# The following two datasets are obtained as output of the "Environmental_data_processing.R" script

summarized_climate_data <- fread("summarized_climate_data.csv")
continuous_proportion_data <- fread("climate_profile_temperature_probability.csv")

# Here we split the temperature time series by climate profile
cold_stable_climate_time_series <- summarized_climate_data %>% filter(`Climate profile` == "Cold stable")
hot_stable_climate_time_series <- summarized_climate_data %>% filter(`Climate profile` == "Hot stable")
cold_variable_climate_time_series <- summarized_climate_data %>% filter(`Climate profile` == "Cold variable")
hot_variable_climate_time_series <- summarized_climate_data %>% filter(`Climate profile` == "Hot variable")

# We will visualize time series of temperature by site and climate profile during the growth season in some figures
# We will also visualize the probability of different temperatures across climate profiles during the growth season
# The following two datasets are obtained as output of the "Growth_season_calculationg.R" script

summarized_climate_data_growing_season <- fread("summarized_climate_data_growing_season.csv")
continuous_proportion_data_growing_season <- fread("climate_profile_temperature_probability_growing_season.csv")

continuous_proportion_data_growing_season_site <- fread("site_temperature_probability_growing_season.csv")

# Here we split the temperature time series by climate profile
cold_stable_climate_time_series_growing_season <- summarized_climate_data_growing_season %>% filter(`Climate profile` == "Cold stable")
hot_stable_climate_time_series_growing_season <- summarized_climate_data_growing_season %>% filter(`Climate profile` == "Hot stable")
cold_variable_climate_time_series_growing_season <- summarized_climate_data_growing_season %>% filter(`Climate profile` == "Cold variable")
hot_variable_climate_time_series_growing_season <- summarized_climate_data_growing_season %>% filter(`Climate profile` == "Hot variable")


# We can find the Phylogenetic tree file obtained from the KBase pipeline and associated metadata in the files below
tpc_tree <- read.tree(file="tpc_tree_kbase_relaxed.tre")

# We then subset the trait metadata to map on our phylogenetic tree tips by metric 
tpc_tree_metadata_auc <- ave_params_gcplyr_combined_common_taxa_auc %>% filter(as.character(seqID) %in% tpc_tree$tip.label & metric == "auc") %>% droplevels()
tpc_tree_metadata_rate <- ave_params_gcplyr_combined_common_taxa_rate %>% filter(as.character(seqID) %in% tpc_tree$tip.label & metric == "rate") %>% droplevels()

# Import and process matrices of phylogenetic distances 
phylogenetic_distance_matrix_auc <- as.matrix(read.table("phylogenetic_distances_auc_kbase.txt"))
phylogenetic_distance_matrix_rate <- as.matrix(read.table("phylogenetic_distances_rate_kbase.txt"))

# We save the tree metadata
write.csv(tpc_tree_metadata_auc, "tpc_tree_kbase_relaxed_metadata_auc.csv",row.names = FALSE)
write.csv(tpc_tree_metadata_rate, "tpc_tree_kbase_relaxed_metadata_rate.csv",row.names = FALSE)

# And we also need the growth curve data used as input to fit the TPCs. This file was generated in the first section of the "TCS_TPC.R" script,
# filtering out low quality TPCs for PCA plotting

tpc_input_gcplyr <- fread("tpc_input_gcplyr_average.csv")

ave_params_gcplyr_combined_auc <- ave_params_gcplyr_combined %>% filter(metric == "auc")
ave_params_gcplyr_combined_rate <- ave_params_gcplyr_combined %>% filter(metric == "rate")
geo_data_auc <- ave_params_gcplyr_combined_auc %>% dplyr::select(c(seqID,`Climate profile`))
geo_data_rate <- ave_params_gcplyr_combined_rate %>% dplyr::select(c(seqID,`Climate profile`))

tpc_input_gcplyr$seqID <- as.character(tpc_input_gcplyr$curve_id)

geo_data_auc$seqID <- as.character(geo_data_auc$seqID)
geo_data_rate$seqID <- as.character(geo_data_rate$seqID)

tpc_input_gcplyr_auc <- tpc_input_gcplyr %>% filter(seqID %in% geo_data_auc$seqID) %>% left_join(geo_data_auc)
tpc_input_gcplyr_rate <- tpc_input_gcplyr %>% filter(seqID %in% geo_data_rate$seqID) %>% left_join(geo_data_rate)

#### Creating the custom color palettes for the main figures ####

# This is the original Spectral palette of the "colornames" package
palette_colors_Spectral <-c("#A71B4B","#DE5925","#F39B29","#FBD476","#A5E8AD","#22C4B3","#0090B5","#584B9F")
# We then obtain a subset of colors for the Orders in our full dataset
#palette_colors_Spectral_Order <-c("#F39B29","#FBD476","#DE5925","#A71B4B","#0090B5","#584B9F")
# We then obtain a subset of colors for the Orders in our full dataset
palette_colors_Spectral_Order <-c("#FBD476","#DE5925","#0090B5","#A71B4B")
# And four colors to indicate the different climate profiles
palette_colors_Spectral_Climate <-c("#584B9F","#22C4B3","#A71B4B","#F39B29")

# And create a color palette for this supplementary figure 
full_sites_palette <- c("darkgreen", "violet", "darkgoldenrod", "grey40",
                        "darkturquoise", "goldenrod2", "black",
                        "darkseagreen4", "turquoise2", "plum2",
                        "darkseagreen2", "turquoise4", "goldenrod4", "grey80", "plum4")

cold_stable_sites_palette <- c("darkgreen", "violet", "darkgoldenrod", "grey40")
cold_variable_sites_palette <- c("darkturquoise", "goldenrod2", "black")
hot_stable_sites_palette <- c("darkseagreen4", "turquoise2", "plum2")
hot_variable_sites_palette <- c("darkseagreen2", "turquoise4", "goldenrod4", "grey80", "plum4")




###### Distribution and variance of growth curves ######
#
#
# Here we continue evaluating the growth curve data across temperatures
# to see the shared constraints to growth variability at the extremes for
# the two main Orders: Enterobacterales and Pseudomonadales
#
#

#### Summarizing the growth curve data by Order, temperature, and metric ####
tpc_input_gcplyr_variance_auc <- tpc_input_gcplyr_auc %>%
  group_by(temp, Order)%>%
  summarize(sd_auc = sd(auc),
            sd_rate = sd(rate),
            auc = mean(auc),
            rate = mean(rate)
  )

tpc_input_gcplyr_variance_rate <- tpc_input_gcplyr_rate %>%
  group_by(temp, Order)%>%
  summarize(sd_auc = sd(auc),
            sd_rate = sd(rate),
            auc = mean(auc),
            rate = mean(rate)
  )


#### S9.top - Plotting the distribution and variance of AUC growth curves ####

SF1C1.auc <- tpc_input_gcplyr_variance_auc %>% 
  filter(Order %in% c("Enterobacterales"))%>%
  ggplot(aes(x = temp,y = sd_auc))+
  stat_smooth(method = "loess", aes(group = Order,fill = Order),color = "white",alpha = 0.3)+
  geom_point(aes(color = Order,size = auc))+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "sd(AUC)"
  ) +
  ylim(0,2.15)+
  scale_color_manual(values = c("#DE5925"))+
  scale_fill_manual(values = c("#DE5925"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")+
  guides(color = "none")+
  guides(fill = "none")


SF1C2.auc <- tpc_input_gcplyr_variance_auc %>% 
  filter(Order %in% c("Pseudomonadales"))%>%
  ggplot(aes(x = temp,y = sd_auc))+
  stat_smooth(method = "loess", aes(group = Order,fill = Order),color = "white",alpha = 0.3)+
  geom_point(aes(color = Order,size = auc))+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "sd(AUC)"
  ) +
  ylim(0,2.15)+
  #scale_color_viridis_d(option = "turbo")+
  scale_color_manual(values = c("#0090B5"))+
  scale_fill_manual(values = c("#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")+
  guides(color = "none")+
  guides(fill = "none")

SF1C3.auc <- tpc_input_gcplyr_auc %>% 
  filter(Order %in% c("Enterobacterales"))%>%
  ggplot(aes(x = temp,y = auc))+
  stat_smooth(method = "loess", aes(group = Order,fill = Order),color = "white",alpha = 0.3)+
  geom_point(aes(color = Order))+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "AUC"
  ) +
  facet_wrap(~Order,ncol = 1)+
  ylim(0,4.5)+
  scale_color_manual(values = c("#DE5925"))+
  scale_fill_manual(values = c("#DE5925"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")+
  annotation_custom(ggplotGrob(SF1C1.auc),xmin = 30, xmax = 44, 
                    ymin = 2.8, ymax = 4.5)


SF1C4.auc <- tpc_input_gcplyr_auc %>% 
  filter(Order %in% c("Pseudomonadales"))%>%
  ggplot(aes(x = temp,y = auc))+
  stat_smooth(method = "loess", aes(group = Order,fill = Order),color = "white",alpha = 0.3)+
  geom_point(aes(color = Order))+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "AUC"
  ) +
  facet_wrap(~Order,ncol = 1)+
  ylim(0,4.5)+
  scale_color_manual(values = c("#0090B5"))+
  scale_fill_manual(values = c("#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        #axis.title.y = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")+
  annotation_custom(ggplotGrob(SF1C2.auc),xmin = 30, xmax = 44, 
                    ymin = 2.8, ymax = 4.5)

# We then arrange the plots into a main figure
S9.top <- plot_grid(SF1C3.auc + theme(legend.position = "none"),SF1C4.auc + theme(legend.position = "none"), ncol = 2, align = "h", rel_widths = c(1, 1),axis = "l")



#### S9.bot - Plotting the distribution and variance of Rate growth curves ####

SF1C1.rate <- tpc_input_gcplyr_variance_rate %>% 
  filter(Order %in% c("Enterobacterales"))%>%
  ggplot(aes(x = temp,y = sd_rate))+
  stat_smooth(method = "loess", aes(group = Order,fill = Order),color = "white",alpha = 0.3)+
  geom_point(aes(color = Order,size = rate))+
  #geom_point(data = tpc_input_gcplyr_variance_general,size = 10,aes(x = temp,y = sd_rate) )+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "sd(rate)"
  ) +
  ylim(0,0.3)+
  scale_color_manual(values = c("#DE5925"))+
  scale_fill_manual(values = c("#DE5925"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")+
  guides(color = "none")+
  guides(fill = "none")


SF1C2.rate <- tpc_input_gcplyr_variance_rate %>% 
  filter(Order %in% c("Pseudomonadales"))%>%
  ggplot(aes(x = temp,y = sd_rate))+
  stat_smooth(method = "loess", aes(group = Order,fill = Order),color = "white",alpha = 0.3)+
  geom_point(aes(color = Order,size = rate))+
  #geom_point(data = tpc_input_gcplyr_variance_general,size = 10,aes(x = temp,y = sd_rate) )+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "sd(rate)"
  ) +
  ylim(0,0.3)+
  #scale_color_viridis_d(option = "turbo")+
  scale_color_manual(values = c("#0090B5"))+
  scale_fill_manual(values = c("#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")+
  guides(color = "none")+
  guides(fill = "none")

SF1C3.rate <- tpc_input_gcplyr_rate %>% 
  filter(Order %in% c("Enterobacterales"))%>%
  ggplot(aes(x = temp,y = rate))+
  stat_smooth(method = "loess", aes(group = Order,fill = Order),color = "white",alpha = 0.3)+
  geom_point(aes(color = Order))+
  #geom_point(data = tpc_input_gcplyr_variance_general,size = 10,aes(x = temp,y = sd_rate) )+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "rate"
  ) +
  facet_wrap(~Order,ncol = 1)+
  ylim(0,1.5)+
  scale_color_manual(values = c("#DE5925"))+
  scale_fill_manual(values = c("#DE5925"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")+
  annotation_custom(ggplotGrob(SF1C1.rate),xmin = 30, xmax = 44, 
                    ymin = 1.0, ymax = 1.5)


SF1C4.rate <- tpc_input_gcplyr_rate %>% 
  filter(Order %in% c("Pseudomonadales"))%>%
  ggplot(aes(x = temp,y = rate))+
  stat_smooth(method = "loess", aes(group = Order,fill = Order),color = "white",alpha = 0.3)+
  geom_point(aes(color = Order))+
  #geom_point(data = tpc_input_gcplyr_variance_general,size = 10,aes(x = temp,y = sd_rate) )+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "rate"
  ) +
  facet_wrap(~Order,ncol = 1)+
  ylim(0,1.5)+
  scale_color_manual(values = c("#0090B5"))+
  scale_fill_manual(values = c("#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #axis.title.x = element_blank(),
        #axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        #axis.title.y = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")+
  annotation_custom(ggplotGrob(SF1C2.rate),xmin = 30, xmax = 44, 
                    ymin = 1, ymax = 1.5)

# We then arrange the plots into a main figure
S9.bot <- plot_grid(SF1C3.rate + theme(legend.position = "none"),SF1C4.rate + theme(legend.position = "none"), ncol = 2, align = "h", rel_widths = c(1, 1),axis = "l")


###### F1A.left - California map ######
#
#
# Here we generate a preliminary figure that shows the location of the sampled sites in a map of California
# We also map basic bioclimatic variables onto it
#
#
# We load a shapefile of California
california_map <- st_read("CA_State.shp")
california_geometry <- st_geometry(california_map)
# We obtain the coordinates of each sampled site from a database
environment_database_noprec <- fread("environment_database_tcs_noprec_250512.csv")
environment_database_full <- fread("environment_database_tcs_full_250512.csv")

reserves_coordinates <- fread("environment_database_tcs_240731.csv")
# We then convert the data to sf object, 
# extract the coordinates from geometry column,
# and combine coordinates with the site names
reserves_coordinates_sf <- st_as_sf(reserves_coordinates, coords = c("lon", "lat"), crs = 4326)
coordinates <- st_coordinates(reserves_coordinates_sf)
data_with_coords <- cbind(reserves_coordinates_sf$Site, coordinates)
# We convert it to a data frame and rename some columns
data_with_coords <- as.data.frame(data_with_coords)
names(data_with_coords) <- c("Site", "lon", "lat")

# We add a column with Site initials to the reserve databases

reserves_coordinates_sf$Site_initials <- factor(reserves_coordinates_sf$Site_initials,
                                                levels = c("B","P","V","Y",
                                                           "SJ","S","W",
                                                           "ME","J","M",
                                                           "AB","DC","QR","PF","EC"))

reserves_coordinates$Site_initials <- factor(reserves_coordinates_sf$Site_initials,
                                             levels = c("B","P","V","Y",
                                                        "SJ","S","W",
                                                        "ME","J","M",
                                                        "AB","DC","QR","PF","EC"))


environment_database_noprec$Site_initials <- factor(environment_database_noprec$Site_initials,
                                             levels = c("B","P","V","Y",
                                                        "SJ","S","W",
                                                        "ME","J","M",
                                                        "AB","DC","QR","PF","EC"))


environment_database_full$Site_initials <- factor(environment_database_full$Site_initials,
                                                    levels = c("B","P","V","Y",
                                                               "SJ","S","W",
                                                               "ME","J","M",
                                                               "AB","DC","QR","PF","EC"))



cold_stable_sites_palette <- c("darkgreen", "violet", "darkgoldenrod", "grey40")
cold_variable_sites_palette <- c("darkturquoise", "goldenrod2", "black")
hot_stable_sites_palette <- c("darkseagreen4", "turquoise2", "plum2")
hot_variable_sites_palette <- c("darkseagreen2", "turquoise4", "goldenrod4", "grey80", "plum4")


# We then plot the map of California with the reserve data, 
# mapping the climate profile as the shape aesthetic

F1A.left <- ggplot() +
  geom_sf(data = california_map, fill = "darkgray") +
  #geom_sf(data = reserves_coordinates_sf, aes(size = 4 , color = Site_initials,shape = `Climate profile`)) +
  geom_sf(data = reserves_coordinates_sf, aes(size = 20 , color = Site_initials)) +
  geom_sf_text(data = reserves_coordinates_sf,aes(label = Site_initials), size = 6,color = "white")+
  #scale_color_scico(palette = "vik", guide = "none")+
  scale_color_manual(values = full_sites_palette)+
  scale_size_continuous(range = c(6, 1), guide = "none") +
  #scale_shape_manual(values = c(15,16,17,18), guide = "none")+# Adjust the size range of points as needed
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    panel.background = element_rect(color = "black"), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'))

###### F1A.right - Dendrograms ######
# 
# 
# This section contains the code to plot the phylogenetic tree obtained from the Kbase pipeline
# as a simplified, ultrametricized, dendrogram version. 
# We use the tpc_tree and tpc_tree_metadata files
#
#


# Plot the original tree
plotTree(tpc_tree,offset=1,type="cladogram")

# We label the nodes to find the one we want to use as root
labelnodes(1:(Ntip(tpc_tree)+tpc_tree$Nnode),
           1:(Ntip(tpc_tree)+tpc_tree$Nnode),
           interactive=FALSE,cex=0.1)




# Reroot the tree at the selected node
#tpc_tree_rooted <- reroot(tpc_tree,node.number = 285)

# Ultrametricize the tree (this is useful to apply Phylogenetic Comparative Methods to bacterial trees,
# and also for plotting purposes, but it is not essential in this case
tpc_tree <- multi2di(tpc_tree)
tpc_tree <- force.ultrametric(tpc_tree)
is.ultrametric(tpc_tree)

# We make the tree and metadata compatible
tpc_tree_metadata_auc$seqID <- as.factor(tpc_tree_metadata_auc$seqID)

rownames(tpc_tree_metadata_auc) <- tpc_tree_metadata_auc$seqID
tpc_tree_metadata_auc <- tpc_tree_metadata_auc[match(tpc_tree$tip.label, row.names(tpc_tree_metadata_auc)), ]


F1A.right <- ggtree(tpc_tree, layout = "circular", branch.length = "none") %<+% tpc_tree_metadata_auc +
  #geom_tiplab(size = 1.5, # color for label font
  #            geom = "label",  # labels not text
  #            offset = 0.3,
  #            label.padding = unit(0, "lines"), # amount of padding around the labels
  #            label.size = 0) + # size of label border
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(legend.position = "none")+ #+
  geom_tippoint(size = 4, aes(color = Site_initials))+
  #geom_tiplab(offset = 0.5, size = 1)+
  #scale_shape_manual(values = c(16,15,18,17))+
  #scale_color_viridis_d(option = "inferno", guide = "none")+
  scale_color_manual(values = full_sites_palette)+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(legend.position = "none", 
        legend.background = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank())


# Add unique Order labels to the outer ring
F1A.right <- F1A.right +
  #geom_cladelabel(node = 275, label = "Weeksellaceae", offset = 6, offset.text = 2, align = FALSE, barsize = 0.5, fontsize = 3.88, angle = 0, geom = "text", hjust = 0, fill = NA, family = "sans")+
  #geom_cladelabel(node = 276, label = "Bacillaceae_H", offset = 6, offset.text = 2, align = FALSE, barsize = 0.5, fontsize = 3.88, angle = 0, geom = "text", hjust = 0, fill = NA, family = "sans")+
  #geom_cladelabel(node = 281, label = "Xanthomonadaceae", offset = 6, offset.text = 2, align = FALSE, barsize = 0.5, fontsize = 3.88, angle = 0, geom = "text", hjust = 0, fill = NA, family = "sans")+
  #geom_cladelabel(node = 286, label = "Burkholderiaceae", offset = 6, offset.text = 2, align = FALSE, barsize = 0.5, fontsize = 3.88, angle = 0, geom = "text", hjust = 0, fill = NA, family = "sans")+
  #geom_cladelabel(node = 292, label = "Moraxellaceae", offset = 6, offset.text = 2, align = FALSE, barsize = 0.5, fontsize = 3.88, angle = 0, geom = "text", hjust = 0, fill = NA, family = "sans")+
  #geom_cladelabel(node = 295, label = "Enterobacteriaceae", offset = 0, offset.text = 0, align = FALSE, barsize = 0.5, fontsize = 3.88, angle = 0, geom = "text", hjust = 0, fill = NA, family = "sans")+
  #geom_cladelabel(node = 295, label = "Enterobacteriaceae", offset = 6, offset.text = 2, align = FALSE, barsize = 0.5, fontsize = 3.88, angle = 0, geom = "text", hjust = 0, fill = NA, family = "sans")+
  #geom_cladelabel(node = 323, label = "Pseudomonadaceae", offset = 0, offset.text = 0, align = FALSE, barsize = 0.5, fontsize = 3.88, angle = 0, geom = "text", hjust = 0, fill = NA, family = "sans")+
  #geom_cladelabel(node = 350, label = "Pseudomonadaceae", offset = 6, offset.text = 20, align = FALSE, barsize = 0.5, fontsize = 3.88, angle = 0, geom = "text", hjust = 0, fill = NA, family = "sans")+
  geom_strip(4, 5, barsize=1, offset = 1.5, label = "Pseudomonadales",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(2, 4, barsize=1, offset = 1.5, label = "Burkholderiales",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(5, 45, barsize=1, offset = 1.5, label = "Enterobacterales",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(56, 265, barsize=1, offset = 1.5, label = "",offset.text = 1, align = TRUE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT")+
  geom_strip(1, 2, barsize=1, offset = 1.5, label = "Xanthomonadales",offset.text = 1, align = TRUE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT")+
  geom_strip(46, 84, barsize=1, offset = 1.5, label = "Pseudomonadales",offset.text = 2, align = TRUE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT")



###### F2A - Visual comparison of thermal traits across Orders ######
#
#
# In this code section we combine the TPC fits and parameters at the Order level
# colored by climate profile for a subset of relatively common clades (Genera >1%)
#
#
#### F1C AUC and S6 ####

# This is a figure displaying the TPCs of the genera present in our collection,
# colored by site

S6 <- ave_preds_gcplyr_combined_major_taxa_auc %>% 
  ggplot( aes(x = temp, y = .fitted)) +
  geom_line(lwd = 0.75,alpha = 0.8,aes( group = seqID,color = Site_initials)) +
  facet_wrap(~factor(Genus),ncol=8,strip.position = "top")+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)+
  xlab("Temperature [°C]")+
  xlim(13,45)+
  ylim(0,7)+
  ylab(expression(AUC))+
  #scale_color_viridis_d()+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        #axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        #axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "bottom"
  )

# This is a figure displaying the TPCs of the two main Orders,
# colored by site

F1C1.auc <- ave_preds_gcplyr_combined_major_taxa_auc %>% 
  ggplot( aes(x = temp, y = .fitted)) +
  geom_line(lwd = 0.75,alpha = 0.8,aes( group = seqID,color = Site_initials)) +
  geom_ribbon(data = ave_preds_gcplyr_combined_major_taxa_auc %>% 
                group_by(Order,temp) %>%
                summarize(sd_fitted = sd(.fitted),
                          .fitted = mean(.fitted)),
              aes(x = temp, ymin = .fitted - sd_fitted, 
                  ymax = .fitted + sd_fitted ),lwd = 2,fill = "grey40",alpha = 0.4)+
  geom_line(data = ave_preds_gcplyr_combined_major_taxa_auc %>% 
              group_by(Order,temp) %>%
              summarize(.fitted = mean(.fitted)),
            aes(x = temp, y = .fitted),lwd = 1.3,color = "grey95")+
  facet_wrap(~factor(Order),ncol=2)+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)+
  xlim(13,45)+
  ylim(0,7)+
  ylab(expression(AUC))+
  #scale_color_viridis_d()+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        #axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none"
  )

F1C2.auc <- ave_params_gcplyr_combined_major_taxa_auc %>% 
  ggplot(aes(x=Order, y=rmax))+
  #scale_color_viridis_b()+
  geom_jitter(width = 0.2,size=4,aes(color=Site_initials))+
  geom_boxplot(alpha = 0.8,outlier.shape = NA)+
  #facet_wrap(~CS)+
  geom_smooth(method = "lm", se = TRUE) +
  #ylim(0,5)+
  #xlim(0,2)+
  #scale_fill_manual(values = order_color)+
  #xlab("Family")+
  ylab(expression(max[AUC]))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none")+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)

F1C3.auc <- ave_params_gcplyr_combined_major_taxa_auc %>%
  ggplot(aes(x=Order, y=topt))+
  #scale_color_viridis_b()+
  geom_jitter(width = 0.2,size=4, aes(color=Site_initials))+
  geom_boxplot(alpha = 0.8,outlier.shape = NA)+
  #facet_wrap(~CS)+
  geom_smooth(method = "lm", se = TRUE) +
  #scale_y_continuous(breaks = 5)+
  #ylim(0,5)+
  #xlim(0,2)+
  #scale_fill_manual(values = order_color)+
  #xlab("Family")+
  scale_y_continuous(labels = scales::label_number(accuracy = 1))+
  ylab(expression(topt[AUC]~(degree*C)))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none")+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)


F1C4.auc <- ave_params_gcplyr_combined_major_taxa_auc %>%
  ggplot(aes(x=Order, y= ctmax))+
  #scale_color_viridis_b()+
  geom_jitter(width = 0.2,size=4, aes(color=Site_initials))+
  geom_boxplot(alpha = 0.8,outlier.shape = NA)+
  #facet_wrap(~CS)+
  geom_smooth(method = "lm", se = TRUE) +
  #ylim(0,5)+
  #xlim(0,2)+
  #scale_fill_manual(values = order_color)+
  #xlab("Family")+
  ylab(expression(CTmax[AUC]~(degree*C)))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 0),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none")+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)

F1C5.auc <- ave_params_gcplyr_combined_major_taxa_auc %>%
  ggplot(aes(x=Order, y= ctmin))+
  #scale_color_viridis_b()+
  geom_jitter(width = 0.2,size=4, aes(color=Site_initials))+
  geom_boxplot(alpha = 0.8,outlier.shape = NA)+
  #facet_wrap(~CS)+
  geom_smooth(method = "lm", se = TRUE) +
  #ylim(0,5)+
  #xlim(0,2)+
  #scale_fill_manual(values = order_color)+
  #xlab("Family")+
  ylab(expression(CTmin[AUC]~(degree*C)))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 0),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none")+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)

# Arrange F2 #

# Combine plots
F1C.auc <- plot_grid(F1C1.auc, F1C2.auc,F1C3.auc,F1C4.auc,F1C5.auc, ncol = 1, align = "v", rel_heights = c(1.5, 1.3,1.3,1.3,1.4),axis = "l")
# Show combined plot
F1C.auc

#### S14B.left ####

F1C1.rate <- ave_preds_gcplyr_combined_major_taxa_rate %>% 
  ggplot( aes(x = temp, y = .fitted)) +
  geom_line(lwd = 0.75,alpha = 0.6,aes( group = seqID,color = Site_initials)) +
  geom_ribbon(data = ave_preds_gcplyr_combined_major_taxa_rate %>% 
                group_by(Order,temp) %>%
                summarize(sd_fitted = sd(.fitted),
                          .fitted = mean(.fitted)),
              aes(x = temp, ymin = .fitted - sd_fitted, 
                  ymax = .fitted + sd_fitted ),lwd = 2,fill = "grey40",alpha = 0.4)+
  geom_line(data = ave_preds_gcplyr_combined_major_taxa_rate %>% 
              group_by(Order,temp) %>%
              summarize(.fitted = mean(.fitted)),
            aes(x = temp, y = .fitted),lwd = 1.3,color = "grey95")+
  facet_wrap(~factor(Order),ncol=6)+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)+
  xlim(13,45)+
  ylim(0,1)+
  ylab(expression(rate))+
  #scale_color_viridis_d()+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none"
  )

F1C2.rate <- ave_params_gcplyr_combined_major_taxa_rate %>% 
  ggplot(aes(x=Order, y=rmax))+
  #scale_color_viridis_b()+
  geom_jitter(width = 0.2,size=4, aes(color=Site_initials))+
  geom_boxplot(alpha = 0.8,outlier.shape = NA)+
  #facet_wrap(~CS)+
  geom_smooth(method = "lm", se = TRUE) +
  #ylim(0,5)+
  #xlim(0,2)+
  #scale_fill_manual(values = order_color)+
  #xlab("Family")+
  ylab(expression(max[rate]))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none")+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)

F1C3.rate <- ave_params_gcplyr_combined_major_taxa_rate %>%
  ggplot(aes(x=Order, y=topt))+
  #scale_color_viridis_b()+
  geom_jitter(width = 0.2,size=4,aes(color=Site_initials))+
  geom_boxplot(alpha = 0.8,outlier.shape = NA)+
  #facet_wrap(~CS)+
  geom_smooth(method = "lm", se = TRUE) +
  #scale_y_continuous(breaks = 5)+
  #ylim(0,5)+
  #xlim(0,2)+
  #scale_fill_manual(values = order_color)+
  #xlab("Family")+
  scale_y_continuous(labels = scales::label_number(accuracy = 1))+
  ylab(expression(topt[rate]~(degree*C)))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none")+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)


F1C4.rate <- ave_params_gcplyr_combined_major_taxa_rate %>%
  ggplot(aes(x=Order, y= ctmax))+
  #scale_color_viridis_b()+
  geom_jitter(width = 0.2,size=4, aes(color=Site_initials))+
  geom_boxplot(alpha = 0.8,outlier.shape = NA)+
  #facet_wrap(~CS)+
  geom_smooth(method = "lm", se = TRUE) +
  #ylim(0,5)+
  #xlim(0,2)+
  #scale_fill_manual(values = order_color)+
  #xlab("Family")+
  ylab(expression(CTmax[rate]~(degree*C)))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 0),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none")+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)

F1C5.rate <- ave_params_gcplyr_combined_major_taxa_rate %>%
  ggplot(aes(x=Order, y= ctmin))+
  #scale_color_viridis_b()+
  geom_jitter(width = 0.2,size=4,aes(color=Site_initials))+
  geom_boxplot(alpha = 0.8,outlier.shape = NA)+
  #facet_wrap(~CS)+
  geom_smooth(method = "lm", se = TRUE) +
  #ylim(0,5)+
  #xlim(0,2)+
  #scale_fill_manual(values = order_color)+
  #xlab("Family")+
  ylab(expression(CTmin[rate]~(degree*C)))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 0),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none")+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)

# Arrange F2 #

# Combine plots
S14B.left <- plot_grid(F1C1.rate, F1C2.rate,F1C3.rate,F1C4.rate,F1C5.rate, ncol = 1, align = "v", rel_heights = c(1.5, 1.3,1.3,1.3,1.4),axis = "l")


#### Basic statistics for F1C ####
#
#
# Here we run basic one-way ANOVAs to compare thermal traits between the plotted Orders,
# Climate profiles, and post-hoc pairwise comparisons using Tukey 
#
#

# Stats for F1C.auc by Order

anova_model_Order_aucmax <- aov(rmax ~ Order, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_Order_aucmax)
posthoc_Order_aucmax <- TukeyHSD(anova_model_Order_aucmax)
posthoc_Order_aucmax

anova_model_Order_topt_auc <- aov(topt ~ Order, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_Order_topt_auc)
posthoc_Order_topt_auc <- TukeyHSD(anova_model_Order_topt_auc)
posthoc_Order_topt_auc

anova_model_Order_ctmax_auc <- aov(ctmax ~ Order, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_Order_ctmax_auc)
posthoc_Order_ctmax_auc <- TukeyHSD(anova_model_Order_ctmax_auc)
posthoc_Order_ctmax_auc

anova_model_Order_ctmin_auc <- aov(ctmin ~ Order, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_Order_ctmin_auc)
posthoc_Order_ctmin_auc <- TukeyHSD(anova_model_Order_ctmin_auc)
posthoc_Order_ctmin_auc

# Stats for F1C.auc by Climate Profile

# The TukeyHSD function cannot handle variables enclosed by ``, so we create a new climate_profile variable
ave_params_gcplyr_combined_major_taxa_auc$climate_profile <- as.factor(ave_params_gcplyr_combined_major_taxa_auc$`Climate profile`)

anova_model_tprofile_aucmax <- aov(rmax ~ climate_profile, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_tprofile_aucmax)
posthoc_tprofile_aucmax <- TukeyHSD(anova_model_tprofile_aucmax)
posthoc_tprofile_aucmax

anova_model_tprofile_topt_auc <- aov(topt ~ climate_profile, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_tprofile_topt_auc)
posthoc_tprofile_topt_auc <- TukeyHSD(anova_model_tprofile_topt_auc)
posthoc_tprofile_topt_auc

anova_model_tprofile_ctmax_auc <- aov(ctmax ~ climate_profile, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_tprofile_ctmax_auc)
posthoc_tprofile_ctmax_auc <- TukeyHSD(anova_model_tprofile_ctmax_auc)
posthoc_tprofile_ctmax_auc

anova_model_tprofile_ctmin_auc <- aov(ctmin ~ climate_profile, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_tprofile_ctmin_auc)
posthoc_tprofile_ctmin_auc <- TukeyHSD(anova_model_tprofile_ctmin_auc)
posthoc_tprofile_ctmin_auc

# We then run a 2-way ANOVA to evaluate the relative contribution of Order and Climate Profile 
# in explaining the variation in thermal traits 
anova_model_Order_tprofile_aucmax <- aov(rmax ~ Order*climate_profile, data = ave_params_gcplyr_combined_major_taxa_auc %>% 
                                     filter(Order %in% c("Enterobacterales","Pseudomonadales")))
emmeans(anova_model_Order_tprofile_aucmax, pairwise ~ Order * climate_profile)

anova_model_Order_tprofile_topt_auc <- aov(topt ~ Order*climate_profile, data = ave_params_gcplyr_combined_major_taxa_auc %>% 
                                     filter(Order %in% c("Enterobacterales","Pseudomonadales")))
emmeans(anova_model_Order_tprofile_topt_auc, pairwise ~ Order * climate_profile)

anova_model_Order_tprofile_ctmax_auc <- aov(ctmax ~ Order*climate_profile, data = ave_params_gcplyr_combined_major_taxa_auc %>% 
                                       filter(Order %in% c("Enterobacterales","Pseudomonadales")))
emmeans(anova_model_Order_tprofile_ctmax_auc, pairwise ~ Order * climate_profile)

anova_model_Order_tprofile_ctmin_auc <- aov(ctmin ~ Order*climate_profile, data = ave_params_gcplyr_combined_major_taxa_auc %>% 
                                        filter(Order %in% c("Enterobacterales","Pseudomonadales")))
emmeans(anova_model_Order_tprofile_ctmin_auc, pairwise ~ Order * climate_profile)

# Stats for F2A.rate by Order

anova_model_Order_ratemax <- aov(rmax ~ Order, data = ave_params_gcplyr_combined_common_taxa_rate)
summary(anova_model_Order_ratemax)
posthoc_Order_ratemax <- TukeyHSD(anova_model_Order_ratemax)
posthoc_Order_ratemax

anova_model_Order_topt_rate <- aov(topt ~ Order, data = ave_params_gcplyr_combined_common_taxa_rate)
summary(anova_model_Order_topt_rate)
posthoc_Order_topt_rate <- TukeyHSD(anova_model_Order_topt_rate)
posthoc_Order_topt_rate

anova_model_Order_ctmax_rate <- aov(ctmax ~ Order, data = ave_params_gcplyr_combined_common_taxa_rate)
summary(anova_model_Order_ctmax_rate)
posthoc_Order_ctmax_rate <- TukeyHSD(anova_model_Order_ctmax_rate)
posthoc_Order_ctmax_rate

anova_model_Order_ctmin_rate <- aov(ctmin ~ Order, data = ave_params_gcplyr_combined_common_taxa_rate)
summary(anova_model_Order_ctmin_rate)
posthoc_Order_ctmin_rate <- TukeyHSD(anova_model_Order_ctmin_rate)
posthoc_Order_ctmin_rate

# Stats for F1C.rate by Climate Profile

# The TukeyHSD function cannot handle variables enclosed by ``, so we create a new climate_profile variable
ave_params_gcplyr_combined_common_taxa_rate$climate_profile <- as.factor(ave_params_gcplyr_combined_common_taxa_rate$`Climate profile`)

anova_model_tprofile_ratemax <- aov(rmax ~ climate_profile, data = ave_params_gcplyr_combined_common_taxa_rate)
summary(anova_model_tprofile_ratemax)
posthoc_tprofile_ratemax <- TukeyHSD(anova_model_tprofile_ratemax)
posthoc_tprofile_ratemax

anova_model_tprofile_topt_rate <- aov(topt ~ climate_profile, data = ave_params_gcplyr_combined_common_taxa_rate)
summary(anova_model_tprofile_topt_rate)
posthoc_tprofile_topt_rate <- TukeyHSD(anova_model_tprofile_topt_rate)
posthoc_tprofile_topt_rate

anova_model_tprofile_ctmax_rate <- aov(ctmax ~ climate_profile, data = ave_params_gcplyr_combined_common_taxa_rate)
summary(anova_model_tprofile_ctmax_rate)
posthoc_tprofile_ctmax_rate <- TukeyHSD(anova_model_tprofile_ctmax_rate)
posthoc_tprofile_ctmax_rate

anova_model_tprofile_ctmin_rate <- aov(ctmin ~ climate_profile, data = ave_params_gcplyr_combined_common_taxa_rate)
summary(anova_model_tprofile_ctmin_rate)
posthoc_tprofile_ctmin_rate <- TukeyHSD(anova_model_tprofile_ctmin_rate)
posthoc_tprofile_ctmin_rate

# We then run a 2-way ANOVA to evaluate the relative contribution of Order and Climate Profile 
# in explaining the variation in thermal traits 
anova_model_Order_tprofile_ratemax <- aov(rmax ~ Order*climate_profile, data = ave_params_gcplyr_combined_common_taxa_rate %>% 
                                           filter(Order %in% c("Enterobacterales","Pseudomonadales")))
emmeans(anova_model_Order_tprofile_ratemax, pairwise ~ Order * climate_profile)

anova_model_Order_tprofile_topt_rate <- aov(topt ~ Order*climate_profile, data = ave_params_gcplyr_combined_common_taxa_rate %>% 
                                             filter(Order %in% c("Enterobacterales","Pseudomonadales")))
emmeans(anova_model_Order_tprofile_topt_rate, pairwise ~ Order * climate_profile)

anova_model_Order_tprofile_ctmax_rate <- aov(ctmax ~ Order*climate_profile, data = ave_params_gcplyr_combined_common_taxa_rate %>% 
                                              filter(Order %in% c("Enterobacterales","Pseudomonadales")))
emmeans(anova_model_Order_tprofile_ctmax_rate, pairwise ~ Order * climate_profile)

anova_model_Order_tprofile_ctmin_rate <- aov(ctmin ~ Order*climate_profile, data = ave_params_gcplyr_combined_common_taxa_rate %>% 
                                              filter(Order %in% c("Enterobacterales","Pseudomonadales")))
emmeans(anova_model_Order_tprofile_ctmin_rate, pairwise ~ Order * climate_profile)
###### F1D.auc - Calculating Phylogenetic Signal for major taxa (>1% of collection) AUC ######

#
#
# Here we apply the methods detailed in Lennon et al (2016) to calculate the phylogenetic
# signal of bacterial traits to our thermal traits, by fitting hierarchical linear models
# for each trait and metric.
#
#


n_permutations <- 10000

phylogenetic_signal_mixed_model_f <- function(trait_name, data) {
  
  # Fit hierarchical linear model with Satterthwaite df
  formula <- as.formula(paste(trait_name, "~ PC1 + PC2 + temperature_of_isolation + (1|Order/Genus)"))
  model <- lmer(formula, data = data)
  
  # Partial marginal R² partitioning
  part_model <- partR2(model, partvars = c("PC1","PC2","temperature_of_isolation"), data = data)
  
  # Extract fixed effect table (now includes p-values via lmerTest)
  fixed_summary <- as.data.frame(summary(model)$coefficients)
  fixed_summary$grp <- rownames(fixed_summary)
  rownames(fixed_summary) <- NULL
  
  # Keep only fixed effects present in part_model (drop intercept)
  fixed_var <- fixed_summary[fixed_summary$grp %in% part_model$R2$term, ]
  
  # Add variance share by matching term names
  fixed_var$share <- part_model$R2$estimate[match(fixed_var$grp, part_model$R2$term)]
  
  # Random effects / residual variance: phylogenetic signal
  observed_var <- as.data.frame(VarCorr(model))
  observed_total <- sum(observed_var$vcov)
  observed_var$proportion <- observed_var$vcov / observed_total
  
  # Run permutation test on phylogenetic signal
  permuted_proportions <- matrix(NA, nrow = n_permutations, ncol = nrow(observed_var))
  colnames(permuted_proportions) <- observed_var$grp
  
  for(i in 1:n_permutations){
    data[[trait_name]] <- sample(data[[trait_name]]) # permute trait
    perm_model <- lmer(formula, data = data)         
    perm_var <- as.data.frame(VarCorr(perm_model))   
    permuted_proportions[i, ] <- perm_var$vcov / sum(perm_var$vcov)
  }
  
  ci_lower <- apply(permuted_proportions, 2, quantile, probs = 0.025)
  ci_upper <- apply(permuted_proportions, 2, quantile, probs = 0.975)
  
  confint_df <- data.frame(
    observed = observed_var$proportion,
    lower_ci = ci_lower,
    upper_ci = ci_upper,
    `Taxonomic Level` = factor(c("Genus","Order","Residual")),
    Trait = trait_name,
    significant = observed_var$proportion < ci_lower | observed_var$proportion > ci_upper
  )
  
  list(fixed_var = fixed_var, phylo_signal = confint_df)
}

# Apply to all traits
results_aucmax <- phylogenetic_signal_mixed_model_f("rmax", ave_params_gcplyr_combined_major_taxa_auc)
results_topt   <- phylogenetic_signal_mixed_model_f("topt", ave_params_gcplyr_combined_major_taxa_auc)
results_ctmax  <- phylogenetic_signal_mixed_model_f("ctmax", ave_params_gcplyr_combined_major_taxa_auc)
results_ctmin  <- phylogenetic_signal_mixed_model_f("ctmin", ave_params_gcplyr_combined_major_taxa_auc)

# Combine phylogenetic signal data
phylogenetic_signal_thermal_traits_auc <- bind_rows(
  results_aucmax$phylo_signal,
  results_topt$phylo_signal,
  results_ctmax$phylo_signal,
  results_ctmin$phylo_signal
)

# Combine fixed effect contributions (keep p-values)
fixed_effects_all_auc <- bind_rows(
  cbind(Trait="rmax", results_aucmax$fixed_var),
  cbind(Trait="topt", results_topt$fixed_var),
  cbind(Trait="ctmax", results_ctmax$fixed_var),
  cbind(Trait="ctmin", results_ctmin$fixed_var)
)

# Clean + reshape for final table
random_effects_clean_auc <- phylogenetic_signal_thermal_traits_auc %>%
  mutate(
    metric = str_extract(Trait, "(?<=_)\\w+$"),
    Trait = str_remove(Trait, "_auc$"),
    Variable = `Taxonomic.Level`,
    Effect = "random"
  ) %>%
  rename(
    Variance_observed = observed,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    Significant = significant
  ) %>%
  dplyr::select(Trait, metric, Variable, Effect, Variance_observed, Lower_CI, Upper_CI, Significant)

fixed_effects_clean_auc <- fixed_effects_all_auc %>%
  mutate(
    metric = "auc",
    Trait = Trait,
    Variable = grp,
    Effect = "fixed",
    p_value = `Pr(>|t|)`   # Satterthwaite p-value approximation
  ) %>%
  rename(
    Estimate = Estimate,
    Std_Error = `Std. Error`,
    t_value = `t value`,
    Partial_R2 = share
  ) %>%
  dplyr::select(Trait, metric, Variable, Effect, Estimate, Std_Error, t_value, p_value, Partial_R2)

# Final combined table
combined_effects_lmer_auc <- bind_rows(fixed_effects_clean_auc, random_effects_clean_auc)

phylogenetic_signal_thermal_traits_auc$Trait <- paste0(
  phylogenetic_signal_thermal_traits_auc$Trait, "_auc"
)

fixed_effects_all_auc$Trait <- paste0(
  fixed_effects_all_auc$Trait, "_auc"
)

combined_effects_lmer_auc <- combined_effects_lmer_auc %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

F2B.auc <- fixed_effects_all_auc %>%
  ggplot(aes(x = `t value`, y = grp,size = share*100))+
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "grey40")+
  geom_point()+
  geom_segment(xend = 0,linewidth = 0.01,linetype = "dashed")+
  #ylab("PC2 -------> Temperature Mean and Minima")+
  xlab("t-score")+
  #xlim(-1,1)+
  #ylim(-1,1)+
  #coord_flip()+
  #scale_color_manual(values = palette_colors_Spectral_Order)+
  facet_grid(Trait~.)+
  #scale_shape_manual(values = c(15,16,17,18))+
  scale_size_continuous(range = c(1,5))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.direction = "horizontal",
        legend.text = element_text(angle = 0),
        legend.position = "top",
        axis.title.y = element_blank())+
  guides(size = guide_legend("% Variance Explained"))
  


write.csv(combined_effects_lmer_auc, "combined_effects_lmer_auc.csv", row.names = FALSE)
write.csv(phylogenetic_signal_thermal_traits_auc,"phylogenetic_signal_thermal_traits_auc.csv",row.names = FALSE)
write.csv(fixed_effects_all_auc,"fixed_effects_lmer_auc.csv",row.names = FALSE)


combined_effects_lmer_auc <- fread("combined_effects_lmer_auc.csv")

phylogenetic_signal_thermal_traits_auc <- fread("phylogenetic_signal_thermal_traits_auc.csv")



phylogenetic_signal_thermal_traits_auc$Trait <- factor(phylogenetic_signal_thermal_traits_auc$Trait, 
                                                       levels = c("rmax_auc","topt_auc","ctmax_auc","ctmin_auc"))


phylogenetic_signal_thermal_traits_auc$Taxonomic.Level <- factor(phylogenetic_signal_thermal_traits_auc$Taxonomic.Level, 
                                                                 levels = c("Order","Genus","Residual"))


F1D.auc <- phylogenetic_signal_thermal_traits_auc %>% 
  #filter(!`Taxonomic Level` == "Family") %>% 
  mutate(Trait = fct_relevel(Trait,"ctmin_auc","ctmax_auc","topt_auc","rmax_auc")) %>%
  #mutate(Trait = fct_relevel(Trait,"ctmin_auc","ctmax_auc","topt_auc","aucmax" )) %>%
  ggplot(aes(x = Trait, y = observed)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), linewidth = 0.5,width = 0.2, color = "gray50") +
  geom_point(aes(color = significant), size = 6) +
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  ylim(0,1)+
  scale_y_continuous(breaks = seq(0, 1, 0.5))+
  coord_flip()+
  facet_wrap( ~Taxonomic.Level,ncol = 3)+
  labs(
    y = "Proportion of Variance Explained",
    color = "Significant"
  ) +
  scale_color_manual(values = c("grey60","black"))+
  theme(strip.background = element_blank(),
        axis.line.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        #axis.title.x = element_blank(),
        #strip.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "top",
        legend.direction = "horizontal"
  )


###### S14B.right - Calculating Phylogenetic Signal for major taxa (>1% of collection) Rate ######

#
#
# Here we apply the methods detailed in Lennon et al (2016) to calculate the phylogenetic
# signal of bacterial traits to our thermal traits, by fitting hierarchical linear models
# for each trait and metric.
#
#


n_permutations <- 10000

phylogenetic_signal_mixed_model_f <- function(trait_name, data) {
  
  # Fit hierarchical linear model with Satterthwaite df
  formula <- as.formula(paste(trait_name, "~ PC1 + PC2 + temperature_of_isolation + (1|Order/Genus)"))
  model <- lmer(formula, data = data)
  
  # Partial marginal R² partitioning
  part_model <- partR2(model, partvars = c("PC1","PC2","temperature_of_isolation"), data = data)
  
  # Extract fixed effect table (now includes p-values via lmerTest)
  fixed_summary <- as.data.frame(summary(model)$coefficients)
  fixed_summary$grp <- rownames(fixed_summary)
  rownames(fixed_summary) <- NULL
  
  # Keep only fixed effects present in part_model (drop intercept)
  fixed_var <- fixed_summary[fixed_summary$grp %in% part_model$R2$term, ]
  
  # Add variance share by matching term names
  fixed_var$share <- part_model$R2$estimate[match(fixed_var$grp, part_model$R2$term)]
  
  # Random effects / residual variance: phylogenetic signal
  observed_var <- as.data.frame(VarCorr(model))
  observed_total <- sum(observed_var$vcov)
  observed_var$proportion <- observed_var$vcov / observed_total
  
  # Run permutation test on phylogenetic signal
  permuted_proportions <- matrix(NA, nrow = n_permutations, ncol = nrow(observed_var))
  colnames(permuted_proportions) <- observed_var$grp
  
  for(i in 1:n_permutations){
    data[[trait_name]] <- sample(data[[trait_name]]) # permute trait
    perm_model <- lmer(formula, data = data)         
    perm_var <- as.data.frame(VarCorr(perm_model))   
    permuted_proportions[i, ] <- perm_var$vcov / sum(perm_var$vcov)
  }
  
  ci_lower <- apply(permuted_proportions, 2, quantile, probs = 0.025)
  ci_upper <- apply(permuted_proportions, 2, quantile, probs = 0.975)
  
  confint_df <- data.frame(
    observed = observed_var$proportion,
    lower_ci = ci_lower,
    upper_ci = ci_upper,
    `Taxonomic Level` = factor(c("Genus","Order","Residual")),
    Trait = trait_name,
    significant = observed_var$proportion < ci_lower | observed_var$proportion > ci_upper
  )
  
  list(fixed_var = fixed_var, phylo_signal = confint_df)
}

# Apply to all traits
results_ratemax <- phylogenetic_signal_mixed_model_f("rmax", ave_params_gcplyr_combined_major_taxa_rate)
results_topt   <- phylogenetic_signal_mixed_model_f("topt", ave_params_gcplyr_combined_major_taxa_rate)
results_ctmax  <- phylogenetic_signal_mixed_model_f("ctmax", ave_params_gcplyr_combined_major_taxa_rate)
results_ctmin  <- phylogenetic_signal_mixed_model_f("ctmin", ave_params_gcplyr_combined_major_taxa_rate)

# Combine phylogenetic signal data
phylogenetic_signal_thermal_traits_rate <- bind_rows(
  results_ratemax$phylo_signal,
  results_topt$phylo_signal,
  results_ctmax$phylo_signal,
  results_ctmin$phylo_signal
)

# Combine fixed effect contributions (keep p-values)
fixed_effects_all_rate <- bind_rows(
  cbind(Trait="rmax", results_ratemax$fixed_var),
  cbind(Trait="topt", results_topt$fixed_var),
  cbind(Trait="ctmax", results_ctmax$fixed_var),
  cbind(Trait="ctmin", results_ctmin$fixed_var)
)

# Clean + reshape for final table
random_effects_clean_rate <- phylogenetic_signal_thermal_traits_rate %>%
  mutate(
    metric = str_extract(Trait, "(?<=_)\\w+$"),
    Trait = str_remove(Trait, "_rate$"),
    Variable = `Taxonomic.Level`,
    Effect = "random"
  ) %>%
  rename(
    Variance_observed = observed,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    Significant = significant
  ) %>%
  dplyr::select(Trait, metric, Variable, Effect, Variance_observed, Lower_CI, Upper_CI, Significant)

fixed_effects_clean_rate <- fixed_effects_all_rate %>%
  mutate(
    metric = "rate",
    Trait = Trait,
    Variable = grp,
    Effect = "fixed",
    p_value = `Pr(>|t|)`   # Satterthwaite p-value approximation
  ) %>%
  rename(
    Estimate = Estimate,
    Std_Error = `Std. Error`,
    t_value = `t value`,
    Partial_R2 = share
  ) %>%
  dplyr::select(Trait, metric, Variable, Effect, Estimate, Std_Error, t_value, p_value, Partial_R2)

# Final combined table
combined_effects_lmer_rate <- bind_rows(fixed_effects_clean_rate, random_effects_clean_rate)

phylogenetic_signal_thermal_traits_rate$Trait <- paste0(
  phylogenetic_signal_thermal_traits_rate$Trait, "_rate"
)

fixed_effects_all_rate$Trait <- paste0(
  fixed_effects_all_rate$Trait, "_rate"
)

combined_effects_lmer_rate <- combined_effects_lmer_rate %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

S15.bot <- fixed_effects_all_rate %>%
  ggplot(aes(x = `t value`, y = grp,size = share*100))+
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "grey40")+
  geom_point()+
  geom_segment(xend = 0,linewidth = 0.01,linetype = "dashed")+
  #ylab("PC2 -------> Temperature Mean and Minima")+
  xlab("t-score")+
  #xlim(-1,1)+
  #ylim(-1,1)+
  #coord_flip()+
  #scale_color_manual(values = palette_colors_Spectral_Order)+
  facet_grid(Trait~.)+
  #scale_shape_manual(values = c(15,16,17,18))+
  scale_size_continuous(range = c(1,5))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.direction = "horizontal",
        legend.text = element_text(angle = 0),
        legend.position = "top",
        axis.title.y = element_blank())+
  guides(size = guide_legend("% Variance Explained"))


write.csv(combined_effects_lmer_rate, "combined_effects_lmer_rate.csv", row.names = FALSE)
write.csv(phylogenetic_signal_thermal_traits_rate,"phylogenetic_signal_thermal_traits_rate.csv",row.names = FALSE)
write.csv(fixed_effects_all_rate,"fixed_effects_lmer_rate.csv",row.names = FALSE)


phylogenetic_signal_thermal_traits_rate <- fread("phylogenetic_signal_thermal_traits_rate.csv")



phylogenetic_signal_thermal_traits_rate$Trait <- factor(phylogenetic_signal_thermal_traits_rate$Trait, 
                                                       levels = c("rmax_rate","topt_rate","ctmax_rate","ctmin_rate"))


phylogenetic_signal_thermal_traits_rate$Taxonomic.Level <- factor(phylogenetic_signal_thermal_traits_rate$Taxonomic.Level, 
                                                                 levels = c("Order","Genus","Residual"))


S14B.right <- phylogenetic_signal_thermal_traits_rate %>% 
  #filter(!`Taxonomic Level` == "Family") %>% 
  mutate(Trait = fct_relevel(Trait,"ctmin_rate","ctmax_rate","topt_rate","rmax_rate")) %>%
  #mutate(Trait = fct_relevel(Trait,"ctmin_rate","ctmax_rate","topt_rate","ratemax" )) %>%
  ggplot(aes(x = Trait, y = observed)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), linewidth = 0.5,width = 0.2, color = "gray50") +
  geom_point(aes(color = significant), size = 6) +
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  ylim(0,1)+
  scale_y_continuous(breaks = seq(0, 1, 0.5))+
  coord_flip()+
  facet_wrap( ~Taxonomic.Level,ncol = 3)+
  labs(
    y = "Proportion of Variance Explained",
    color = "Significant"
  ) +
  scale_color_manual(values = c("grey60","black"))+
  theme(strip.background = element_blank(),
        axis.line.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        #axis.title.x = element_blank(),
        #strip.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "top",
        legend.direction = "horizontal"
  )




#### Correlation with climatic variables for AUC ####
#
#
# Here we test the correlation between thermal traits and individual bioclimatic
# variables (instead of broad Climate Profiles) using generalized linear models.
# We test all possible additive combinations of several orthogonal variables and 
# compare their AIC values to identify the most likely influential variables on 
# specific thermal traits
#
#


#### S13.top.left - Correlation of thermal traits and bioclimatic variables for auc ####
correlation_matrix_auc <- ave_params_gcplyr_combined_common_taxa_auc %>% 
  dplyr::select(c("seqID","Order","Family","Genus","temperature_of_isolation","Site_initials","rmax",     
                  "topt" ,"ctmin","ctmax","thermal_safety_margin","thermal_tolerance",
                  "Seasonality","MAT", "MATR","MMTR","MATmax","MMTmax","MATmin","MMTmin","Isothermality","Climate profile",
                  "Annual_precipitation","Precipitation_seasonality","MAPmax","GC%","Total length"))


correlation_matrix_auc_3 <- correlation_matrix_auc %>% dplyr::select(c("rmax", "topt" ,"ctmin","ctmax",
                                                                       "Seasonality","MAT", "MATR","MMTR","MATmax","MMTmax","MATmin","MMTmin","Isothermality",
                                                                       "Annual_precipitation","Precipitation_seasonality","MAPmax"))

correlation_matrix_auc_3 <- correlation_matrix_auc_3 %>% 
  rename(
    "Rmax" = "rmax",
    "Precipitation seasonality" = "Precipitation_seasonality",
    "Annual precipitation" = "Annual_precipitation",
    "Temperature Seasonality" = "Seasonality",
    "CTmax" = "ctmax",
    "CTmin" = "ctmin",
    "Topt" = "topt"
  )

f_auc <- cor(correlation_matrix_auc_3)


library(ggcorrplot)

S13.top.left <- ggcorrplot(f_auc, hc.order = FALSE, type = "upper",
                                                     outline.col = "white",
                                                     ggtheme = ggplot2::theme_classic(
                                                       base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                                                     colors = c("#6D9EC1", "white", "#E46726"),
                                                     lab = TRUE,
                                                     insig = "blank"
)+
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = c(0.9,0.3),
        legend.direction = "vertical",
        legend.text = element_text(angle = 0))+
  annotate(
    "rect",
    xmin = 0.5, xmax = 4.5,
    ymin = 3.5, ymax = 15.5,
    color = "black", fill = "transparent",
    linewidth = 1)




#### S13.top.center and right - Correlation of thermal traits and bioclimatic variables for auc by Order ####

correlation_matrix_auc_Enterobacterales <- correlation_matrix_auc %>% filter(Order == "Enterobacterales")
correlation_matrix_auc_Pseudomonadales <- correlation_matrix_auc %>% filter(Order == "Pseudomonadales")

correlation_matrix_auc_3_Enterobacterales <- correlation_matrix_auc_Enterobacterales %>% dplyr::select(c("rmax", "topt" ,"ctmin","ctmax",
                                                                       "Seasonality","MAT", "MATR","MMTR","MATmax","MMTmax","MATmin","MMTmin","Isothermality",
                                                                       "Annual_precipitation","Precipitation_seasonality","MAPmax"))

correlation_matrix_auc_3_Enterobacterales <- correlation_matrix_auc_3_Enterobacterales %>% 
  rename(
    "Rmax" = "rmax",
    "Precipitation seasonality" = "Precipitation_seasonality",
    "Annual precipitation" = "Annual_precipitation",
    "Temperature Seasonality" = "Seasonality",
    "CTmax" = "ctmax",
    "CTmin" = "ctmin",
    "Topt" = "topt"
  )


f_auc_Enterobacterales <- cor(correlation_matrix_auc_3_Enterobacterales)


# Get the upper triangle
S13.top.center <- ggcorrplot(f_auc_Enterobacterales, hc.order = FALSE, type = "upper",
                                                     outline.col = "white",
                                                     ggtheme = ggplot2::theme_classic(
                                                       base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                                                     colors = c("#6D9EC1", "white", "#E46726"),
                                                     lab = TRUE,
                                                     insig = "blank"
)+
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = c(0.9,0.3),
        legend.direction = "vertical",
        legend.text = element_text(angle = 0))+
  annotate(
    "rect",
    xmin = 0.5, xmax = 4.5,
    ymin = 3.5, ymax = 15.5,
    color = "black", fill = "transparent",
    linewidth = 1)


correlation_matrix_auc_3_Pseudomonadales <- correlation_matrix_auc_Pseudomonadales %>% dplyr::select(c("rmax", "topt" ,"ctmin","ctmax",
                                                                                                           "Seasonality","MAT", "MATR","MMTR","MATmax","MMTmax","MATmin","MMTmin","Isothermality",
                                                                                                           "Annual_precipitation","Precipitation_seasonality","MAPmax"))

correlation_matrix_auc_3_Pseudomonadales <- correlation_matrix_auc_3_Pseudomonadales %>% 
  rename(
    "Rmax" = "rmax",
    "Precipitation seasonality" = "Precipitation_seasonality",
    "Annual precipitation" = "Annual_precipitation",
    "Temperature Seasonality" = "Seasonality",
    "CTmax" = "ctmax",
    "CTmin" = "ctmin",
    "Topt" = "topt"
  )

f_auc_Pseudomonadales <- cor(correlation_matrix_auc_3_Pseudomonadales)



# Get the upper triangle
S13.top.right <- ggcorrplot(f_auc_Pseudomonadales, hc.order = FALSE, type = "upper",
                                                                      outline.col = "white",
                                                                      ggtheme = ggplot2::theme_classic(
                                                                        base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                                                                      colors = c("#6D9EC1", "white", "#E46726"),
                                                                      lab = TRUE,
                                                                      insig = "blank"
)+
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = c(0.9,0.3),
        legend.direction = "vertical",
        legend.text = element_text(angle = 0))+
  annotate(
    "rect",
    xmin = 0.5, xmax = 4.5,
    ymin = 3.5, ymax = 15.5,
    color = "black", fill = "transparent",
    linewidth = 1)


#### Correlation matrix for AUC ####
correlation_matrix_auc <- ave_params_gcplyr_combined_major_taxa_auc %>% 
  dplyr::select(c("seqID","Order","Family","Genus","temperature_of_isolation","Site_initials","rmax",     
                  "topt" ,"ctmin","ctmax","thermal_safety_margin","thermal_tolerance","PC1","PC2"))



correlation_matrix_auc_2 <- correlation_matrix_auc %>% dplyr::select(c("temperature_of_isolation","Order","Genus","rmax", "topt" ,"ctmin","ctmax","thermal_safety_margin","thermal_tolerance",
                                                                       "PC1","PC2"))


hist(log(correlation_matrix_auc_2$rmax))
hist(correlation_matrix_auc_2$topt)
hist(correlation_matrix_auc_2$ctmax)
hist(correlation_matrix_auc_2$ctmin)


#### Running LMERs for AUC ####
#
#
# Here we will fit hierarchical linear mixed effect models for each thermal trait
# We will evaluate the random effects of nested phylogenetic groupings 
# and the fixed effects of the PCs obtained from the ordination of environmental variables
# with a formula of the type "Trait_lmer = Trait ~ PC1 + PC2 + (1 | Order/Genus)
#
#

# Then we run a model for each trait. For rmax we first transform its distribution to a more normalized shape using log10. 
# This transformation gets us closer to a normal shape, but the Shapiro-Wilk test is very sensitive, so it still fails

lmer_auc_rmax <- correlation_matrix_auc_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmin, ctmax))

model_lmer_auc_rmax <- lmer(log10(rmax) ~ PC1 + PC2 + Order + temperature_of_isolation + (1 | Genus/Order), data = lmer_auc_rmax)
plot(model_lmer_auc_rmax)              # residuals vs fitted
qqnorm(resid(model_lmer_auc_rmax))     # normality of residuals
qqline(resid(model_lmer_auc_rmax))
summary(model_lmer_auc_rmax)
                     

lmer_auc_topt <- correlation_matrix_auc_2 %>% dplyr::select(!c(rmax, thermal_tolerance,thermal_safety_margin, ctmin, ctmax))

model_lmer_auc_topt <- lmer(topt ~ PC1 + PC2 + Order  + temperature_of_isolation + (1 | Genus/Order), data = lmer_auc_topt)
plot(model_lmer_auc_topt)              # residuals vs fitted
qqnorm(resid(model_lmer_auc_topt))     # normality of residuals
qqline(resid(model_lmer_auc_topt))
summary(model_lmer_auc_topt)


lmer_auc_ctmax <- correlation_matrix_auc_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmin, rmax))

model_lmer_auc_ctmax <- lmer(ctmax ~ PC1 + PC2 + Order + temperature_of_isolation + (1 | Genus/Order), data = lmer_auc_ctmax)
plot(model_lmer_auc_ctmax)              # residuals vs fitted
qqnorm(resid(model_lmer_auc_ctmax))     # normality of residuals
qqline(resid(model_lmer_auc_ctmax))
summary(model_lmer_auc_ctmax)


lmer_auc_ctmin <- correlation_matrix_auc_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, rmax, ctmax))

model_lmer_auc_ctmin <- lmer(ctmin ~ PC1 + PC2 + Order + temperature_of_isolation + (1 | Genus/Order), data = lmer_auc_ctmin)
plot(model_lmer_auc_ctmin)              # residuals vs fitted
qqnorm(resid(model_lmer_auc_ctmin))     # normality of residuals
qqline(resid(model_lmer_auc_ctmin))
summary(model_lmer_auc_ctmin)

# Function to extract coefficients for a given model and trait name
extract_coefficients_lmer <- function(model, trait_name) {
  coefs <- summary(model)$coefficients
  c(
    Trait = trait_name,
    PC1_mean = coefs["PC1", "Estimate"],
    PC1_sd = coefs["PC1", "Std. Error"],
    PC2_mean = coefs["PC2", "Estimate"],
    PC2_sd = coefs["PC2", "Std. Error"],
    Order_mean = coefs["OrderPseudomonadales", "Estimate"],
    Order_sd = coefs["OrderPseudomonadales", "Std. Error"],
    T_of_iso_mean = coefs["temperature_of_isolation", "Estimate"],
    T_of_iso_sd = coefs["temperature_of_isolation", "Std. Error"]
  )
}

# Extract from each model
lmer_auc_ctmin_coefficients <- extract_coefficients_lmer(model_lmer_auc_ctmin, "CTmin")
lmer_auc_ctmax_coefficients <- extract_coefficients_lmer(model_lmer_auc_ctmax, "CTmax")
lmer_auc_topt_coefficients  <- extract_coefficients_lmer(model_lmer_auc_topt,  "Topt")
lmer_auc_rmax_coefficients  <- extract_coefficients_lmer(model_lmer_auc_rmax,  "rmax")


# Combine into a data frame
lmer_auc_output <- data.frame(bind_rows(lmer_auc_rmax_coefficients, 
                                    lmer_auc_topt_coefficients, 
                                    lmer_auc_ctmax_coefficients, 
                                    lmer_auc_ctmin_coefficients), 
                          stringsAsFactors = FALSE)

# Convert numeric columns from character (after cbind) back to numeric
lmer_auc_output[, 2:9] <- lapply(lmer_auc_output[, 2:9], as.numeric)


library(partR2)
partR2(model_lmer_auc_rmax, partvars = c("PC1", "PC2", "Order", "temperature_of_isolation"))
partR2(model_lmer_auc_topt, partvars = c("PC1", "PC2", "Order", "temperature_of_isolation"))
partR2(model_lmer_auc_ctmax, partvars = c("PC1", "PC2", "Order", "temperature_of_isolation"))
partR2(model_lmer_auc_ctmin, partvars = c("PC1", "PC2", "Order", "temperature_of_isolation"))


#### Running GLMs for AUC ####
#
# Here we run a similar model, but using a generalized linear model syntax
# This model does not allow for nested variables and random effects, but the results are
# very similar to he ones we just obtained 
#
#
#### GLM max auc ####
glm_aucmax <- correlation_matrix_auc_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmin, ctmax))

# We first get a model with all the variables
full_model_rmax_auc_full_data <- glm(rmax ~ ., data = glm_aucmax, na.action = na.fail)

# And then we perform model selection with the function dredge
all_models_rmax_auc_full_data <- dredge(full_model_rmax_auc_full_data)

# We finally rank all models by AICc
best_models_rmax_auc_full_data <- subset(all_models_rmax_auc_full_data, delta < 2)  

# And check the best model/s
best_models_rmax_auc_full_data

summary(glm(rmax ~ Order + PC1 + PC2 + temperature_of_isolation + Genus, data = glm_aucmax, na.action = na.fail))

glm_max_results_auc <- (glm(rmax ~ Order + PC1 + PC2 + temperature_of_isolation + Genus, data = glm_aucmax, na.action = na.fail))
rsq.partial(glm_max_results_auc)


df_glm_max_results_auc <- data.frame(
  "metric" = "auc",
  "Trait" = "rmax",
  "Variable" = rsq.partial(glm_max_results_auc)$variable,
  "Partial R2" = rsq.partial(glm_max_results_auc)$partial.rsq
)


#### GLM ctmax auc ####
glm_ctmax <- correlation_matrix_auc_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmin, rmax ))

full_model_ctmax_auc_full_data <- glm(ctmax ~ ., data = glm_ctmax, na.action = na.fail)

# And then we perform model selection with the function dredge
all_models_ctmax_auc_full_data <- dredge(full_model_ctmax_auc_full_data)

# We finally rank all models by AICc
best_models_ctmax_auc_full_data <- subset(all_models_ctmax_auc_full_data, delta < 2)  
# And check the best model/s
best_models_ctmax_auc_full_data

summary(glm(ctmax ~ Order + PC1 + PC2 + temperature_of_isolation + Genus, data = glm_ctmax, na.action = na.fail))
summary(glm(ctmax ~ Order + PC2 + temperature_of_isolation, data = glm_ctmax, na.action = na.fail))

glm_ctmax_results_auc <- (glm(ctmax ~ Order + PC1 + PC2 + temperature_of_isolation + Genus, data = glm_ctmax, na.action = na.fail))
rsq.partial(glm_ctmax_results_auc)

df_glm_ctmax_results_auc <- data.frame(
  "metric" = "auc",
  "Trait" = "ctmax",
  "Variable" = rsq.partial(glm_ctmax_results_auc)$variable,
  "Partial R2" = rsq.partial(glm_ctmax_results_auc)$partial.rsq
)

#### GLM ctmin auc ####
glm_ctmin <- correlation_matrix_auc_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmax, rmax ))

full_model_ctmin_auc_full_data <- glm(ctmin ~ ., data = glm_ctmin, na.action = na.fail)

# And then we perform model selection with the function dredge
all_models_ctmin_auc_full_data <- dredge(full_model_ctmin_auc_full_data)

# We finally rank all models by AICc
best_models_ctmin_auc_full_data <- subset(all_models_ctmin_auc_full_data, delta < 2)  

# And check the best model/s
best_models_ctmin_auc_full_data

summary(glm(ctmin ~ Order + PC1 + PC2 + temperature_of_isolation + Genus, data = glm_ctmin, na.action = na.fail))
glm_ctmin_results_auc <- (glm(ctmin ~ Order + PC1 + PC2 + temperature_of_isolation + Genus, data = glm_ctmin, na.action = na.fail))
rsq.partial(glm_ctmin_results_auc)

df_glm_ctmin_results_auc <- data.frame(
  "metric" = "auc",
  "Trait" = "ctmin",
  "Variable" = rsq.partial(glm_ctmin_results_auc)$variable,
  "Partial R2" = rsq.partial(glm_ctmin_results_auc)$partial.rsq
)

#### GLM topt auc ####
glm_topt <- correlation_matrix_auc_2 %>% dplyr::select(!c(ctmin, thermal_tolerance,thermal_safety_margin, ctmax, rmax ))
full_model_topt_auc_full_data <- glm(topt ~ ., data = glm_topt, na.action = na.fail)

# And then we perform model selection with the function dredge
all_models_topt_auc_full_data <- dredge(full_model_topt_auc_full_data)

# We finally rank all models by AICc
best_models_topt_auc_full_data <- subset(all_models_topt_auc_full_data, delta < 2)  

# And check the best model/s
best_models_topt_auc_full_data

summary(glm(topt ~ Order + PC1 + PC2 + temperature_of_isolation + Genus, data = glm_topt, na.action = na.fail))

glm_topt_results_auc <- (glm(topt ~ Order + PC1 + PC2 + temperature_of_isolation + Genus, data = glm_topt, na.action = na.fail))
rsq.partial(glm_topt_results_auc)

df_glm_topt_results_auc <- data.frame(
  "metric" = "auc",
  "Trait" = "topt",
  "Variable" = rsq.partial(glm_topt_results_auc)$variable,
  "Partial R2" = rsq.partial(glm_topt_results_auc)$partial.rsq
)



glm_table_auc <- bind_rows(df_glm_max_results_auc,
                           df_glm_topt_results_auc,
                           df_glm_ctmax_results_auc,
                           df_glm_ctmin_results_auc)


combined_effects_lmer_auc <- fread("combined_effects_lmer_auc.csv")



combined_effects_glm_lmer_auc <- left_join(combined_effects_lmer_auc,glm_table_auc, by = c("Trait","metric","Variable"))
combined_effects_glm_lmer_auc$Partial_R2 <- as.numeric(combined_effects_glm_lmer_auc$Partial_R2)
combined_effects_glm_lmer_auc$Variance_observed <- as.numeric(combined_effects_glm_lmer_auc$Variance_observed)


combined_effects_glm_lmer_auc$pR2_dif <- round(combined_effects_glm_lmer_auc$Partial_R2 - combined_effects_glm_lmer_auc$Partial.R2,digits = 3)
combined_effects_glm_lmer_auc$Var_dif <- round(combined_effects_glm_lmer_auc$Variance_observed - combined_effects_glm_lmer_auc$Partial.R2,digits = 3)


#### S13.bot.left - Correlation of thermal traits and bioclimatic variables for rate ####
#
#
# Here we test the correlation between thermal traits and individual bioclimatic
# variables (instead of broad Climate Profiles) using generalized linear models.
# We test all possible additive combinations of several orthogonal variables and 
# compare their AIC values to identify the most likely influential variables on 
# specific thermal traits
#
#


correlation_matrix_rate <- ave_params_gcplyr_combined_common_taxa_rate %>% 
  dplyr::select(c("seqID","Order","Family","Genus","temperature_of_isolation","Site_initials","rmax",     
                  "topt" ,"ctmin","ctmax","thermal_safety_margin","thermal_tolerance",
                  "Seasonality","MAT", "MATR","MMTR","MATmax","MMTmax","MATmin","MMTmin","Isothermality","Climate profile",
                  "Annual_precipitation","Precipitation_seasonality","MAPmax","GC%","Total length"))


correlation_matrix_rate_3 <- correlation_matrix_rate %>% dplyr::select(c("rmax", "topt" ,"ctmin","ctmax",
                                                                       "Seasonality","MAT", "MATR","MMTR","MATmax","MMTmax","MATmin","MMTmin","Isothermality",
                                                                       "Annual_precipitation","Precipitation_seasonality","MAPmax"))

correlation_matrix_rate_3 <- correlation_matrix_rate_3 %>% 
  rename(
    "Rmax" = "rmax",
    "Precipitation seasonality" = "Precipitation_seasonality",
    "Annual precipitation" = "Annual_precipitation",
    "Temperature Seasonality" = "Seasonality",
    "CTmax" = "ctmax",
    "CTmin" = "ctmin",
    "Topt" = "topt"
  )

f_rate <- cor(correlation_matrix_rate_3)


library(ggcorrplot)

S13.bot.left <- ggcorrplot(f_rate, hc.order = FALSE, type = "upper",
                        outline.col = "white",
                        ggtheme = ggplot2::theme_classic(
                          base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                        colors = c("#6D9EC1", "white", "#E46726"),
                        lab = TRUE,
                        insig = "blank"
)+
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = c(0.9,0.3),
        legend.direction = "vertical",
        legend.text = element_text(angle = 0))+
  annotate(
    "rect",
    xmin = 0.5, xmax = 4.5,
    ymin = 3.5, ymax = 15.5,
    color = "black", fill = "transparent",
    linewidth = 1)



#### S13.bot.center and right - Correlation of thermal traits and bioclimatic variables for rate by Order ####

correlation_matrix_rate_Enterobacterales <- correlation_matrix_rate %>% filter(Order == "Enterobacterales")
correlation_matrix_rate_Pseudomonadales <- correlation_matrix_rate %>% filter(Order == "Pseudomonadales")

correlation_matrix_rate_3_Enterobacterales <- correlation_matrix_rate_Enterobacterales %>% dplyr::select(c("rmax", "topt" ,"ctmin","ctmax",
                                                                                                         "Seasonality","MAT", "MATR","MMTR","MATmax","MMTmax","MATmin","MMTmin","Isothermality",
                                                                                                         "Annual_precipitation","Precipitation_seasonality","MAPmax"))

correlation_matrix_rate_3_Enterobacterales <- correlation_matrix_rate_3_Enterobacterales %>% 
  rename(
    "Rmax" = "rmax",
    "Precipitation seasonality" = "Precipitation_seasonality",
    "Annual precipitation" = "Annual_precipitation",
    "Temperature Seasonality" = "Seasonality",
    "CTmax" = "ctmax",
    "CTmin" = "ctmin",
    "Topt" = "topt"
  )


f_rate_Enterobacterales <- cor(correlation_matrix_rate_3_Enterobacterales)


# Get the upper triangle
S13.bot.center <- ggcorrplot(f_rate_Enterobacterales, hc.order = FALSE, type = "upper",
                        outline.col = "white",
                        ggtheme = ggplot2::theme_classic(
                          base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                        colors = c("#6D9EC1", "white", "#E46726"),
                        lab = TRUE,
                        insig = "blank"
)+
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = c(0.9,0.3),
        legend.direction = "vertical",
        legend.text = element_text(angle = 0))+
  annotate(
    "rect",
    xmin = 0.5, xmax = 4.5,
    ymin = 3.5, ymax = 15.5,
    color = "black", fill = "transparent",
    linewidth = 1)


correlation_matrix_rate_3_Pseudomonadales <- correlation_matrix_rate_Pseudomonadales %>% dplyr::select(c("rmax", "topt" ,"ctmin","ctmax",
                                                                                                       "Seasonality","MAT", "MATR","MMTR","MATmax","MMTmax","MATmin","MMTmin","Isothermality",
                                                                                                       "Annual_precipitation","Precipitation_seasonality","MAPmax"))

correlation_matrix_rate_3_Pseudomonadales <- correlation_matrix_rate_3_Pseudomonadales %>% 
  rename(
    "Rmax" = "rmax",
    "Precipitation seasonality" = "Precipitation_seasonality",
    "Annual precipitation" = "Annual_precipitation",
    "Temperature Seasonality" = "Seasonality",
    "CTmax" = "ctmax",
    "CTmin" = "ctmin",
    "Topt" = "topt"
  )

f_rate_Pseudomonadales <- cor(correlation_matrix_rate_3_Pseudomonadales)



# Get the upper triangle
S13.bot.right <- ggcorrplot(f_rate_Pseudomonadales, hc.order = FALSE, type = "upper",
                        outline.col = "white",
                        ggtheme = ggplot2::theme_classic(
                          base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                        colors = c("#6D9EC1", "white", "#E46726"),
                        lab = TRUE,
                        insig = "blank"
)+
  theme(panel.border = element_rect(color = "black", fill = NA),
        legend.position = c(0.9,0.3),
        legend.direction = "vertical",
        legend.text = element_text(angle = 0))+
  annotate(
    "rect",
    xmin = 0.5, xmax = 4.5,
    ymin = 3.5, ymax = 15.5,
    color = "black", fill = "transparent",
    linewidth = 1)



#### Correlation matrix for Rate ####
correlation_matrix_rate <- ave_params_gcplyr_combined_major_taxa_rate %>% 
  dplyr::select(c("seqID","Order","Family","Genus","temperature_of_isolation","Site_initials","rmax",     
                  "topt" ,"ctmin","ctmax","thermal_safety_margin","thermal_tolerance","PC1","PC2"))



correlation_matrix_rate_2 <- correlation_matrix_rate %>% dplyr::select(c("Order","Genus","rmax", "topt" ,"ctmin","ctmax","thermal_safety_margin","thermal_tolerance",
                                                                       "PC1","PC2"))


hist(log(correlation_matrix_rate_2$rmax))
hist(correlation_matrix_rate_2$topt)
hist(correlation_matrix_rate_2$ctmax)
hist(correlation_matrix_rate_2$ctmin)


#### Running LMERs for Rate ####
#
#
# Here we will fit hierarchical linear mixed effect models for each thermal trait
# We will evaluate the random effects of nested phylogenetic groupings 
# and the fixed effects of the PCs obtained from the ordination of environmental variables
# with a formula of the type "Trait_lmer = Trait ~ PC1 + PC2 + (1 | Order/Genus)
#
#

# First we build an empty dataframe to store our results 

glm_rate_output <- data.frame(
  Trait = character(),
  PC1_mean = numeric(),
  PC1_sd = numeric(),
  PC2_mean = numeric(),
  PC2_sd = numeric(),
  stringsAsFactors = FALSE
)


# Then we run a model for each trait. For rmax we first transform its distribution to a more normalized shape using log10. 
# This transformation gets us closer to a normal shape, but the Shapiro-Wilk test is very sensitive, so it still fails

glm_rate_rmax <- correlation_matrix_rate_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmin, ctmax))

model_glm_rate_rmax <- lmer(log10(rmax) ~ PC1 + PC2 + (1 | Order/Genus), data = glm_rate_rmax)
plot(model_glm_rate_rmax)              # residuals vs fitted
qqnorm(resid(model_glm_rate_rmax))     # normality of residuals
qqline(resid(model_glm_rate_rmax))
summary(model_glm_rate_rmax)


glm_rate_topt <- correlation_matrix_rate_2 %>% dplyr::select(!c(rmax, thermal_tolerance,thermal_safety_margin, ctmin, ctmax))

model_glm_rate_topt <- lmer(topt ~ PC1 + PC2 + (1 | Order/Genus), data = glm_rate_topt)
plot(model_glm_rate_topt)              # residuals vs fitted
qqnorm(resid(model_glm_rate_topt))     # normality of residuals
qqline(resid(model_glm_rate_topt))
summary(model_glm_rate_topt)


glm_rate_ctmax <- correlation_matrix_rate_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmin, rmax))

model_glm_rate_ctmax <- lmer(ctmax ~ PC1 + PC2 + (1 | Order/Genus), data = glm_rate_ctmax)
plot(model_glm_rate_ctmax)              # residuals vs fitted
qqnorm(resid(model_glm_rate_ctmax))     # normality of residuals
qqline(resid(model_glm_rate_ctmax))
summary(model_glm_rate_ctmax)


glm_rate_ctmin <- correlation_matrix_rate_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, rmax, ctmax))

model_glm_rate_ctmin <- lmer(ctmin ~ PC1 + PC2 + (1 | Order/Genus), data = glm_rate_ctmin)
plot(model_glm_rate_ctmin)              # residuals vs fitted
qqnorm(resid(model_glm_rate_ctmin))     # normality of residuals
qqline(resid(model_glm_rate_ctmin))
summary(model_glm_rate_ctmin)

# Function to extract coefficients for a given model and trait name
extract_coefficients_lmer <- function(model, trait_name) {
  coefs <- summary(model)$coefficients
  c(
    Trait = trait_name,
    PC1_mean = coefs["PC1", "Estimate"],
    PC1_sd = coefs["PC1", "Std. Error"],
    PC2_mean = coefs["PC2", "Estimate"],
    PC2_sd = coefs["PC2", "Std. Error"]
  )
}

# Extract from each model
glm_rate_ctmin_coefficients <- extract_coefficients_lmer(model_glm_rate_ctmin, "CTmin")
glm_rate_ctmax_coefficients <- extract_coefficients_lmer(model_glm_rate_ctmax, "CTmax")
glm_rate_topt_coefficients  <- extract_coefficients_lmer(model_glm_rate_topt,  "Topt")
glm_rate_rmax_coefficients  <- extract_coefficients_lmer(model_glm_rate_rmax,  "rmax")


# Combine into a data frame
glm_rate_output <- data.frame(bind_rows(glm_rate_rmax_coefficients, 
                                       glm_rate_topt_coefficients, 
                                       glm_rate_ctmax_coefficients, 
                                       glm_rate_ctmin_coefficients), 
                             stringsAsFactors = FALSE)

# Convert numeric columns from character (after cbind) back to numeric
glm_rate_output[, 2:5] <- lapply(glm_rate_output[, 2:5], as.numeric)






#### Running GLMs for Rate ####
#
# Here we run a similar model, but using a generalized linear model syntax
# This model does not allow for nested variables and random effects, but the results are
# very similar to he ones we just obtained 
#
#
#### GLM max rate ####
glm_ratemax <- correlation_matrix_rate_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmin, ctmax))

# We first get a model with all the variables
full_model_rmax_rate_full_data <- glm(rmax ~ ., data = glm_ratemax, na.action = na.fail)

# And then we perform model selection with the function dredge
all_models_rmax_rate_full_data <- dredge(full_model_rmax_rate_full_data)

# We finally rank all models by AICc
best_models_rmax_rate_full_data <- subset(all_models_rmax_rate_full_data, delta < 2)  

# And check the best model/s
best_models_rmax_rate_full_data


#### GLM ctmax rate ####
glm_ctmax <- correlation_matrix_rate_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmin, rmax ))

full_model_ctmax_rate_full_data <- glm(ctmax ~ ., data = glm_ctmax, na.action = na.fail)

# And then we perform model selection with the function dredge
all_models_ctmax_rate_full_data <- dredge(full_model_ctmax_rate_full_data)

# We finally rank all models by AICc
best_models_ctmax_rate_full_data <- subset(all_models_ctmax_rate_full_data, delta < 2)  

# And check the best model/s
best_models_ctmax_rate_full_data

#### GLM ctmin rate ####
glm_ctmin <- correlation_matrix_rate_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmax, rmax ))

full_model_ctmin_rate_full_data <- glm(ctmin ~ ., data = glm_ctmin, na.action = na.fail)

# And then we perform model selection with the function dredge
all_models_ctmin_rate_full_data <- dredge(full_model_ctmin_rate_full_data)

# We finally rank all models by AICc
best_models_ctmin_rate_full_data <- subset(all_models_ctmin_rate_full_data, delta < 2)  

# And check the best model/s
best_models_ctmin_rate_full_data


#### GLM topt rate ####
glm_topt <- correlation_matrix_rate_2 %>% dplyr::select(!c(ctmin, thermal_tolerance,thermal_safety_margin, ctmax, rmax ))
full_model_topt_rate_full_data <- glm(topt ~ ., data = glm_topt, na.action = na.fail)

# And then we perform model selection with the function dredge
all_models_topt_rate_full_data <- dredge(full_model_topt_rate_full_data)

# We finally rank all models by AICc
best_models_topt_rate_full_data <- subset(all_models_topt_rate_full_data, delta < 2)  

# And check the best model/s
best_models_topt_rate_full_data



#### Getting Areas of TPC triangles  ####
#
#
# Here we simplify the TPCs by only selecting their limits (CTmin and CTmax) and peak value (rmax) 
# as vertices of a theoretical performance triangle. We then calculate their area and summarize by 
# metric and climate profile
#
#
ave_params_gcplyr_combined_major_taxa$upper_triangle <- ((ave_params_gcplyr_combined_major_taxa$ctmax - ave_params_gcplyr_combined_major_taxa$topt)*ave_params_gcplyr_combined_major_taxa$rmax)/2 
ave_params_gcplyr_combined_major_taxa$lower_triangle <- ((ave_params_gcplyr_combined_major_taxa$topt - ave_params_gcplyr_combined_major_taxa$ctmin)*ave_params_gcplyr_combined_major_taxa$rmax)/2 
ave_params_gcplyr_combined_major_taxa$triangle <- ((ave_params_gcplyr_combined_major_taxa$ctmax - ave_params_gcplyr_combined_major_taxa$ctmin)*ave_params_gcplyr_combined_major_taxa$rmax)/2 


#### S8 - Plotting the TPC triangles by Climate Profile ####

# Here we need to obtain the triangle shapes by mapping the coordinates 
# of each vertex to x and y columns of a new dataframe 

full_triangle_data <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,Family,Genus,`Climate profile`,topt,rmax,ctmin,ctmax,metric) %>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(seqID,metric) %>%
  summarize(
    x = list(c(ctmax,ctmin, topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    `Climate profile` = first(`Climate profile`)
  ) %>%
  unnest(cols = c(x, y)) %>%
  mutate(group = interaction(seqID,metric))  

upper_triangle_data <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,Family,Genus,`Climate profile`,rmax,topt,ctmin,ctmax,metric) %>%
  distinct() %>%  
  group_by(seqID,metric) %>%
  summarize(
    x = list(c(ctmax, topt, topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    `Climate profile` = first(`Climate profile`)
  ) %>%
  unnest(cols = c(x, y)) %>%  
  mutate(group = interaction(seqID,metric))  

lower_triangle_data <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,Family,Genus,`Climate profile`,rmax,topt,ctmin,ctmax,metric) %>%
  distinct() %>%  
  group_by(seqID,metric) %>%
  summarize(
    x = list(c(ctmin, topt, topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),    `Climate profile` = first(`Climate profile`)
  ) %>%
  unnest(cols = c(x, y)) %>%  
  mutate(group = interaction(seqID,metric))  


full_triangle_data_Order <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,`Climate profile`, rmax, topt,ctmax,ctmin,metric) %>%
  group_by(`Climate profile`,Order,metric) %>%
  summarize(sd_topt = sd(topt)/mean(topt),
            sd_max = sd(rmax)/mean(rmax),
            sd_ctmax = sd(ctmax)/mean(ctmax),
            sd_ctmin = sd(ctmin)/mean(ctmin),
            topt = mean(topt),
            rmax = mean(rmax),
            ctmax = mean(ctmax),
            ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(`Climate profile`,Order,metric) %>%
  summarize(
    x = list(c(ctmax,ctmin,topt,topt)),
    y = list(c(0,0,0, rmax)),
    z = list(c(sd_ctmax,sd_ctmin,sd_topt,sd_max)),
    Order = first(Order),
    `Climate profile` = first(`Climate profile`)
  ) %>%
  unnest(cols = c(x,y,z)) %>%  # Unnest to get a tidy dataframe for plotting
  mutate(group = interaction(Order,metric))  

upper_triangle_data_Order <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,`Climate profile`,rmax, topt,ctmax,ctmin,metric) %>%
  group_by(`Climate profile`,Order,metric) %>%
  summarize(sd_topt = sd(topt)/mean(topt),
            sd_max = sd(rmax)/mean(rmax),
            sd_ctmax = sd(ctmax)/mean(ctmax),
            sd_ctmin = sd(ctmin)/mean(ctmin),
            topt = mean(topt),
            rmax = mean(rmax),
            ctmax = mean(ctmax),
            ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(`Climate profile`,Order,metric) %>%
  summarize(
    x = list(c(ctmax,topt,topt)),
    y = list(c(0,0,rmax)),
    Order = first(Order),
    `Climate profile` = first(`Climate profile`)
  ) %>%
  unnest(cols = c(x,y)) %>%  # Unnest to get a tidy dataframe for plotting
  mutate(group = interaction(Order,metric))  

lower_triangle_data_Order <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,`Climate profile`, rmax, topt, ctmax,ctmin,metric) %>%
  group_by(`Climate profile`,Order,metric) %>%
  summarize(sd_topt = sd(topt)/mean(topt),
            sd_max = sd(rmax)/mean(rmax),
            sd_ctmax = sd(ctmax)/mean(ctmax),
            sd_ctmin = sd(ctmin)/mean(ctmin),
            topt = mean(topt),
            rmax = mean(rmax),
            ctmax = mean(ctmax),
            ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(`Climate profile`,Order,metric) %>%
  summarize(
    x = list(c(ctmin,topt,topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order),
    `Climate profile` = first(`Climate profile`)
  ) %>%
  unnest(cols = c(x,y)) %>%  # Unnest to get a tidy dataframe for plotting
  mutate(group = interaction(Order,metric))  


full_triangle_data$`Climate profile` <- factor(full_triangle_data$`Climate profile`)


S8.left <- full_triangle_data %>% 
  filter(metric == "auc") %>%
  ggplot(aes(x = x, y = y, group = group, fill = Order)) +
  geom_vline(xintercept = 42, linetype = "dashed", color = "grey60", linewidth = 0.5)+
  geom_polygon(alpha = 0.2) +
  geom_polygon(data = upper_triangle_data_Order %>% 
                 filter(metric == "auc"),
               aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  geom_polygon(data = lower_triangle_data_Order %>% 
                 filter(metric == "auc"),
               aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  geom_point(data = full_triangle_data_Order %>% 
               filter(metric == "auc"),
             shape = 21,color = "white",aes(x = x, y = y,size = z))+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "AUC"
  ) +
  #facet_wrap(`Climate profile`~Order,ncol= 2)+
  facet_grid(`Climate profile`~Order)+
  #scale_fill_viridis_d(option = "turbo")+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")


S8.right <- full_triangle_data %>% 
  filter(metric == "rate") %>%
  ggplot(aes(x = x, y = y, group = group, fill = Order)) +
  geom_vline(xintercept = 42, linetype = "dashed", color = "grey60", linewidth = 0.5)+
  geom_polygon(alpha = 0.2) +
  geom_polygon(data = upper_triangle_data_Order %>% 
                 filter(metric == "rate"),
               aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  geom_polygon(data = lower_triangle_data_Order %>% 
                 filter(metric == "rate"),
               aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  geom_point(data = full_triangle_data_Order %>% 
               filter(metric == "rate"),
             shape = 21,color = "white",aes(x = x, y = y,size = z))+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "rate"
  ) +
  #facet_wrap(`Climate profile`~Order,ncol= 2)+
  facet_grid(`Climate profile`~Order)+
  #scale_fill_viridis_d(option = "turbo")+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")


#### F3A - Plotting the TPC triangles ####

# Here we need to obtain the triangle shapes by mapping the coordinates 
# of each vertex to x and y columns of a new dataframe 

full_triangle_data_noenv <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,Family,Genus,topt,rmax,ctmin,ctmax,metric) %>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(seqID,metric) %>%
  summarize(
    x = list(c(ctmax,ctmin, topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus)
  ) %>%
  unnest(cols = c(x, y)) %>%
  mutate(group = interaction(seqID,metric))  

upper_triangle_data_noenv <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,Family,Genus,rmax,topt,ctmin,ctmax,metric) %>%
  distinct() %>%  
  group_by(seqID,metric) %>%
  summarize(
    x = list(c(ctmax, topt, topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus)
  ) %>%
  unnest(cols = c(x, y)) %>%  
  mutate(group = interaction(seqID,metric))  

lower_triangle_data_noenv <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,Family,Genus,rmax,topt,ctmin,ctmax,metric) %>%
  distinct() %>%  
  group_by(seqID,metric) %>%
  summarize(
    x = list(c(ctmin, topt, topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus)
  ) %>%
  unnest(cols = c(x, y)) %>%  
  mutate(group = interaction(seqID,metric))  


full_triangle_data_noenv_Order <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order, rmax, topt,ctmax,ctmin,metric) %>%
  group_by(Order,metric) %>%
  summarize(sd_topt = sd(topt)/mean(topt),
            sd_max = sd(rmax)/mean(rmax),
            sd_ctmax = sd(ctmax)/mean(ctmax),
            sd_ctmin = sd(ctmin)/mean(ctmin),
            topt = mean(topt),
            rmax = mean(rmax),
            ctmax = mean(ctmax),
            ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(Order,metric) %>%
  summarize(
    x = list(c(ctmax,ctmin,topt,topt)),
    y = list(c(0,0,0, rmax)),
    z = list(c(sd_ctmax,sd_ctmin,sd_topt,sd_max)),
    Order = first(Order)
  ) %>%
  unnest(cols = c(x,y,z)) %>%  # Unnest to get a tidy dataframe for plotting
  mutate(group = interaction(Order,metric))  

upper_triangle_data_noenv_Order <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,rmax, topt,ctmax,ctmin,metric) %>%
  group_by(Order,metric) %>%
  summarize(sd_topt = sd(topt)/mean(topt),
            sd_max = sd(rmax)/mean(rmax),
            sd_ctmax = sd(ctmax)/mean(ctmax),
            sd_ctmin = sd(ctmin)/mean(ctmin),
            topt = mean(topt),
            rmax = mean(rmax),
            ctmax = mean(ctmax),
            ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(Order,metric) %>%
  summarize(
    x = list(c(ctmax,topt,topt)),
    y = list(c(0,0,rmax)),
    Order = first(Order)
  ) %>%
  unnest(cols = c(x,y)) %>%  # Unnest to get a tidy dataframe for plotting
  mutate(group = interaction(Order,metric))  

lower_triangle_data_noenv_Order <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order, rmax, topt, ctmax,ctmin,metric) %>%
  group_by(Order,metric) %>%
  summarize(sd_topt = sd(topt)/mean(topt),
            sd_max = sd(rmax)/mean(rmax),
            sd_ctmax = sd(ctmax)/mean(ctmax),
            sd_ctmin = sd(ctmin)/mean(ctmin),
            topt = mean(topt),
            rmax = mean(rmax),
            ctmax = mean(ctmax),
            ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(Order,metric) %>%
  summarize(
    x = list(c(ctmin,topt,topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order)
  ) %>%
  unnest(cols = c(x,y)) %>%  # Unnest to get a tidy dataframe for plotting
  mutate(group = interaction(Order,metric))  


F3A.auc_noenv <- full_triangle_data_noenv %>% 
  filter(metric == "auc") %>%
  ggplot(aes(x = x, y = y, group = group, fill = Order)) +
  geom_vline(xintercept = 42, linetype = "dashed", color = "grey60", linewidth = 0.5)+
  geom_polygon(alpha = 0.2) +
  geom_polygon(data = upper_triangle_data_noenv_Order %>% 
                 filter(metric == "auc"),
               aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  geom_polygon(data = lower_triangle_data_noenv_Order %>% 
                 filter(metric == "auc"),
               aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  geom_point(data = full_triangle_data_noenv_Order %>% 
               filter(metric == "auc"),
             shape = 21,color = "white",aes(x = x, y = y,size = z))+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "AUC"
  ) +
  facet_wrap(~Order,nrow= 2)+
  #scale_fill_viridis_d(option = "turbo")+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")


F3A.auc_noenv

F3A.rate_noenv <- full_triangle_data_noenv %>% 
  filter(metric == "rate") %>%
  ggplot(aes(x = x, y = y, group = group, fill = Order)) +
  geom_vline(xintercept = 42, linetype = "dashed", color = "grey60", linewidth = 0.5)+
  geom_polygon(alpha = 0.2) +
  geom_polygon(data = upper_triangle_data_noenv_Order %>% 
                 filter(metric == "rate"),
               aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  geom_polygon(data = lower_triangle_data_noenv_Order %>% 
                 filter(metric == "rate"),
               aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  geom_point(data = full_triangle_data_noenv_Order %>% 
               filter(metric == "rate"),
             shape = 21,color = "white",aes(x = x, y = y,size = z))+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "rate"
  ) +
  facet_wrap(~Order,nrow= 2)+
  #scale_fill_viridis_d(option = "turbo")+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(size = "none")


F3A.rate_noenv

#### F3C (Topt - MAT) ####

ave_params_gcplyr_combined_major_taxa <- ave_params_gcplyr_combined_major_taxa %>%
  mutate("Variability" = case_when(
    `Climate profile` == "Cold stable" ~ "Stable",
    `Climate profile` == "Hot stable" ~ "Stable",
    `Climate profile` == "Cold variable" ~ "Variable",
    `Climate profile` == "Hot variable" ~ "Variable"
    ))

ave_params_gcplyr_combined_major_taxa <- ave_params_gcplyr_combined_major_taxa %>%
  mutate("Mean" = case_when(
    `Climate profile` == "Cold stable" ~ "Cold",
    `Climate profile` == "Hot stable" ~ "Hot",
    `Climate profile` == "Cold variable" ~ "Cold",
    `Climate profile` == "Hot variable" ~ "Hot"
  ))

F3C.auc <- ave_params_gcplyr_combined_major_taxa %>%
  filter(metric == "auc")%>%
  ggplot(aes(x=`Climate profile`, y=(topt - MAT)))+
  #scale_color_viridis_b()+
  geom_hline(yintercept = 0,linetype = "dashed",color = "grey60",linewidth = 0.75)+
  geom_jitter(alpha = 0.8,width = 0.2,size= 5,aes(color = Order))+
  geom_boxplot(alpha = 0.8,outliers = FALSE)+
  #facet_wrap(~Order,ncol = 1)+
  #ylim(16,40)+
  #ylim(0.5,3.5)+
  #scale_fill_manual(values = order_color)+
  ylab(expression(Topt[AUC]~" - MAT"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "turbo")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)
  )

S16C <- ave_params_gcplyr_combined_major_taxa %>%
  filter(metric == "rate")%>%
  ggplot(aes(x=`Climate profile`, y=(topt - MAT)))+
  #scale_color_viridis_b()+
  geom_hline(yintercept = 0,linetype = "dashed",color = "grey60",linewidth = 0.75)+
  geom_jitter(alpha = 0.8,width = 0.2,size= 5,aes(color = Order))+
  geom_boxplot(alpha = 0.8,outliers = FALSE)+
  #facet_wrap(~Order,ncol = 1)+
  #ylim(16,40)+
  #ylim(0.5,3.5)+
  #scale_fill_manual(values = order_color)+
  ylab(expression(Topt[rate]~" - MAT"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "turbo")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)
  )

#### F3B top - Randomizing thermal traits for AUC and plotting slope above Topt ####

# Number of simulations
n_sim <- 1000

# Store all results in a list
simulated_slopes_auc <- vector("list", n_sim)

set.seed(123)  # For reproducibility
for (i in 1:n_sim) {
  shuffled_df <- ave_params_gcplyr_combined_major_taxa_auc
  shuffled_df$topt <- sample(shuffled_df$topt)
  shuffled_df$ctmax <- sample(shuffled_df$ctmax)
  shuffled_df$ctmin <- sample(shuffled_df$ctmin)  # optional if you're not using ctmin in slope
  shuffled_df$sim_id <- i
  
  # Calculate slope
  shuffled_df$superopt_slope <- -1 * (shuffled_df$rmax / (shuffled_df$ctmax - shuffled_df$topt))
  
  # Keep only the columns needed for plotting
  simulated_slopes_auc[[i]] <- shuffled_df[, c("rmax", "superopt_slope", "sim_id")]
}

# Combine into one dataframe
simulated_slopes_auc_df <- do.call(rbind, simulated_slopes_auc)

# Remove impossible slopes
simulated_slopes_auc_df <- simulated_slopes_auc_df %>% filter(superopt_slope <= 0) 


model_observed_slopes <- lm(superopt_slope ~ rmax,data = ave_params_gcplyr_combined_major_taxa_auc)
summary(model_observed_slopes)

simulated_slopes_auc_df_max <- simulated_slopes_auc_df %>%
  group_by(rmax)%>%
  summarize(superopt_slope = max(superopt_slope))

model_max_simulated_slopes <- lm(superopt_slope ~ rmax,data = simulated_slopes_auc_df_max)
summary(model_max_simulated_slopes)

model_simulated_slopes <- lm(superopt_slope ~ rmax,data = simulated_slopes_auc_df)
summary(model_simulated_slopes)

F3Btop.auc <- ave_params_gcplyr_combined_major_taxa_auc %>%
  ggplot(aes(x=rmax, y= superopt_slope))+
  geom_point(data = simulated_slopes_auc_df, aes(x = rmax, y = superopt_slope),color = "grey85")+
  #scale_color_viridis_b()+
  geom_smooth(method = "lm", se = TRUE, fill = "grey60", color="grey90") +
  geom_point(alpha = 0.9,aes(color = Order, size = ctmax))+
  #facet_wrap(~`Climate profile`,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  ylim(-1,0)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  ylab(expression("Performance slope above Topt"))+
  xlab(expression("AUCmax"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)
  )

F3Btop.auc


#### F3B bot - Comparing slope ranks AUC ####

slope_data <- ave_params_gcplyr_combined_major_taxa_auc %>%
  dplyr::select(c(superopt_slope,Order,rmax,seqID,topt,ctmax,MMTR,MATR,`Climate profile`))%>%
  rename("observed_slope" = "superopt_slope")

slope_data_ranks <- ave_params_gcplyr_combined_major_taxa_auc %>%
  dplyr::select(c(superopt_slope,Order,rmax,seqID,topt,ctmax,MMTR,MATR,`Climate profile`))%>%
  mutate("sim_id" = 0)

simulated_slopes_auc_df_3 <- simulated_slopes_auc_df %>% 
  left_join(slope_data,by = "rmax")%>%
  dplyr::select(!c(observed_slope))

slope_data_ranks_full <- bind_rows(slope_data_ranks,simulated_slopes_auc_df_3)

slope_data_ranks_full_2 <- slope_data_ranks_full %>%
  group_by(Order, rmax,seqID)%>%
  mutate("slope_rank" = rank(superopt_slope),
         "mean_slope" = mean(superopt_slope),
         "median_slope" = median(superopt_slope),
         "slope_ratio_mean" = superopt_slope/mean_slope,
         "slope_ratio_median" = superopt_slope/median_slope,
         sim_id = sim_id)


slope_data_ranks_full_2

slope_data_ranks_full_3 <- slope_data_ranks_full_2 %>%
  filter(superopt_slope < 0) %>%
  filter(sim_id == 0)

slope_ratio_mean_t <- t.test(slope_data_ranks_full_3$slope_ratio_mean,
       mu = 1,
       alternative = "less")


F3Bbot.auc <- slope_data_ranks_full_2 %>% 
  filter(superopt_slope < 0) %>%
  filter(sim_id == 0) %>%
  ggplot(aes(x = rmax, y = slope_ratio_mean,color = Order)) +
  #geom_smooth()+
  ylim(0,2.5)+
  geom_hline(yintercept = 1,linetype = "dashed",linewidth = 0.5)+
  geom_point(size = 4,alpha = 0.8)+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  ylab(expression("Observed Slope / Predicted Slope"))+
  xlab(expression("AUCmax"))+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)
  )

F3B.auc <- plot_grid(F3Btop.auc + theme(legend.position = "none"), F3Bbot.auc, ncol = 1, align = "v", rel_heights = c(1,1),axis = "l")


#### S16B top - Randomizing thermal traits for rate and plotting slope above Topt ####

# Number of simulations
n_sim <- 1000

# Store all results in a list
simulated_slopes_rate <- vector("list", n_sim)

set.seed(123)  # For reproducibility
for (i in 1:n_sim) {
  shuffled_df <- ave_params_gcplyr_combined_major_taxa_rate
  shuffled_df$topt <- sample(shuffled_df$topt)
  shuffled_df$ctmax <- sample(shuffled_df$ctmax)
  shuffled_df$ctmin <- sample(shuffled_df$ctmin)  # optional if you're not using ctmin in slope
  shuffled_df$sim_id <- i
  
  # Calculate slope
  shuffled_df$superopt_slope <- -1 * (shuffled_df$rmax / (shuffled_df$ctmax - shuffled_df$topt))
  
  # Keep only the columns needed for plotting
  simulated_slopes_rate[[i]] <- shuffled_df[, c("rmax", "superopt_slope", "sim_id")]
}

# Combine into one dataframe
simulated_slopes_rate_df <- do.call(rbind, simulated_slopes_rate)

# Remove impossible slopes
simulated_slopes_rate_df <- simulated_slopes_rate_df %>% filter(superopt_slope <= 0) 

S16B.top <- ave_params_gcplyr_combined_major_taxa_rate %>%
  ggplot(aes(x=rmax, y= superopt_slope))+
  geom_point(data = simulated_slopes_rate_df, aes(x = rmax, y = superopt_slope),color = "grey85")+
  #scale_color_viridis_b()+
  geom_smooth(method = "lm", se = TRUE, fill = "grey60", color="grey90") +
  geom_point(alpha = 0.9,aes(color = Order, size = ctmax))+
  #facet_wrap(~`Climate profile`,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  ylim(-0.5,0)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  ylab(expression("Performance slope above Topt"))+
  xlab(expression("ratemax"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)
  )

#### S16B bot - Comparing slope ranks rate ####

slope_data <- ave_params_gcplyr_combined_major_taxa_rate %>%
  dplyr::select(c(superopt_slope,Order,rmax,seqID))%>%
  rename("observed_slope" = "superopt_slope")

slope_data_ranks <- ave_params_gcplyr_combined_major_taxa_rate %>%
  dplyr::select(c(superopt_slope,Order,rmax,seqID))%>%
  mutate("sim_id" = 0)

simulated_slopes_rate_df_3 <- simulated_slopes_rate_df %>% 
  left_join(slope_data,by = "rmax")%>%
  dplyr::select(!c(observed_slope))

slope_data_ranks_full <- bind_rows(slope_data_ranks,simulated_slopes_rate_df_3)

slope_data_ranks_full_2 <- slope_data_ranks_full %>%
  group_by(Order, rmax,seqID)%>%
  mutate("slope_rank" = rank(superopt_slope),
         "mean_slope" = mean(superopt_slope),
         "median_slope" = median(superopt_slope),
         "slope_ratio_mean" = superopt_slope/mean_slope,
         "slope_ratio_median" = superopt_slope/median_slope,
         sim_id = sim_id)


S16B.bot <- slope_data_ranks_full_2 %>% 
  filter(superopt_slope < 0) %>%
  filter(sim_id == 0) %>%
  ggplot(aes(x = rmax, y = slope_ratio_mean,color = Order)) +
  #geom_smooth()+
  ylim(0,2.5)+
  geom_hline(yintercept = 1,linetype = "dashed",linewidth = 0.5)+
  geom_point(size = 4,alpha = 0.8)+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  ylab(expression("Observed Slope / Average Slope"))+
  xlab(expression("ratemax"))+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)
  )

S16B <- plot_grid(S16B.top + theme(legend.position = "none"), S16.bot, ncol = 1, align = "v", rel_heights = c(1,1),axis = "l")


#### Calculating the linear correlation between the slope above Topt and peak performances ####

lm_aucmax_superopt_slope_auc <- lm(superopt_slope ~ rmax, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(lm_aucmax_superopt_slope_auc)
sigma(lm_aucmax_superopt_slope_auc)

lm_ratemax_superopt_slope_rate <- lm(superopt_slope ~ rmax, data = ave_params_gcplyr_combined_major_taxa_rate)
summary(lm_ratemax_superopt_slope_rate)
sigma(lm_ratemax_superopt_slope_rate)


lm_aucmax_superopt_slope_auc_shuffled <- lm(superopt_slope ~ rmax, data = simulated_slopes_auc_df)
summary(lm_aucmax_superopt_slope_auc_shuffled)
sigma(lm_aucmax_superopt_slope_auc_shuffled)

lm_ratemax_superopt_slope_rate_shuffled <- lm(superopt_slope ~ rmax, data = simulated_slopes_rate_df)
summary(lm_ratemax_superopt_slope_rate_shuffled)
sigma(lm_ratemax_superopt_slope_rate_shuffled)



# Extract slopes and standard errors
slope_auc <- coef(summary(lm_aucmax_superopt_slope_auc))["rmax", "Estimate"]
slope_auc_shuffled <- coef(summary(lm_aucmax_superopt_slope_auc_shuffled))["rmax", "Estimate"]
se_slope_auc<- coef(summary(lm_aucmax_superopt_slope_auc))["rmax", "Std. Error"]
se_slope_auc_shuffled <- coef(summary(lm_aucmax_superopt_slope_auc_shuffled))["rmax", "Std. Error"]

# Calculate Welch's T-test statistic
t_stat_Welch_auc <- (slope_auc - slope_auc_shuffled) / sqrt(se_slope_auc^2 + se_slope_auc_shuffled^2)

# Degrees of freedom
Welch_auc <- ((se_slope_auc^2 + se_slope_auc_shuffled^2)^2) / (((se_slope_auc^2)^2 / (266 - 2)) + ((se_slope_auc_shuffled^2)^2 / (133000 - 2)))

# Calculate p-value
p_value_Welch_auc <- 2 * pt(-abs(t_stat_Welch_auc), df)

# Display the test statistic and p-value
t_stat_Welch_auc
p_value_Welch_auc

# Extract slopes and standard errors
slope_rate <- coef(summary(lm_ratemax_superopt_slope_rate))["rmax", "Estimate"]
slope_rate_shuffled <- coef(summary(lm_ratemax_superopt_slope_rate_shuffled))["rmax", "Estimate"]
se_slope_rate<- coef(summary(lm_ratemax_superopt_slope_rate))["rmax", "Std. Error"]
se_slope_rate_shuffled <- coef(summary(lm_ratemax_superopt_slope_rate_shuffled))["rmax", "Std. Error"]

# Calculate Welch's T-test statistic
t_stat_Welch_rate <- (slope_rate - slope_rate_shuffled) / sqrt(se_slope_rate^2 + se_slope_rate_shuffled^2)

# Degrees of freedom
Welch_rate <- ((se_slope_rate^2 + se_slope_rate_shuffled^2)^2) / (((se_slope_rate^2)^2 / (266 - 2)) + ((se_slope_rate_shuffled^2)^2 / (133000 - 2)))

# Calculate p-value
p_value_Welch_rate <- 2 * pt(-abs(t_stat_Welch_rate), df)

# Display the test statistic and p-value
t_stat_Welch_rate
p_value_Welch_rate


#### Processing data to plot F4 by Climate Profile AUC ####
#
#
# In this figure we will see the TPCs overlapped with the temperature distributions of sites and 
# climate profiles, and we will calculate the thermal Fitness of the two main Orders in each
# environment
#
#

# Enterobacterales does not have representatives in site W, which is why
# we are excluding it from the analyses by Site

# Unifying temperature scales with curve df AUC 
# The temperature scales between trait data and temperature time series differ, and 
# thus we need to harmonize them

temperature_probability_df <- continuous_proportion_data_growing_season_site %>%
  rename("Temperature" = "Temp_continuous")



# Scale the datasets to a (0,1) range
continuous_proportion_data_growing_season_site <- continuous_proportion_data_growing_season_site %>%
  mutate(Proportion_scaled = (Proportion_continuous - min(Proportion_continuous)) /
           (max(Proportion_continuous) - min(Proportion_continuous)))

# Summarizing data by taxon, climate profile, and curve AUC 

# Summarize by curve
means_by_curve_climate_auc <- ave_preds_gcplyr_combined_major_taxa %>% filter(metric == "auc")

means_by_curve_climate_auc <- means_by_curve_climate_auc %>%
  #group_by(`Climate profile`) %>%
  mutate(Fitted_scaled = (.fitted - min(.fitted)) / (max(.fitted) - min(.fitted)))%>%
  ungroup()

means_by_curve_climate_auc <- means_by_curve_climate_auc %>%
  rename("Temperature" = "temp")


# We define the common temperature scale
common_temperature_auc <- seq(
  min(c(means_by_curve_climate_auc$Temperature, temperature_probability_df$Temperature)),
  max(c(means_by_curve_climate_auc$Temperature, temperature_probability_df$Temperature)),
  by = 0.1
)

# Interpolation function for a single group
interpolate_group <- function(data, common_temperature, cols_to_interp) {
  interpolated <- data.frame(Temperature = common_temperature)
  for (col in cols_to_interp) {
    interpolated[[col]] <- approx(
      data$Temperature, data[[col]], common_temperature_auc, rule = 2
    )$y
  }
  return(interpolated)
}

# We interpolate the temperature_probability_df for each  Site
probability_interpolated_auc <- temperature_probability_df %>%
  group_by(`Climate profile`,Site_initials) %>%
  group_split() %>%
  lapply(function(group) {
    interpolated <- interpolate_group(group, common_temperature_auc, c("Proportion_continuous", "Proportion_scaled"))
    interpolated$Site_initials <- unique(group$Site_initials)
    interpolated
  }) %>%
  bind_rows()



# We interpolate the means_by_curve_climate for each Climate profile
curve_means_interpolated_auc <- means_by_curve_climate_auc %>%
  group_by(seqID,Order, `Climate profile`,Site_initials) %>%
  group_split() %>%
  lapply(function(group) {
    interpolated <- interpolate_group(group, common_temperature_auc, c(".fitted", "Fitted_scaled"))
    interpolated$seqID <- unique(group$seqID)
    interpolated$Order <- unique(group$Order)
    interpolated$`Climate profile` <- unique(group$`Climate profile`)
    interpolated$Site_initials<- unique(group$Site_initials)
    interpolated
  }) %>%
  bind_rows()


probability_interpolated_auc <- probability_interpolated_auc %>%
  mutate(`Climate profile` = 
           case_when(Site_initials %in% c("W","S","SJ") ~ "Cold variable",
                     Site_initials %in% c("Y","B","V","P") ~ "Cold stable",
                     Site_initials %in% c("AB","DC","QR","PF","EC") ~ "Hot variable",
                     Site_initials %in% c("M","J","ME") ~ "Hot stable"))


curve_means_interpolated_auc <- curve_means_interpolated_auc %>% 
  rename("Original_climate" = "Climate profile")

probability_interpolated_auc <- probability_interpolated_auc %>% 
  rename("New_climate" = "Climate profile")


curve_means_interpolated_auc <- curve_means_interpolated_auc %>% 
  rename("Original_site" = "Site_initials")

probability_interpolated_auc <- probability_interpolated_auc %>% 
  rename("New_site" = "Site_initials")

# We finally merge the two interpolated dataframes by Temperature and Climate profile
Fitness_by_curve_climate_auc <- full_join(
  curve_means_interpolated_auc,
  probability_interpolated_auc,
  by = c("Temperature")) %>% 
  filter(Temperature >= 13 & Temperature <= 43)

Fitness_by_curve_climate_auc$Fitness <- Fitness_by_curve_climate_auc$Fitted_scaled * Fitness_by_curve_climate_auc$Proportion_scaled

Fitness_by_order_climate_auc <- Fitness_by_curve_climate_auc %>% 
  group_by(Order,Temperature,Original_climate,New_climate) %>%
  summarize(
    Fitness = mean(Fitness),
    Fitted_scaled = mean(Fitted_scaled)
  )%>%
  filter(Original_climate == New_climate)


# We integrate it by curve and Climate profile

Fitness_by_curve_climate_auc_summarized <- Fitness_by_curve_climate_auc %>%
  group_by(seqID,Order,Original_climate,New_climate,Original_site,New_site) %>%
  summarize("Integrated Fitness" = gcplyr::auc(Temperature,Fitness)) 


##### Calculating metrics of maladaptation by curve and Site AUC ####

Fitness_by_curve_climate_auc_summarized_Entero <- Fitness_by_curve_climate_auc_summarized %>% 
  filter(Order == "Enterobacterales") %>%
  ungroup() %>%  # ungroup AFTER filtering
  dplyr::select(-Order) %>%  # cleaner syntax to drop column
  mutate(across(where(is.numeric), ~ round(., 2)))


# We add the original_fitness as a column

original_fitness <- Fitness_by_curve_climate_auc_summarized_Entero %>%
  filter(Original_site == New_site) %>%
  rename(Original_fitness = `Integrated Fitness`)%>%
  dplyr::select(seqID, Original_fitness) 

Fitness_by_curve_climate_auc_summarized_Entero <- Fitness_by_curve_climate_auc_summarized_Entero %>%
  left_join(original_fitness, by = "seqID")


native_fitness <- Fitness_by_curve_climate_auc_summarized_Entero %>%
  filter(Original_site == New_site) %>%
  rename(Native_fitness = `Integrated Fitness`)%>%
  group_by(New_site) %>%
  summarize(sd_Native_fitness = sd(Native_fitness),
            Native_fitness = mean(Native_fitness))%>% 
  dplyr::select(New_site,Native_fitness) 

Fitness_by_curve_climate_auc_summarized_Entero <- Fitness_by_curve_climate_auc_summarized_Entero %>%
  left_join(native_fitness, by = "New_site")

Fitness_by_curve_climate_auc_summarized_Entero <- Fitness_by_curve_climate_auc_summarized_Entero %>%
  group_by(seqID,Original_site,New_site,Original_climate,New_climate)%>%
  summarize(delta_A_H = `Integrated Fitness` - Original_fitness,
            ratio_A_H = `Integrated Fitness` / Original_fitness,
            delta_F_L = `Integrated Fitness` - Native_fitness,
            ratio_F_L = `Integrated Fitness` / Native_fitness,
            `Integrated Fitness` = `Integrated Fitness`,
            Original_fitness = Original_fitness,
            Native_fitness = Native_fitness)%>%
  mutate(across(where(is.numeric), ~ round(., 2)))

Fitness_by_curve_climate_auc_summarized_Pseudo <- Fitness_by_curve_climate_auc_summarized %>% 
  filter(Order == "Pseudomonadales") %>%
  ungroup() %>%  # ungroup AFTER filtering
  dplyr::select(-Order) %>%  # cleaner syntax to drop column
  mutate(across(where(is.numeric), ~ round(., 2)))


# We add the original_fitness as a column

original_fitness <- Fitness_by_curve_climate_auc_summarized_Pseudo %>%
  filter(Original_site == New_site) %>%
  rename(Original_fitness = `Integrated Fitness`)%>%
  dplyr::select(seqID, Original_fitness) 

Fitness_by_curve_climate_auc_summarized_Pseudo <- Fitness_by_curve_climate_auc_summarized_Pseudo %>%
  left_join(original_fitness, by = "seqID")


native_fitness <- Fitness_by_curve_climate_auc_summarized_Pseudo %>%
  filter(Original_site == New_site) %>%
  rename(Native_fitness = `Integrated Fitness`)%>%
  group_by(New_site) %>%
  summarize(sd_Native_fitness = sd(Native_fitness),
            Native_fitness = mean(Native_fitness))%>% 
  dplyr::select(New_site,Native_fitness) 

Fitness_by_curve_climate_auc_summarized_Pseudo <- Fitness_by_curve_climate_auc_summarized_Pseudo %>%
  left_join(native_fitness, by = "New_site")

Fitness_by_curve_climate_auc_summarized_Pseudo <- Fitness_by_curve_climate_auc_summarized_Pseudo %>%
  group_by(seqID,Original_site,New_site,Original_climate,New_climate)%>%
  summarize(delta_A_H = `Integrated Fitness` - Original_fitness,
            ratio_A_H = `Integrated Fitness` / Original_fitness,
            delta_F_L = `Integrated Fitness` - Native_fitness,
            ratio_F_L = `Integrated Fitness` / Native_fitness,
            `Integrated Fitness` = `Integrated Fitness`,
            Original_fitness = Original_fitness,
            Native_fitness = Native_fitness)%>%
  mutate(across(where(is.numeric), ~ round(., 2)))





##### Summarizing metrics of maladaptation by Order and Site AUC ####

Fitness_by_order_climate_auc_summarized_Entero <- Fitness_by_curve_climate_auc_summarized_Entero %>%
  group_by(Original_climate,New_climate,Original_site,New_site)%>%
  summarize(sd_delta_A_H = sd(delta_A_H),
            sd_ratio_A_H = sd(ratio_A_H),
            sd_delta_F_L = sd(delta_F_L),
            sd_ratio_F_L = sd(ratio_F_L),
            delta_A_H = mean(delta_A_H),
            ratio_A_H = mean(ratio_A_H),
            delta_F_L = mean(delta_F_L),
            ratio_F_L = mean(ratio_F_L))%>%
  mutate(across(where(is.numeric), ~ round(., 2)))



Fitness_by_order_climate_auc_summarized_Pseudo <- Fitness_by_curve_climate_auc_summarized_Pseudo %>%
  group_by(Original_climate,New_climate,Original_site,New_site)%>%
  summarize(sd_delta_A_H = sd(delta_A_H),
            sd_ratio_A_H = sd(ratio_A_H),
            sd_delta_F_L = sd(delta_F_L),
            sd_ratio_F_L = sd(ratio_F_L),
            delta_A_H = mean(delta_A_H),
            ratio_A_H = mean(ratio_A_H),
            delta_F_L = mean(delta_F_L),
            ratio_F_L = mean(ratio_F_L))%>%
  mutate(across(where(is.numeric), ~ round(., 2)))




Fitness_by_order_climate_auc_summarized_Pseudo$New_site <- factor(Fitness_by_order_climate_auc_summarized_Pseudo$New_site,
                                                                  levels = c("Y","B","V","P","W","S","SJ",
                                                                             "M","J","ME","AB","DC","QR","PF","EC"))

Fitness_by_order_climate_auc_summarized_Pseudo$Original_site <- factor(Fitness_by_order_climate_auc_summarized_Pseudo$Original_site,
                                                                       levels = c("Y","B","V","P","W","S","SJ",
                                                                                  "M","J","ME","AB","DC","QR","PF","EC"))


Fitness_by_order_climate_auc_summarized_Entero$New_site <- factor(Fitness_by_order_climate_auc_summarized_Entero$New_site,
                                                                  levels = c("Y","B","V","P","W","S","SJ",
                                                                             "M","J","ME","AB","DC","QR","PF","EC"))

Fitness_by_order_climate_auc_summarized_Entero$Original_site <- factor(Fitness_by_order_climate_auc_summarized_Entero$Original_site,
                                                                       levels = c("Y","B","V","P","W","S","SJ",
                                                                                  "M","J","ME","AB","DC","QR","PF","EC"))


#### Plotting the temperature distribution by site, faceted by Climate profile, with the TPCs overlaid AUC ####

F4A1A.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Cold stable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Cold stable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = Site_initials),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Cold stable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = Site_initials),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


F4A1B.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Cold variable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Cold variable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = Site_initials),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Cold variable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = Site_initials),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


F4A1C.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Hot stable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Hot stable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = Site_initials),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Hot stable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = Site_initials),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

F4A1D.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Hot variable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Hot variable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = Site_initials),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Hot variable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = Site_initials),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 



F4A2B.auc <- Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Cold variable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("AUCmax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
        
  )

F4A2D.auc <- Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Hot variable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("AUCmax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
  )

F4A2A.auc <- Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Cold stable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("AUCmax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
        
        
  )

F4A2C.auc <- Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Hot stable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("AUCmax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
        
  )

F4AA.auc<- F4A1A.auc+ 
  annotation_custom(ggplotGrob(F4A2A.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AB.auc <- F4A1B.auc + 
  annotation_custom(ggplotGrob(F4A2B.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AC.auc <- F4A1C.auc + 
  annotation_custom(ggplotGrob(F4A2C.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AD.auc <- F4A1D.auc + 
  annotation_custom(ggplotGrob(F4A2D.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4A.auc.site <- plot_grid(F4AB.auc,F4AD.auc,F4AA.auc,F4AC.auc, 
                     ncol = 2,
                     rel_widths = c(1,1,1,1),
                     align = "hv",
                     axis = "l")
F4A.auc.site



#### F4 - Plotting the temperature distribution by Climate profile, with the TPCs overlaid AUC ####

# Scale the datasets to a (0,1) range
continuous_proportion_data_growing_season <- continuous_proportion_data_growing_season %>%
  mutate(Proportion_scaled = (Proportion_continuous - min(Proportion_continuous)) /
           (max(Proportion_continuous) - min(Proportion_continuous)))

F4A1A.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Cold stable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Cold stable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Cold stable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


F4A1B.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Cold variable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Cold variable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Cold variable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


F4A1C.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Hot stable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Hot stable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Hot stable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

F4A1D.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Hot variable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Hot variable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Hot variable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

F4AA.auc<- F4A1A.auc+ 
  annotation_custom(ggplotGrob(F4A2A.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AB.auc <- F4A1B.auc + 
  annotation_custom(ggplotGrob(F4A2B.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AC.auc <- F4A1C.auc + 
  annotation_custom(ggplotGrob(F4A2C.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AD.auc <- F4A1D.auc + 
  annotation_custom(ggplotGrob(F4A2D.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4A.auc <- plot_grid(F4AB.auc,F4AD.auc,F4AA.auc,F4AC.auc, 
                     ncol = 2,
                     rel_widths = c(1,1,1,1),
                     align = "hv",
                     axis = "l")
F4A.auc





#### S12 - Plotting the yearly temperature distribution by Climate profile, with the TPCs overlaid AUC ####
#
#
# In this figure we will see the TPCs overlapped with the YEARLY temperature distributions of sites and 
# climate profiles, and we will calculate the thermal Fitness of the two main Orders in each
# environment
#
#


# Unifying temperature scales with curve df AUC 
# The temperature scales between trait data and temperature time series differ, and 
# thus we need to harmonize them



# Scale the datasets to a (0,1) range
continuous_proportion_data <- continuous_proportion_data %>%
  mutate(Proportion_scaled = (Proportion_continuous - min(Proportion_continuous)) /
           (max(Proportion_continuous) - min(Proportion_continuous)))

yearly_temperature_probability_df <- continuous_proportion_data %>%
  rename("Temperature" = "Temp_continuous")

# Summarizing data by taxon, climate profile, and curve AUC 

# Summarize by curve
means_by_curve_climate_auc <- ave_preds_gcplyr_combined_major_taxa %>% filter(metric == "auc")

means_by_curve_climate_auc <- means_by_curve_climate_auc %>%
  #group_by(`Climate profile`) %>%
  mutate(Fitted_scaled = (.fitted - min(.fitted)) / (max(.fitted) - min(.fitted)))%>%
  ungroup()

means_by_curve_climate_auc <- means_by_curve_climate_auc %>%
  rename("Temperature" = "temp")


# We define the common temperature scale
yearly_common_temperature_auc <- seq(
  min(c(means_by_curve_climate_auc$Temperature, yearly_temperature_probability_df$Temperature)),
  max(c(means_by_curve_climate_auc$Temperature, yearly_temperature_probability_df$Temperature)),
  by = 0.1
)

# Interpolation function for a single group
interpolate_group <- function(data, common_temperature, cols_to_interp) {
  interpolated <- data.frame(Temperature = common_temperature)
  for (col in cols_to_interp) {
    interpolated[[col]] <- approx(
      data$Temperature, data[[col]], yearly_common_temperature_auc, rule = 2
    )$y
  }
  return(interpolated)
}

# We interpolate the temperature_probability_df for each Site
yearly_probability_interpolated_auc <- yearly_temperature_probability_df %>%
  group_by(`Climate profile`) %>%
  group_split() %>%
  lapply(function(group) {
    interpolated <- interpolate_group(group, yearly_common_temperature_auc, c("Proportion_continuous", "Proportion_scaled"))
    interpolated$`Climate profile` <- unique(group$`Climate profile`)
    interpolated
  }) %>%
  bind_rows()



# We interpolate the means_by_curve_climate for each Climate profile
yearly_curve_means_interpolated_auc <- means_by_curve_climate_auc %>%
  group_by(seqID,Order, `Climate profile`) %>%
  group_split() %>%
  lapply(function(group) {
    interpolated <- interpolate_group(group, yearly_common_temperature_auc, c(".fitted", "Fitted_scaled"))
    interpolated$seqID <- unique(group$seqID)
    interpolated$Order <- unique(group$Order)
    interpolated$`Climate profile` <- unique(group$`Climate profile`)
    interpolated
  }) %>%
  bind_rows()



yearly_curve_means_interpolated_auc <- yearly_curve_means_interpolated_auc %>% 
  rename("Original_climate" = "Climate profile")

yearly_probability_interpolated_auc <- yearly_probability_interpolated_auc %>% 
  rename("New_climate" = "Climate profile")


# We finally merge the two interpolated dataframes by Temperature and Climate profile
yearly_Fitness_by_curve_climate_auc <- full_join(
  yearly_curve_means_interpolated_auc,
  yearly_probability_interpolated_auc,
  by = c("Temperature")) %>% 
  filter(Temperature >= 13 & Temperature <= 43)

yearly_Fitness_by_curve_climate_auc$Fitness <- yearly_Fitness_by_curve_climate_auc$Fitted_scaled * yearly_Fitness_by_curve_climate_auc$Proportion_scaled

yearly_Fitness_by_order_climate_auc <- yearly_Fitness_by_curve_climate_auc %>% 
  group_by(Order,Temperature,Original_climate,New_climate) %>%
  summarize(
    Fitness = mean(Fitness),
    Fitted_scaled = mean(Fitted_scaled)
  )%>%
  filter(Original_climate == New_climate)


# We integrate it by curve and Climate profile

yearly_Fitness_by_curve_climate_auc_summarized <- yearly_Fitness_by_curve_climate_auc %>%
  group_by(seqID,Order,Original_climate,New_climate) %>%
  summarize("Integrated Fitness" = gcplyr::auc(Temperature,Fitness)) 


F4A2B.auc <- yearly_Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Cold variable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("AUCmax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
        
  )

F4A2D.auc <- yearly_Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Hot variable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("AUCmax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
  )

F4A2A.auc <- yearly_Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Cold stable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("AUCmax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
        
        
  )

F4A2C.auc <- yearly_Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Hot stable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("AUCmax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
        
  )


F4A1A.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = yearly_Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Cold stable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = yearly_Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = yearly_Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data %>% filter(`Climate profile` %in% c("Cold stable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data %>% filter(`Climate profile` %in% c("Cold stable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


F4A1B.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = yearly_Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Cold variable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = yearly_Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = yearly_Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Cold variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data %>% filter(`Climate profile` %in% c("Cold variable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data %>% filter(`Climate profile` %in% c("Cold variable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


F4A1C.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = yearly_Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Hot stable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = yearly_Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = yearly_Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data %>% filter(`Climate profile` %in% c("Hot stable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data %>% filter(`Climate profile` %in% c("Hot stable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

F4A1D.auc <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = yearly_Fitness_by_curve_climate_auc %>% filter(Original_climate %in% c("Hot variable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = yearly_Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = yearly_Fitness_by_order_climate_auc %>% filter(Original_climate %in% c("Hot variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data %>% filter(`Climate profile` %in% c("Hot variable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data %>% filter(`Climate profile` %in% c("Hot variable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

F4AA.auc<- F4A1A.auc+ 
  annotation_custom(ggplotGrob(F4A2A.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AB.auc <- F4A1B.auc + 
  annotation_custom(ggplotGrob(F4A2B.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AC.auc <- F4A1C.auc + 
  annotation_custom(ggplotGrob(F4A2C.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AD.auc <- F4A1D.auc + 
  annotation_custom(ggplotGrob(F4A2D.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

S12 <- plot_grid(F4AB.auc,F4AD.auc,F4AA.auc,F4AC.auc, 
                ncol = 2,
                rel_widths = c(1,1,1,1),
                align = "hv",
                axis = "l")
S12





#### Calculating the fraction of time bacteria are exposed to temperatures above Topt AUC ####

Fitness_by_order_climate_auc_topt <- Fitness_by_order_climate_auc %>%
  group_by(Original_climate,Order) %>%
  summarize(Topt = Temperature[Fitted_scaled == max(Fitted_scaled)],
            CTmax = Temperature[Fitted_scaled == min(Fitted_scaled [Temperature > Topt])])

growth_season_continuous_proportion_data_auc_topt <- continuous_proportion_data_growing_season %>%
  mutate("Original_climate" = continuous_proportion_data_growing_season$`Climate profile`)%>%
  left_join(Fitness_by_order_climate_auc_topt)

growth_season_supraoptimal_time_auc <- growth_season_continuous_proportion_data_auc_topt %>%
  group_by(Order, `Climate profile`) %>%
  summarize(
    Hard_times = (sum(Proportion_continuous[Temp_continuous > Topt])/sum(Proportion_continuous))*100,
    Harder_times = (sum(Proportion_continuous[Temp_continuous > CTmax])/sum(Proportion_continuous))*100,
    CTmax = first(CTmax),
    Topt = first(Topt)
  )

yearly_continuous_proportion_data_auc_topt <- continuous_proportion_data %>%
  mutate("Original_climate" = continuous_proportion_data$`Climate profile`)%>%
  left_join(Fitness_by_order_climate_auc_topt)

yearly_supraoptimal_time <- yearly_continuous_proportion_data_auc_topt %>%
  group_by(Order, `Climate profile`) %>%
  summarize(
    Hard_times = (sum(Proportion_continuous[Temp_continuous > Topt])/sum(Proportion_continuous))*100,
    Harder_times = (sum(Proportion_continuous[Temp_continuous > CTmax])/sum(Proportion_continuous))*100,
    CTmax = first(CTmax),
    Topt = first(Topt)
  )

#### Correlations between TPC traits AUC ####
#
#
# Here we generate a dataframe with the pairwise correlation indices obtained from 
# the global correlation analysis between TPC traits and climate variables
# We then try to plot the correlations using the simplified TPC triangles as template
#
#

# Combine upper and lower triangle data to get full triangle geometry
triangle_vertices_auc <- bind_rows(
  upper_triangle_data_noenv_Order,
  lower_triangle_data_noenv_Order
) %>% 
  filter(metric == "auc") %>%
  distinct()


arrow_data_auc <- triangle_vertices_auc %>%
  group_by(group) %>%
  summarise(
    Order = first(Order),
    metric = first(metric),
    segment_pairs = list({
      indices <- combn(n(), 2)  
      map_dfr(1:ncol(indices), function(k) {
        i <- indices[1, k]
        j <- indices[2, k]
        tibble(
          x_start = x[i],
          y_start = y[i],
          x_end   = x[j],
          y_end   = y[j]
        )
      })
    }),
    .groups = "drop"
  ) %>%
  unnest(segment_pairs)%>%
  mutate("arrow" = c("Topt_CTmax",
                     "Max_CTmax",
                     "CTmin_CTmax",
                     "Topt_Max",
                     "Topt_CTmin",
                     "Max_CTmin",
                     "Topt_CTmax",
                     "Max_CTmax",
                     "CTmin_CTmax",
                     "Topt_Max",
                     "Topt_CTmin",
                     "Max_CTmin"))


Correlation_between_thermal_traits_auc <- data.frame("arrow" = c("Topt_CTmax",
                                                             "Max_CTmax",
                                                             "CTmin_CTmax",
                                                             "Topt_Max",
                                                             "Topt_CTmin",
                                                             "Max_CTmin",
                                                             "Topt_CTmax",
                                                             "Max_CTmax",
                                                             "CTmin_CTmax",
                                                             "Topt_Max",
                                                             "Topt_CTmin",
                                                             "Max_CTmin"),
                                                 "correlation" = c(0.46,0.17,-0.05,0.28,0.11,0.41,
                                                                   0.27,0.1,-0.44,-0.01,0.52,0.03),
                                                 "Order" = c(rep("Enterobacterales",6),
                                                             rep("Pseudomonadales",6)))


arrow_data_auc <- left_join(arrow_data_auc,Correlation_between_thermal_traits_auc,by = c("Order","arrow"))%>%
  filter(!arrow %in% c("CTmin_CTmax"))

arrow_data_auc <- arrow_data_auc %>%
  mutate(
    x_mid = (x_start + x_end) / 2,
    y_mid = (y_start + y_end) / 2
  )

F3A.auc.correlations <- full_triangle_data_noenv %>% 
  filter(metric == "auc") %>%
  ggplot(aes(x = x, y = y, group = group)) +
  geom_vline(xintercept = 42, linetype = "dashed", color = "grey60", linewidth = 0.5)+
  #geom_polygon(data = upper_triangle_data_noenv_Order %>% 
  #               filter(metric == "auc"),
  #             aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  #geom_polygon(data = lower_triangle_data_noenv_Order %>% 
  #               filter(metric == "auc"),
  #             aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  geom_segment(data = arrow_data_auc,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = correlation),
               arrow = arrow(length = unit(0, "inches"),
                             type = "closed"),
               inherit.aes = FALSE,
               linewidth = 3)+
  geom_label(
    data = arrow_data_auc,
    aes(x = x_mid, y = y_mid, label = round(correlation, 2),color = correlation,fill = correlation),
    size = 5,
    label.size = 0,
    inherit.aes = FALSE,
    family = "Atkinson Hyperlegible Next VF Light")+
  geom_label(
    data = arrow_data_auc,
    aes(x = x_mid, y = y_mid, label = round(correlation, 2)),
    size = 3.5,
    fill = "white",
    color = "black",
    label.size = 0.1,
    inherit.aes = FALSE,
    family = "Atkinson Hyperlegible Next VF Light")+
  geom_point(data = full_triangle_data_noenv_Order %>% 
               filter(metric == "auc"),
             shape = 21,fill = "black",color = "white",aes(x = x, y = y,size = z))+
  labs(
    x = "Temperature (°C)",
    y = "AUC",
    color = "Pearson’s r",
    size = "SD/Mean"
  ) +
  facet_grid(Order~.)+
  scale_color_continuous_diverging(rev = TRUE)+
  scale_fill_continuous_diverging(rev = TRUE)+
  #scale_fill_viridis_d(option = "turbo")+
  #scale_fill_manual(values = c("#DE5925","#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        #legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(fill = "none")

#### Processing data to plot F4 by Climate Profile rate ####
#
#
# In this figure we will see the TPCs overlapped with the temperature distributions of sites and 
# climate profiles, and we will calculate the thermal Fitness of the two main Orders in each
# environment
#
#

# Enterobacterales does not have representatives in site W, which is why
# we are excluding it from the analyses by Site

# Scale the datasets to a (0,1) range
continuous_proportion_data_growing_season_site <- continuous_proportion_data_growing_season_site %>%
  mutate(Proportion_scaled = (Proportion_continuous - min(Proportion_continuous)) /
           (max(Proportion_continuous) - min(Proportion_continuous)))

#### Summarizing data by taxon, climate profile, and curve rate ####

# Summarize by curve
means_by_curve_climate_rate <- ave_preds_gcplyr_combined_major_taxa %>% filter(metric == "rate")

means_by_curve_climate_rate <- means_by_curve_climate_rate %>%
  #group_by(`Climate profile`) %>%
  mutate(Fitted_scaled = (.fitted - min(.fitted)) / (max(.fitted) - min(.fitted)))%>%
  ungroup()


#### Unifying temperature scales with curve df rate ######
# The temperature scales between trait data and temperature time series differ, and 
# thus we need to harmonize them

temperature_probability_df <- continuous_proportion_data_growing_season_site %>%
  rename("Temperature" = "Temp_continuous")

# We define the common temperature scale
common_temperature_rate <- seq(
  min(c(means_by_curve_climate_rate$Temperature, temperature_probability_df$Temperature)),
  max(c(means_by_curve_climate_rate$Temperature, temperature_probability_df$Temperature)),
  by = 0.1
)

# Interpolation function for a single group
interpolate_group <- function(data, common_temperature, cols_to_interp) {
  interpolated <- data.frame(Temperature = common_temperature)
  for (col in cols_to_interp) {
    interpolated[[col]] <- approx(
      data$Temperature, data[[col]], common_temperature_rate, rule = 2
    )$y
  }
  return(interpolated)
}

# We interpolate the temperature_probability_df for each  Site
probability_interpolated_rate <- temperature_probability_df %>%
  group_by(`Climate profile`,Site_initials) %>%
  group_split() %>%
  lapply(function(group) {
    interpolated <- interpolate_group(group, common_temperature_rate, c("Proportion_continuous", "Proportion_scaled"))
    interpolated$Site_initials <- unique(group$Site_initials)
    interpolated
  }) %>%
  bind_rows()


# We interpolate the means_by_curve_climate for each Climate profile
curve_means_interpolated_rate <- means_by_curve_climate_rate %>%
  group_by(seqID,Order, `Climate profile`,Site_initials) %>%
  group_split() %>%
  lapply(function(group) {
    interpolated <- interpolate_group(group, common_temperature_rate, c(".fitted", "Fitted_scaled"))
    interpolated$seqID <- unique(group$seqID)
    interpolated$Order <- unique(group$Order)
    interpolated$`Climate profile` <- unique(group$`Climate profile`)
    interpolated$Site_initials<- unique(group$Site_initials)
    interpolated
  }) %>%
  bind_rows()


probability_interpolated_rate <- probability_interpolated_rate %>%
  mutate(`Climate profile` = 
           case_when(Site_initials %in% c("W","S","SJ") ~ "Cold variable",
                     Site_initials %in% c("Y","B","V","P") ~ "Cold stable",
                     Site_initials %in% c("AB","DC","QR","PF") ~ "Hot variable",
                     Site_initials %in% c("EC","M","J","ME") ~ "Hot stable"))


curve_means_interpolated_rate <- curve_means_interpolated_rate %>% 
  rename("Original_climate" = "Climate profile")

probability_interpolated_rate <- probability_interpolated_rate %>% 
  rename("New_climate" = "Climate profile")


curve_means_interpolated_rate <- curve_means_interpolated_rate %>% 
  rename("Original_site" = "Site_initials")

probability_interpolated_rate <- probability_interpolated_rate %>% 
  rename("New_site" = "Site_initials")

# We finally merge the two interpolated dataframes by Temperature and Climate profile
Fitness_by_curve_climate_rate <- full_join(
  curve_means_interpolated_rate,
  probability_interpolated_rate,
  by = c("Temperature")) %>% 
  filter(Temperature >= 13 & Temperature <= 43)

Fitness_by_curve_climate_rate$Fitness <- Fitness_by_curve_climate_rate$Fitted_scaled * Fitness_by_curve_climate_rate$Proportion_scaled

Fitness_by_order_climate_rate <- Fitness_by_curve_climate_rate %>% 
  group_by(Order,Temperature,Original_climate,New_climate) %>%
  summarize(
    Fitness = mean(Fitness),
    Fitted_scaled = mean(Fitted_scaled)
  )%>%
  filter(Original_climate == New_climate)


# We integrate it by curve and Climate profile

Fitness_by_curve_climate_rate_summarized <- Fitness_by_curve_climate_rate %>%
  group_by(seqID,Order,Original_climate,New_climate,Original_site,New_site) %>%
  summarize("Integrated Fitness" = gcplyr::auc(Temperature,Fitness)) 


#### Calculating metrics of maladaptation by curve and Site rate ####

Fitness_by_curve_climate_rate_summarized_Entero <- Fitness_by_curve_climate_rate_summarized %>% 
  filter(Order == "Enterobacterales") %>%
  ungroup() %>%  # ungroup AFTER filtering
  dplyr::select(-Order) %>%  # cleaner syntax to drop column
  mutate(across(where(is.numeric), ~ round(., 2)))


# We add the original_fitness as a column

original_fitness <- Fitness_by_curve_climate_rate_summarized_Entero %>%
  filter(Original_site == New_site) %>%
  rename(Original_fitness = `Integrated Fitness`)%>%
  dplyr::select(seqID, Original_fitness) 

Fitness_by_curve_climate_rate_summarized_Entero <- Fitness_by_curve_climate_rate_summarized_Entero %>%
  left_join(original_fitness, by = "seqID")


native_fitness <- Fitness_by_curve_climate_rate_summarized_Entero %>%
  filter(Original_site == New_site) %>%
  rename(Native_fitness = `Integrated Fitness`)%>%
  group_by(New_site) %>%
  summarize(sd_Native_fitness = sd(Native_fitness),
            Native_fitness = mean(Native_fitness))%>% 
  dplyr::select(New_site,Native_fitness) 

Fitness_by_curve_climate_rate_summarized_Entero <- Fitness_by_curve_climate_rate_summarized_Entero %>%
  left_join(native_fitness, by = "New_site")

Fitness_by_curve_climate_rate_summarized_Entero <- Fitness_by_curve_climate_rate_summarized_Entero %>%
  group_by(seqID,Original_site,New_site,Original_climate,New_climate)%>%
  summarize(delta_A_H = `Integrated Fitness` - Original_fitness,
            ratio_A_H = `Integrated Fitness` / Original_fitness,
            delta_F_L = `Integrated Fitness` - Native_fitness,
            ratio_F_L = `Integrated Fitness` / Native_fitness,
            `Integrated Fitness` = `Integrated Fitness`,
            Original_fitness = Original_fitness,
            Native_fitness = Native_fitness)%>%
  mutate(across(where(is.numeric), ~ round(., 2)))

Fitness_by_curve_climate_rate_summarized_Pseudo <- Fitness_by_curve_climate_rate_summarized %>% 
  filter(Order == "Pseudomonadales") %>%
  ungroup() %>%  # ungroup AFTER filtering
  dplyr::select(-Order) %>%  # cleaner syntax to drop column
  mutate(across(where(is.numeric), ~ round(., 2)))


# We add the original_fitness as a column

original_fitness <- Fitness_by_curve_climate_rate_summarized_Pseudo %>%
  filter(Original_site == New_site) %>%
  rename(Original_fitness = `Integrated Fitness`)%>%
  dplyr::select(seqID, Original_fitness) 

Fitness_by_curve_climate_rate_summarized_Pseudo <- Fitness_by_curve_climate_rate_summarized_Pseudo %>%
  left_join(original_fitness, by = "seqID")


native_fitness <- Fitness_by_curve_climate_rate_summarized_Pseudo %>%
  filter(Original_site == New_site) %>%
  rename(Native_fitness = `Integrated Fitness`)%>%
  group_by(New_site) %>%
  summarize(sd_Native_fitness = sd(Native_fitness),
            Native_fitness = mean(Native_fitness))%>% 
  dplyr::select(New_site,Native_fitness) 

Fitness_by_curve_climate_rate_summarized_Pseudo <- Fitness_by_curve_climate_rate_summarized_Pseudo %>%
  left_join(native_fitness, by = "New_site")

Fitness_by_curve_climate_rate_summarized_Pseudo <- Fitness_by_curve_climate_rate_summarized_Pseudo %>%
  group_by(seqID,Original_site,New_site,Original_climate,New_climate)%>%
  summarize(delta_A_H = `Integrated Fitness` - Original_fitness,
            ratio_A_H = `Integrated Fitness` / Original_fitness,
            delta_F_L = `Integrated Fitness` - Native_fitness,
            ratio_F_L = `Integrated Fitness` / Native_fitness,
            `Integrated Fitness` = `Integrated Fitness`,
            Original_fitness = Original_fitness,
            Native_fitness = Native_fitness)%>%
  mutate(across(where(is.numeric), ~ round(., 2)))





#### Summarizing metrics of maladaptation by Order and Site rate ####

Fitness_by_order_climate_rate_summarized_Entero <- Fitness_by_curve_climate_rate_summarized_Entero %>%
  group_by(Original_climate,New_climate,Original_site,New_site)%>%
  summarize(sd_delta_A_H = sd(delta_A_H),
            sd_ratio_A_H = sd(ratio_A_H),
            sd_delta_F_L = sd(delta_F_L),
            sd_ratio_F_L = sd(ratio_F_L),
            delta_A_H = mean(delta_A_H),
            ratio_A_H = mean(ratio_A_H),
            delta_F_L = mean(delta_F_L),
            ratio_F_L = mean(ratio_F_L))%>%
  mutate(across(where(is.numeric), ~ round(., 2)))



Fitness_by_order_climate_rate_summarized_Pseudo <- Fitness_by_curve_climate_rate_summarized_Pseudo %>%
  group_by(Original_climate,New_climate,Original_site,New_site)%>%
  summarize(sd_delta_A_H = sd(delta_A_H),
            sd_ratio_A_H = sd(ratio_A_H),
            sd_delta_F_L = sd(delta_F_L),
            sd_ratio_F_L = sd(ratio_F_L),
            delta_A_H = mean(delta_A_H),
            ratio_A_H = mean(ratio_A_H),
            delta_F_L = mean(delta_F_L),
            ratio_F_L = mean(ratio_F_L))%>%
  mutate(across(where(is.numeric), ~ round(., 2)))




Fitness_by_order_climate_rate_summarized_Pseudo$New_site <- factor(Fitness_by_order_climate_rate_summarized_Pseudo$New_site,
                                                                   levels = c("Y","B","V","P","W","S","SJ",
                                                                              "M","J","ME","AB","DC","QR","PF","EC"))
Fitness_by_order_climate_rate_summarized_Pseudo$Original_site <- factor(Fitness_by_order_climate_rate_summarized_Pseudo$Original_site,
                                                                        levels = c("Y","B","V","P","W","S","SJ",
                                                                                   "M","J","ME","AB","DC","QR","PF","EC"))


Fitness_by_order_climate_rate_summarized_Entero$New_site <- factor(Fitness_by_order_climate_rate_summarized_Entero$New_site,
                                                                   levels = c("Y","B","V","P","W","S","SJ",
                                                                              "M","J","ME","AB","DC","QR","PF","EC"))

Fitness_by_order_climate_rate_summarized_Entero$Original_site <- factor(Fitness_by_order_climate_rate_summarized_Entero$Original_site,
                                                                        levels = c("Y","B","V","P","W","S","SJ",
                                                                                   "M","J","ME","AB","DC","QR","PF","EC"))

#### Plotting the temperature profile and scaled TPCs, along with the integrated fitness for the two majoritary Orders by Site rate ####


F4_by_site_rate <- Fitness_by_curve_climate_rate %>%
  filter (Original_site == New_site)%>% rename(`Climate profile` = Original_climate,
                                               Site_initials = Original_site) %>%
  ggplot()  +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  facet_wrap(`Climate profile`~Site_initials,nrow = 4)+
  geom_line(lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_line(
    aes(x = Temperature, y = Proportion_scaled, group = Site_initials),
    size = 0.5
  ) +
  geom_ribbon(
    aes(x = Temperature, ymin = 0, ymax = Proportion_scaled, group = Site_initials),
    alpha = 0.6,fill = "grey85"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

#### Plotting the temperature distribution by site, faceted by Climate profile, with the TPCs overlaid rate ####

F4A1A.rate <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_rate %>% filter(Original_climate %in% c("Cold stable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Cold stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Cold stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Cold stable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = Site_initials),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Cold stable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = Site_initials),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


F4A1B.rate <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_rate %>% filter(Original_climate %in% c("Cold variable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Cold variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Cold variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Cold variable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = Site_initials),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Cold variable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = Site_initials),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


F4A1C.rate <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_rate %>% filter(Original_climate %in% c("Hot stable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Hot stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Hot stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Hot stable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = Site_initials),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Hot stable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = Site_initials),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

F4A1D.rate <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_rate %>% filter(Original_climate %in% c("Hot variable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Hot variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Hot variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Hot variable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = Site_initials),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season_site %>% filter(`Climate profile` %in% c("Hot variable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = Site_initials),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 



F4A2B.rate <- Fitness_by_curve_climate_rate_summarized %>% 
  filter(Original_climate %in% c("Cold variable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("ratemax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
        
  )

F4A2D.rate <- Fitness_by_curve_climate_rate_summarized %>% 
  filter(Original_climate %in% c("Hot variable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("ratemax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
  )

F4A2A.rate <- Fitness_by_curve_climate_rate_summarized %>% 
  filter(Original_climate %in% c("Cold stable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("ratemax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
        
        
  )

F4A2C.rate <- Fitness_by_curve_climate_rate_summarized %>% 
  filter(Original_climate %in% c("Hot stable")) %>%
  ggplot(aes(x=Order, y=`Integrated Fitness`))+
  #scale_color_viridis_b()+
  #geom_violin(alpha = 0.4,outliers = FALSE)+
  geom_jitter(size = 2,width = 0.1,aes(color = Order))+
  #geom_dotplot(binaxis= "y",stackdir = "center",dotsize = 1, aes(color = Order, fill = Order))+ 
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  #facet_wrap(~Original_climate,ncol = 2)+
  #geom_abline(slope = -0.02,intercept = -0.1)+
  #ylim(16,40)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression("Performance slope above Topt"))+
  #xlab(expression("ratemax"))+
  theme_classic(base_size = 12,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill = "transparent", color = NA)
        
  )

F4AA.rate<- F4A1A.rate+ 
  annotation_custom(ggplotGrob(F4A2A.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AB.rate <- F4A1B.rate + 
  annotation_custom(ggplotGrob(F4A2B.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AC.rate <- F4A1C.rate + 
  annotation_custom(ggplotGrob(F4A2C.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AD.rate <- F4A1D.rate + 
  annotation_custom(ggplotGrob(F4A2D.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4A.rate.site <- plot_grid(F4AB.rate,F4AD.rate,F4AA.rate,F4AC.rate, 
                     ncol = 2,
                     rel_widths = c(1,1,1,1),
                     align = "hv",
                     axis = "l")
F4A.rate.site



#### S17 - Plotting the temperature distribution by Climate profile, with the TPCs overlaid rate ####

F4A1A.rate <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_rate %>% filter(Original_climate %in% c("Cold stable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Cold stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Cold stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Cold stable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Cold stable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


F4A1B.rate <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_rate %>% filter(Original_climate %in% c("Cold variable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Cold variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Cold variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Cold variable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Cold variable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


F4A1C.rate <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_rate %>% filter(Original_climate %in% c("Hot stable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Hot stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Hot stable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Hot stable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Hot stable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

F4A1D.rate <- ggplot() +
  annotate(
    "rect",
    xmin = -5, xmax = 55,
    ymin = 0, ymax = 1.25,
    color = "black", fill = "white",
    linewidth = 1)+
  geom_line(data = Fitness_by_curve_climate_rate %>% filter(Original_climate %in% c("Hot variable")),
            lwd = 0.5,alpha = 0.5,aes(x = Temperature, y = Fitted_scaled,group = seqID,color = Order)) +
  geom_area(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Hot variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,fill = Order),
    alpha = 0.2,position = "identity",
  ) +
  geom_line(
    data = Fitness_by_order_climate_rate %>% filter(Original_climate %in% c("Hot variable")),
    aes(x = Temperature, y = Fitted_scaled, group = Order,color = Order),
    alpha = 1,position = "identity",linewidth = 1.5
  ) +
  geom_line(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Hot variable")),
    aes(x = Temp_continuous, y = Proportion_scaled, group = `Climate profile`),
    size = 0.5
  ) +
  geom_ribbon(
    data = continuous_proportion_data_growing_season %>% filter(`Climate profile` %in% c("Hot variable")),
    aes(x = Temp_continuous, ymin = 0, ymax = Proportion_scaled, group = `Climate profile`),
    alpha = 0.3,fill = "grey95"
  ) +
  scale_y_continuous(
    name = "Temperature probability",
    sec.axis = sec_axis(~ ., name = "Scaled Performance")
  ) +
  #scale_color_viridis_d(option = "turbo", direction = 1) +
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  #scale_fill_viridis_d(option = "turbo", direction = 1) +
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  xlim(-5,55)+
  ylim(0,1.1)+
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = NA),
    #axis.title = element_blank(),
    ###axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

F4AA.rate<- F4A1A.rate+ 
  annotation_custom(ggplotGrob(F4A2A.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AB.rate <- F4A1B.rate + 
  annotation_custom(ggplotGrob(F4A2B.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AC.rate <- F4A1C.rate + 
  annotation_custom(ggplotGrob(F4A2C.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F4AD.rate <- F4A1D.rate + 
  annotation_custom(ggplotGrob(F4A2D.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

S17 <- plot_grid(F4AB.rate,F4AD.rate,F4AA.rate,F4AC.rate, 
                     ncol = 2,
                     rel_widths = c(1,1,1,1),
                     align = "hv",
                     axis = "l")





#### Calculating the fraction of time bacteria are exposed to temperatures above Topt rate ####

Fitness_by_order_climate_rate_topt <- Fitness_by_order_climate_rate %>%
  group_by(Original_climate,Order) %>%
  summarize(Topt = Temperature[Fitted_scaled == max(Fitted_scaled)],
            CTmax = Temperature[Fitted_scaled == min(Fitted_scaled [Temperature > Topt])])

growth_season_continuous_proportion_data_rate_topt <- continuous_proportion_data_growing_season %>%
  mutate("Original_climate" = continuous_proportion_data_growing_season$`Climate profile`)%>%
  left_join(Fitness_by_order_climate_rate_topt)

growth_season_supraoptimal_time_rate <- growth_season_continuous_proportion_data_rate_topt %>%
  group_by(Order, `Climate profile`) %>%
  summarize(
    Hard_times = (sum(Proportion_continuous[Temp_continuous > Topt])/sum(Proportion_continuous))*100,
    Harder_times = (sum(Proportion_continuous[Temp_continuous > CTmax])/sum(Proportion_continuous))*100,
    CTmax = first(CTmax),
    Topt = first(Topt)
  )

yearly_continuous_proportion_data_rate_topt <- continuous_proportion_data %>%
  mutate("Original_climate" = continuous_proportion_data$`Climate profile`)%>%
  left_join(Fitness_by_order_climate_rate_topt)

yearly_supraoptimal_time <- yearly_continuous_proportion_data_rate_topt %>%
  group_by(Order, `Climate profile`) %>%
  summarize(
    Hard_times = (sum(Proportion_continuous[Temp_continuous > Topt])/sum(Proportion_continuous))*100,
    Harder_times = (sum(Proportion_continuous[Temp_continuous > CTmax])/sum(Proportion_continuous))*100,
    CTmax = first(CTmax),
    Topt = first(Topt)
  )

#### S16A - Correlations between TPC traits rate ####
#
#
# Here we generate a dataframe with the pairwise correlation indices obtained from 
# the global correlation analysis between TPC traits and climate variables
# We then try to plot the correlations using the simplified TPC triangles as template
#
#

# Combine upper and lower triangle data to get full triangle geometry
triangle_vertices_rate <- bind_rows(
  upper_triangle_data_noenv_Order,
  lower_triangle_data_noenv_Order
) %>% 
  filter(metric == "rate") %>%
  distinct()


arrow_data_rate <- triangle_vertices_rate %>%
  group_by(group) %>%
  summarise(
    Order = first(Order),
    metric = first(metric),
    segment_pairs = list({
      indices <- combn(n(), 2)  
      map_dfr(1:ncol(indices), function(k) {
        i <- indices[1, k]
        j <- indices[2, k]
        tibble(
          x_start = x[i],
          y_start = y[i],
          x_end   = x[j],
          y_end   = y[j]
        )
      })
    }),
    .groups = "drop"
  ) %>%
  unnest(segment_pairs)%>%
  mutate("arrow" = c("Topt_CTmax",
                     "Max_CTmax",
                     "CTmin_CTmax",
                     "Topt_Max",
                     "Topt_CTmin",
                     "Max_CTmin",
                     "Topt_CTmax",
                     "Max_CTmax",
                     "CTmin_CTmax",
                     "Topt_Max",
                     "Topt_CTmin",
                     "Max_CTmin"))


Correlation_between_thermal_traits <- data.frame("arrow" = c("Topt_CTmax",
                                                             "Max_CTmax",
                                                             "CTmin_CTmax",
                                                             "Topt_Max",
                                                             "Topt_CTmin",
                                                             "Max_CTmin",
                                                             "Topt_CTmax",
                                                             "Max_CTmax",
                                                             "CTmin_CTmax",
                                                             "Topt_Max",
                                                             "Topt_CTmin",
                                                             "Max_CTmin"),
                                                 "correlation" = c(0.45,-0.3,-0.61,0.09,-0.12,0.46,
                                                                   0.41,-0.1,-0.46,-0.05,0.41,0.11),
                                                 "Order" = c(rep("Enterobacterales",6),
                                                             rep("Pseudomonadales",6)))


arrow_data_rate <- left_join(arrow_data_rate,Correlation_between_thermal_traits,by = c("Order","arrow"))%>%
  filter(!arrow %in% c("CTmin_CTmax"))

arrow_data_rate <- arrow_data_rate %>%
  mutate(
    x_mid = (x_start + x_end) / 2,
    y_mid = (y_start + y_end) / 2
  )

S16A <- full_triangle_data_noenv %>% 
  filter(metric == "rate") %>%
  ggplot(aes(x = x, y = y, group = group)) +
  geom_vline(xintercept = 42, linetype = "dashed", color = "grey60", linewidth = 0.5)+
  #geom_polygon(data = upper_triangle_data_noenv_Order %>% 
  #               filter(metric == "rate"),
  #             aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  #geom_polygon(data = lower_triangle_data_noenv_Order %>% 
  #               filter(metric == "rate"),
  #             aes(x = x, y = y, group = group,fill = Order),color = "white", alpha = 1, inherit.aes = FALSE) +
  geom_segment(data = arrow_data_rate,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = correlation),
               arrow = arrow(length = unit(0, "inches"),
                             type = "closed"),
               inherit.aes = FALSE,
               linewidth = 3)+
  geom_label(
    data = arrow_data_rate,
    aes(x = x_mid, y = y_mid, label = round(correlation, 2),color = correlation,fill = correlation),
    size = 5,
    label.size = 0,
    inherit.aes = FALSE,
    family = "Atkinson Hyperlegible Next VF Light")+
  geom_label(
    data = arrow_data_rate,
    aes(x = x_mid, y = y_mid, label = round(correlation, 2)),
    size = 3.5,
    fill = "white",
    color = "black",
    label.size = 0.1,
    inherit.aes = FALSE,
    family = "Atkinson Hyperlegible Next VF Light")+
  geom_point(data = full_triangle_data_noenv_Order %>% 
               filter(metric == "rate"),
             shape = 21,fill = "black",color = "white",aes(x = x, y = y,size = z))+
  labs(
    x = "Temperature (°C)",
    y = "rate",
    color = "Pearson’s r",
    size = "SD/Mean"
  ) +
  facet_grid(Order~.)+
  scale_color_continuous_diverging(rev = TRUE)+
  scale_fill_continuous_diverging(rev = TRUE)+
  #scale_fill_viridis_d(option = "turbo")+
  #scale_fill_manual(values = c("#DE5925","#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        #legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )+
  guides(fill = "none")




#### S10 - Coloring slope above Topt by climate profile ####

S10 <- ave_params_gcplyr_combined_major_taxa %>%
  filter(metric == "auc")%>%
  ggplot(aes(x=rmax, y=superopt_slope))+
  #scale_color_viridis_b()+
  #geom_smooth(method = "lm", se = TRUE, fill = "grey40", color="grey90") +
  geom_abline(slope = -0.065,intercept = c(0,0))+
  geom_smooth(method = "lm", se = TRUE, fill = "grey40", aes(color = `Climate profile`)) +
  geom_point(alpha = 0.8,aes(color = `Climate profile`, size = ctmax))+
  facet_wrap(~`Climate profile`,ncol = 2)+
  #ylim(-0.4,0)+
  #xlim(0.75,3.15)+
  #scale_fill_manual(values = order_color)+
  ylab(expression("Performance slope above Topt"))+
  xlab(expression("AUCmax"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "plasma")+
  scale_color_manual(values = palette_colors_Spectral_Climate)+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        #legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)
  )


#### S11 - Calculating the variance of Topt and CTmax for bins of AUC ####

ave_params_gcplyr_combined_major_taxa_auc_Entero <- 
  ave_params_gcplyr_combined_major_taxa_auc %>%
  filter(Order == "Enterobacterales") %>%
  mutate("Growth quartile" = case_when(
    rmax <= quantile(rmax)[2] ~ "q1",
    rmax <= quantile(rmax)[3] & rmax > quantile(rmax)[2] ~ "q2",
    rmax <= quantile(rmax)[4] & rmax > quantile(rmax)[3] ~ "q3",
    rmax <= quantile(rmax)[5] & rmax > quantile(rmax)[4] ~ "q4",
  ))

ave_params_gcplyr_combined_major_taxa_auc_Pseudo <- 
  ave_params_gcplyr_combined_major_taxa_auc %>%
  filter(Order == "Pseudomonadales")%>%
  mutate("Growth quartile" = case_when(
    rmax <= quantile(rmax)[2] ~ "q1",
    rmax <= quantile(rmax)[3] & rmax > quantile(rmax)[2] ~ "q2",
    rmax <= quantile(rmax)[4] & rmax > quantile(rmax)[3] ~ "q3",
    rmax <= quantile(rmax)[5] & rmax > quantile(rmax)[4] ~ "q4",
  ))

ave_params_gcplyr_combined_major_taxa_binned <- 
  bind_rows(ave_params_gcplyr_combined_major_taxa_auc_Entero,
            ave_params_gcplyr_combined_major_taxa_auc_Pseudo)


f114.a <- ave_params_gcplyr_combined_major_taxa_binned %>%
  ggplot(aes(x = `Growth quartile`, y = topt,color = Order))+
  geom_boxplot(alpha = 0.6, outliers = FALSE, aes(group = `Growth quartile`))+
  geom_point()+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  facet_grid(~Order)+
  #scale_fill_gradient2(low="#DE5925",mid = "white", high="#0090B5")+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    #strip.text.x.top = ,
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

f114.b <- ave_params_gcplyr_combined_major_taxa_binned %>%
  ggplot(aes(x = `Growth quartile`, y = ctmax,color = Order))+
  geom_boxplot(alpha = 0.6, outliers = FALSE, aes(group = `Growth quartile`))+
  geom_point()+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  facet_grid(~Order)+
  #scale_fill_gradient2(low="#DE5925",mid = "white", high="#0090B5")+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    #strip.text.x.top = ,
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


ave_params_gcplyr_combined_major_taxa_binned_variance <- ave_params_gcplyr_combined_major_taxa_binned %>%
  group_by(Order, `Growth quartile`) %>%
  summarize(sd_topt = sd(topt),
            sd_ctmax = sd(ctmax),
            sd_ctmin = sd(ctmin))


f114.c <- ave_params_gcplyr_combined_major_taxa_binned_variance %>%
  ggplot(aes(x = `Growth quartile`, y = sd_topt,color = Order))+
  #geom_boxplot(alpha = 0.6, outliers = FALSE, aes(group = max_growth_category))+
  geom_line(aes(group = "Order"))+
  geom_point()+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  facet_grid(~Order)+
  #scale_fill_gradient2(low="#DE5925",mid = "white", high="#0090B5")+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    #strip.text.x.top = ,
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

f114.d <- ave_params_gcplyr_combined_major_taxa_binned_variance %>%
  ggplot(aes(x = `Growth quartile`, y = sd_ctmax,color = Order))+
  #geom_boxplot(alpha = 0.6, outliers = FALSE, aes(group = max_growth_category))+
  geom_line(aes(group = "Order"))+
  geom_point()+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_fill_manual(values = c("#DE5925","#0090B5"))+
  facet_grid(~Order)+
  #scale_fill_gradient2(low="#DE5925",mid = "white", high="#0090B5")+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    #strip.text.x.top = ,
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


S11 <- plot_grid(f114.a,f114.b,f114.c,f114.d, ncol = 2,rel_widths = c(1,1,1,1),align = "hv")



#### S2 - Additional assay to evaluate the influence of carbon concentration on TPC shape ####
#
# 
# We used a subset of the isolate collection, repeated growth measurements
# and TPC fitting in slightly different conditions (e.g., 
# some of our original microplate readers were replaced). Although the values obtained here
# might not be directly comparable to the results obtained with our original methods and equipment
# (hence, we did not add them as replicate measurements of the original dataset),
# we can expect the results of comparing different carbon concentrations within experimental settings
# to be applicable to our original experiment.
# 
#
# Below, we have two clean datasets with TPC fits and parameters

ave_params_CCA_auc_gcplyr <- fread("250915_ave_params_CCA_auc_gcplyr.csv")
ave_preds_CCA_auc_gcplyr <- fread("250915_ave_preds_CCA_auc_gcplyr.csv")



S2A <- ave_preds_CCA_gcplyr_combined %>%
  filter(!seqID %in% c(280)) %>% # This isolate did not revive after -80C storage, but was included in the platemap
  filter(metric != "od") %>%
  ggplot(aes(x = temp, y = .fitted, group = interaction(seqID,Glucose_stock_X),color = as.factor(Glucose_stock_X)))+
  geom_line()+
  scale_color_manual(values = c("purple3","pink3"))+
  facet_grid(metric~.,scales = "free_y")+
  theme_classic(base_size = 12, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


# We add identifiers to each concentration value

ave_params_CCA_gcplyr_combined_wide <- ave_params_CCA_gcplyr_combined %>%
  dplyr::select(!c(e,eh,q10,thermal_safety_margin)) %>%
  pivot_wider(names_from = Glucose_stock_X, values_from = c(rmax,topt,ctmin,ctmax,thermal_tolerance,breadth,skewness),names_sep = "_")



S2B <- ave_params_CCA_gcplyr_combined_wide %>%
  filter(!seqID %in% c(280)) %>%
  filter(metric != "od") %>%
  ggplot(aes(x = rmax_10, y = rmax_100, group = seqID,color = Order))+
  #geom_abline(intercept = c(0,0),slope = 1)+
  geom_smooth(aes(group = model_name),method = "lm")+
  geom_point(aes(size = (rmax_10 + rmax_100) ))+
  scale_color_manual(values = c("orange3","blue3"))+
  facet_grid(~metric,scales = "free")+
  theme_classic(base_size = 14, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

S2C <- ave_params_CCA_gcplyr_combined_wide %>%
  filter(!seqID %in% c(280)) %>%
  filter(metric != "od") %>%
  ggplot(aes(x = topt_10, y = topt_100, group = seqID,color = Order))+
  geom_abline(intercept = c(0,0),slope = 1)+
  geom_smooth(aes(group = model_name),method = "lm")+
  geom_point(aes(size = (rmax_10 + rmax_100) ))+
  #ylim(23,33)+
  #xlim(23,33)+
  scale_color_manual(values = c("orange3","blue3"))+
  facet_grid(metric~.,scales = "free")+
  theme_classic(base_size = 12, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

S2D <- ave_params_CCA_gcplyr_combined_wide %>%
  filter(!seqID %in% c(280)) %>%
  filter(metric != "od") %>%
  ggplot(aes(x = ctmin_10, y = ctmin_100, group = seqID,color = Order))+
  geom_abline(intercept = c(0,0),slope = 1)+
  geom_smooth(aes(group = model_name),method = "lm")+
  geom_point(aes(size = (rmax_10 + rmax_100) ))+
  #ylim(23,33)+
  #xlim(23,33)+
  scale_color_manual(values = c("orange3","blue3"))+
  facet_grid(metric~.,scales = "free")+
  theme_classic(base_size = 12, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

S2E <- ave_params_CCA_gcplyr_combined_wide %>%
  filter(!seqID %in% c(280)) %>%
  filter(metric != "od") %>%
  ggplot(aes(x = ctmax_10, y = ctmax_100, group = seqID,color = Order))+
  geom_abline(intercept = c(0,0),slope = 1)+
  geom_smooth(aes(group = model_name),method = "lm")+
  geom_point(aes(size = (rmax_10 + rmax_100) ))+
  ylim(30,55)+
  xlim(30,55)+
  scale_color_manual(values = c("orange3","blue3"))+
  facet_grid(metric~.,scales = "free")+
  theme_classic(base_size = 12, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    #strip.text = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


S2F <- ave_params_CCA_gcplyr_combined_wide %>%
  filter(!seqID %in% c(280)) %>%
  filter(metric != "od") %>%
  ggplot(aes(x = thermal_tolerance_10, y = thermal_tolerance_100, group = seqID,color = Order))+
  geom_abline(intercept = c(0,0),slope = 1)+
  geom_smooth(aes(group = model_name),method = "lm")+
  geom_point(aes(size = (rmax_10 + rmax_100) ))+
  #ylim(35,55)+
  #xlim(35,55)+
  scale_color_manual(values = c("orange3","blue3"))+
  facet_grid(metric~.,scales = "free")+
  theme_classic(base_size = 12, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    legend.position = "none",
    panel.background = element_rect(color = "black"),
    plot.background = element_rect(fill = 'transparent', color = "white",linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    strip.background = element_blank(),
    #strip.text = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

library(cowplot)

S2A.1 <- S2A

S2A.2 <- plot_grid(S2C,S2D,S2E, 
                   ncol = 3,
                   rel_widths = c(1,1,1.2),
                   rel_heights = c(1,1,1),
                   align = "v",
                   axis = "l")
S2 <- plot_grid(S2A.1,S2A.2, ncol = 2, rel_widths = c(0.5,1))

S2



