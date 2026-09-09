
###### Introduction ######
#
#
# This script is used to generate the main and supplementary figures of the paper
# All needed data is provided in the form of raw files and intermediate, pre-processed files, 
# for parts of the code that are computationally intensive.
# Additionally, output figures and datasets are also provided.
# To run any part of the script, the following three sections should be run first
# Main sections are indicated with six # symbols, as ###### Name of the section 

#### Set directory to parent folder - files from here are called with relative paths

#### Load libraries 
library(extrafont)
library(data.table)
library(tidyverse)
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
library(Hmisc)
library(car)
library(harrietr)

###### Open relevant files and format climate profiles
# Open intermediate files with taxa and thermal traits
ave_params_gcplyr_combined_common_taxa_auc <- fread("05_TCS_phylogenetic_models/03_Output/ave_params_gcplyr_combined_common_taxa_auc.csv.gz")
ave_params_gcplyr_combined_common_taxa_rate <- fread("05_TCS_phylogenetic_models/03_Output/ave_params_gcplyr_combined_common_taxa_rate.csv.gz")

ave_params_gcplyr_combined_major_taxa_auc <- fread("05_TCS_phylogenetic_models/03_Output/ave_params_gcplyr_combined_major_taxa_auc.csv.gz")
ave_params_gcplyr_combined_major_taxa_rate <- fread("05_TCS_phylogenetic_models/03_Output/ave_params_gcplyr_combined_major_taxa_rate.csv.gz")

ave_preds_gcplyr_combined_common_taxa_auc <- fread("05_TCS_phylogenetic_models/03_Output/ave_preds_gcplyr_combined_common_taxa_auc.csv.gz")
ave_preds_gcplyr_combined_common_taxa_rate <- fread("05_TCS_phylogenetic_models/03_Output/ave_preds_gcplyr_combined_common_taxa_rate.csv.gz")

ave_preds_gcplyr_combined_major_taxa_auc <- fread("05_TCS_phylogenetic_models/03_Output/ave_preds_gcplyr_combined_major_taxa_auc.csv.gz")
ave_preds_gcplyr_combined_major_taxa_rate <- fread("05_TCS_phylogenetic_models/03_Output/ave_preds_gcplyr_combined_major_taxa_rate.csv.gz")


# Reading tree file and tree metadata 
tpc_tree <- read.tree(file="05_TCS_phylogenetic_models/03_Output/tpc_tree_kbase_relaxed.tre")
tpc_tree_metadata_auc <- fread("05_TCS_phylogenetic_models/03_Output/tpc_tree_metadata_auc.csv.gz")
tpc_tree_metadata_rate <- fread("05_TCS_phylogenetic_models/03_Output/tpc_tree_metadata_rate.csv.gz")

# We will visualize time series of temperature by site and climate profile in some figures
# We will also visualize the probability of different temperatures across climate profiles
# The following two datasets are obtained as output of the "Environmental_data_processing.R" script

summarized_climate_data <- fread("01_TCS_bioclimatic_variables/03_Output/summarized_climate_data.csv")
continuous_proportion_data <- fread("01_TCS_bioclimatic_variables/03_Output/climate_profile_temperature_probability.csv")

# Here we split the temperature time series by climate profile
cold_stable_climate_time_series <- summarized_climate_data %>% filter(`Climate profile` == "Cold stable")
hot_stable_climate_time_series <- summarized_climate_data %>% filter(`Climate profile` == "Hot stable")
cold_variable_climate_time_series <- summarized_climate_data %>% filter(`Climate profile` == "Cold variable")
hot_variable_climate_time_series <- summarized_climate_data %>% filter(`Climate profile` == "Hot variable")

# We will visualize time series of temperature by site and climate profile during the growth season in some figures
# We will also visualize the probability of different temperatures across climate profiles during the growth season
# The following two datasets are obtained as output of the "Growth_season_calculationg.R" script

summarized_climate_data_growing_season <- fread("01_TCS_bioclimatic_variables/03_Output/site_temperature_probability_growing_season.csv")
continuous_proportion_data_growing_season <- fread("01_TCS_bioclimatic_variables/03_Output/climate_profile_temperature_probability_growing_season.csv")

continuous_proportion_data_growing_season_site <- fread("01_TCS_bioclimatic_variables/03_Output/site_temperature_probability_growing_season.csv")

# Here we split the temperature time series by climate profile
cold_stable_climate_time_series_growing_season <- summarized_climate_data_growing_season %>% filter(`Climate profile` == "Cold stable")
hot_stable_climate_time_series_growing_season <- summarized_climate_data_growing_season %>% filter(`Climate profile` == "Hot stable")
cold_variable_climate_time_series_growing_season <- summarized_climate_data_growing_season %>% filter(`Climate profile` == "Cold variable")
hot_variable_climate_time_series_growing_season <- summarized_climate_data_growing_season %>% filter(`Climate profile` == "Hot variable")

##### Fig S4. Calculate phylogenetic distances
phylogenetic_distance_matrix_auc <- cophenetic.phylo(tpc_tree)
phylogenetic_distance_matrix_rate <- cophenetic.phylo(tpc_tree)

phylogenetic_distance_matrix_auc_df <- melt_dist(phylogenetic_distance_matrix_auc,dist_name = "distance")
phylogenetic_distance_matrix_auc_df$iso1 <- as.numeric(phylogenetic_distance_matrix_auc_df$iso1)
phylogenetic_distance_matrix_auc_df$iso2 <- as.numeric(phylogenetic_distance_matrix_auc_df$iso2)

phylogenetic_distance_matrix_auc_df_2 <- left_join(phylogenetic_distance_matrix_auc_df, tpc_tree_metadata_auc %>%
                                                     dplyr::select(iso1,Class,Order,Family,Genus,Species),
                                                   by = "iso1")%>%
  rename("Class_iso1" = "Class",
         "Order_iso1" = "Order",
         "Family_iso1" = "Family",
         "Genus_iso1" = "Genus",
         "Species_iso1" = "Species")

phylogenetic_distance_matrix_auc_df_2 <- left_join(phylogenetic_distance_matrix_auc_df_2, tpc_tree_metadata_auc %>%
                                                     dplyr::select(iso2,Class,Order,Family,Genus,Species),
                                                   by = "iso2")%>%
  rename("Class_iso2" = "Class",
         "Order_iso2" = "Order",
         "Family_iso2" = "Family",
         "Genus_iso2" = "Genus",
         "Species_iso2" = "Species")

phylogenetic_distance_matrix_auc_df_3 <- phylogenetic_distance_matrix_auc_df_2 %>%
  mutate("Same_Class" = case_when(Class_iso1 == Class_iso2 ~ "Y",
                                  TRUE ~ "N"),
         "Same_Order" = case_when(Order_iso1 == Order_iso2 ~ "Y",
                                  TRUE ~ "N"),
         "Same_Family" = case_when(Family_iso1 == Family_iso2 ~ "Y",
                                  TRUE ~ "N"),
         "Same_Genus" = case_when(Genus_iso1 == Genus_iso2 ~ "Y",
                                  TRUE ~ "N"),
         "Same_Species" = case_when(Species_iso1 == Species_iso2 ~ "Y",
                                  TRUE ~ "N"))


min_distance_dif_Class <- phylogenetic_distance_matrix_auc_df_3 %>%
  filter(Same_Class == "N")%>%
  dplyr::summarize(min_distance = 
                     min(distance))

min_distance_dif_Order <- phylogenetic_distance_matrix_auc_df_3 %>%
  filter(Same_Order == "N")%>%
  dplyr::summarize(min_distance = 
                     min(distance))

min_distance_dif_Family <- phylogenetic_distance_matrix_auc_df_3 %>%
  filter(Same_Family == "N")%>%
  dplyr::summarize(min_distance = 
                     min(distance))

min_distance_dif_Genus <- phylogenetic_distance_matrix_auc_df_3 %>%
  filter(Same_Genus == "N")%>%
  dplyr::summarize(min_distance = 
                     min(distance))

min_distance_dif_Species <- phylogenetic_distance_matrix_auc_df_3 %>%
  filter(Same_Species == "N" & Same_Genus == "Y")%>%
  dplyr::summarize(min_distance = 
                     min(distance))


S4 <- phylogenetic_distance_matrix_auc_df %>%
  ggplot(aes(x = distance))+
  annotate("rect", 
           xmin = 0.0, xmax = min_distance_dif_Species[[1]],   # Part of the x-axis
           ymin = -Inf, ymax = Inf, # All of the y-axis
           fill = "grey60", alpha = 0.2)+ 
  annotate("rect", 
           xmin = min_distance_dif_Genus[[1]], xmax = min_distance_dif_Family[[1]],   # Part of the x-axis
           ymin = -Inf, ymax = Inf, # All of the y-axis
           fill = "grey60", alpha = 0.2)+ 
  annotate("rect", 
           xmin = min_distance_dif_Order[[1]], xmax = Inf,   # Part of the x-axis
           ymin = -Inf, ymax = Inf, # All of the y-axis
           fill = "grey60", alpha = 0.2)+ 
  #annotate(geom = "text", x = min_distance_dif_Species[[1]]/2, y = 1600, label = "Strain")+
  annotate(geom = "text", x = -0.05 + (min_distance_dif_Genus[[1]]-min_distance_dif_Species[[1]])/2, y = 1600, label = "Sp.")+
  annotate(geom = "text", x = min_distance_dif_Genus[[1]]+ (min_distance_dif_Family[[1]]-min_distance_dif_Genus[[1]])/2, y = 1600, label = "Genus")+
  annotate(geom = "text", x = min_distance_dif_Family[[1]]+ (min_distance_dif_Order[[1]]-min_distance_dif_Family[[1]])/2, y = 1600, label = "Family")+
  annotate(geom = "text", x = min_distance_dif_Order[[1]]+ (1.2 -min_distance_dif_Order[[1]])/2, y = 1600, label = "Order")+
  geom_histogram(bins = 400)+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        #axis.title.x = element_blank(),
        #axis.line.x = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )


#### Loading growth curve data used to fit the TPCs ####
# And we also need the growth curve data used as input to fit the TPCs. 
# This file was generated in the first section of the "TCS_TPCs.R" script,
# filtering out low quality TPCs for PCA plotting

tpc_input_gcplyr <- fread("04_TCS_thermal_performance_curves/02_Output/tpc_input_gcplyr_average.csv.gz")

ave_params_gcplyr_combined <- fread("04_TCS_thermal_performance_curves/02_Output/251029_tcs_params_gcplyr.csv.gz")
ave_preds_gcplyr_combined <- fread("04_TCS_thermal_performance_curves/02_Output/251029_tcs_fits_gcplyr.csv.gz")

ave_preds_gcplyr_combined <- ave_preds_gcplyr_combined %>%
  mutate(.fitted = if_else(.fitted < 0, 0, .fitted))
  

ave_params_gcplyr_combined_major_taxa <- ave_params_gcplyr_combined %>% 
  filter(Order %in% c("Enterobacterales","Pseudomonadales")) %>%
  droplevels()

ave_preds_gcplyr_combined_major_taxa <- ave_preds_gcplyr_combined %>% 
  filter(Order %in% c("Enterobacterales","Pseudomonadales")) %>%
  droplevels()

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

#### Summarizing the growth curve data by Order, temperature, and metric ####

tpc_input_gcplyr_variance_auc <- tpc_input_gcplyr_auc %>%
  group_by(temp, Order)%>%
  dplyr::summarize(n = n(),
                   sd_auc = sd(auc),
                   sd_rate = sd(rate),
                   auc = mean(auc),
                   rate = mean(rate)
  )


wilcoxon_tests_variance_auc <- tpc_input_gcplyr_auc %>%
  filter(Order %in% c("Enterobacterales","Pseudomonadales")) %>%
  group_by(Order, temp) %>%
  dplyr::summarize(
    median_auc = median(auc, na.rm = TRUE),
    n = sum(!is.na(auc)),
    auc_values = list(auc),
    .groups = "drop"
  ) %>%
  group_by(Order) %>%
  mutate(
    mu = 0.10 * max(median_auc, na.rm = TRUE) # Obtaining the median value
  ) %>%
  
  rowwise() %>% # running the Wilcoxon tests per temperature and Order
  mutate(
    
    lower_than_mu = sum(unlist(auc_values) < mu), # values lower than the threshold
    
    test = list(
      wilcox.test(
        unlist(auc_values),
        mu = mu,
        alternative = "less",
        conf.int = TRUE
      )
    ),
    
    p_value = test$p.value,
    conf_low = test$conf.int[1],
    conf_high = test$conf.int[2]
  ) %>%
  
  ungroup() %>%
  select(-auc_values, -test)



tpc_input_gcplyr_variance_auc <- tpc_input_gcplyr_auc %>%
  group_by(temp, Order)%>%
  dplyr::summarize(n = n(),
                   sd_auc = sd(auc),
            sd_rate = sd(rate),
            auc = mean(auc),
            rate = mean(rate)
  )

tpc_input_gcplyr_variance_rate <- tpc_input_gcplyr_rate %>%
  group_by(temp, Order)%>%
  dplyr::summarize(n = n(),
                   sd_auc = sd(auc),
            sd_rate = sd(rate),
            auc = mean(auc),
            rate = mean(rate)
  )


#### S10.top - Plotting the distribution and variance of AUC growth curves ####

S10C1.auc <- tpc_input_gcplyr_variance_auc %>% 
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


S10C2.auc <- tpc_input_gcplyr_variance_auc %>% 
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

S10C3.auc <- tpc_input_gcplyr_auc %>% 
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
  annotation_custom(ggplotGrob(S10C1.auc),xmin = 30, xmax = 44, 
                    ymin = 2.8, ymax = 4.5)


S10C4.auc <- tpc_input_gcplyr_auc %>% 
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
  annotation_custom(ggplotGrob(S10C2.auc),xmin = 30, xmax = 44, 
                    ymin = 2.8, ymax = 4.5)

# We then arrange the plots into a main figure
S10.top <- plot_grid(S10C3.auc + 
                       theme(legend.position = "none"),S10C4.auc + 
                       theme(legend.position = "none"), ncol = 2, align = "h", rel_widths = c(1, 1),axis = "l")



#### S10.bot - Plotting the distribution and variance of Rate growth curves ####

S10C1.rate <- tpc_input_gcplyr_variance_rate %>% 
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


S10C2.rate <- tpc_input_gcplyr_variance_rate %>% 
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

S10C3.rate <- tpc_input_gcplyr_rate %>% 
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
  annotation_custom(ggplotGrob(S10C1.rate),xmin = 30, xmax = 44, 
                    ymin = 1.0, ymax = 1.5)


S10C4.rate <- tpc_input_gcplyr_rate %>% 
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
  annotation_custom(ggplotGrob(S10C2.rate),xmin = 30, xmax = 44, 
                    ymin = 1, ymax = 1.5)

# We then arrange the plots into a main figure
S10.bot <- plot_grid(S10C3.rate + 
                       theme(legend.position = "none"),S10C4.rate + 
                       theme(legend.position = "none"), ncol = 2, align = "h", rel_widths = c(1, 1),axis = "l")


###### F2A.left - California map ######
#
# Here we generate a preliminary figure that shows the location of the sampled sites in a map of California
# We also map basic bioclimatic variables onto it
#
#
# We load a shapefile of California
california_map <- st_read("06_TCS_figures/02_Data/CA_State.shp")
california_geometry <- st_geometry(california_map)
# We obtain the coordinates of each sampled site from a database
environment_database_noprec <- fread("01_TCS_bioclimatic_variables/03_Output/environment_database_tcs_noprec_250512.csv")
environment_database_full <- fread("01_TCS_bioclimatic_variables/03_Output/environment_database_tcs_full_250512.csv.gz")

reserves_coordinates <- fread("06_TCS_figures/02_Data/SitesMetadata.csv")
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

F2A.left <- ggplot() +
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

###### F2A.right - Dendrograms ######
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


F2A.right <- ggtree(tpc_tree, layout = "circular", branch.length = "none") %<+% tpc_tree_metadata_auc +
  #geom_tiplab(size = 1.5, # color for label font
  #            geom = "label",  # labels not text
  #            offset = 0.3,
  #            label.padding = unit(0, "lines"), # amount of padding around the labels
  #            label.size = 0) + # size of label border
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  geom_tippoint(size = 4, aes(color = Site_initials))+
  #geom_tiplab(offset = 0.5, size = 1)+
  #scale_shape_manual(values = c(16,15,18,17))+
  #scale_color_viridis_d(option = "inferno", guide = "none")+
  scale_color_manual(values = full_sites_palette)+
  theme(legend.position = "none", 
        legend.background = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank())


# Add unique Order labels to the outer ring
F2A.right <- F2A.right +
  geom_strip(4, 6, barsize=1, offset = 2, label = "Pseudomonadales",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(2, 4, barsize=1, offset = 1.5, label = "Burkholderiales",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(6, 49, barsize=1, offset = 1.5, label = "Enterobacterales",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(1, 2, barsize=1, offset = 1, label = "Xanthomonadales",offset.text = 1, align = TRUE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT")+
  geom_strip(50, 252, barsize=1, offset = 1.5, label = "Pseudomonadales",offset.text = 15, align = TRUE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT")


F2A.right


S1 <- ggtree(tpc_tree, layout = "circular", branch.length = "none") %<+% 
  tpc_tree_metadata_auc +
  #geom_tiplab(size = 1.5, # color for label font
  #            geom = "label",  # labels not text
  #            offset = 0.3,
  #            label.padding = unit(0, "lines"), # amount of padding around the labels
  #            label.size = 0) + # size of label border
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  geom_tippoint(size = 2, aes(color = Site_initials))+
  #geom_tiplab(offset = 0.5, size = 1)+
  #scale_shape_manual(values = c(16,15,18,17))+
  #scale_color_viridis_d(option = "inferno", guide = "none")+
  scale_color_manual(values = full_sites_palette)+
  theme(legend.position = "none", 
        legend.background = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank())

S1 <- ggtree(tpc_tree, layout = "circular", branch.length = "none") %<+% tpc_tree_metadata_auc +
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
  theme(legend.position = "bottom", 
        legend.background = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank())


# Add unique Order labels to the outer ring
S1 <- S1 +
  geom_strip(1, 2, barsize=1, offset = 1, label = "Stenotrophomonas",offset.text = 1, align = TRUE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT")+
  geom_strip(2, 3, barsize=1, offset = 2, label = "Telluria",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT",color = "grey50") +
  geom_strip(3, 4, barsize=1, offset = 1, label = "Paraburkholderia",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(4, 6, barsize=1, offset = 2, label = "Acinetobacter",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT",color = "grey50") +
  geom_strip(6, 7, barsize=1, offset = 1, label = "Aeromonas",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(7, 8, barsize=1, offset = 2, label = "Shewanella",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT",color = "grey50") +
  geom_strip(8, 15, barsize=1, offset = 1, label = "Pantoea",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(15,20, barsize=1, offset = 2, label = "Mixta",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT",color = "grey50") +
  geom_strip(20, 22, barsize=1, offset = 1, label = "Erwinia",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(22, 24, barsize=1, offset = 2, label = "Rahnella",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT",color = "grey50") +
  geom_strip(24, 30, barsize=1, offset = 1, label = "Serratia_J",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(30, 31, barsize=1, offset = 2, label = "Buttiauxella",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT",color = "grey50") +
  geom_strip(31, 32, barsize=1, offset = 1, label = "Cronobacter",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(32, 34, barsize=1, offset = 2, label = "Siccibacter",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT",color = "grey50") +
  geom_strip(34, 35, barsize=1, offset = 1, label = "Scandinavium",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(35, 38, barsize=1, offset = 2, label = "UBA7405",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT",color = "grey50") +
  geom_strip(38, 40, barsize=1, offset = 1, label = "Lelliottia",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT") +
  geom_strip(40, 50, barsize=1, offset = 2, label = "Enterobacter",offset.text = 1, align = FALSE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT",color = "grey50") +
  geom_strip(50, 190, barsize=1, offset = 1, label = "Pseudomonas_E",offset.text = 15, align = TRUE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT")+
  geom_strip(190, 191, barsize=1, offset = 2, label = "Pseudomonas_B",offset.text = 15, align = TRUE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT",color = "grey50")+
  geom_strip(191, 205, barsize=1, offset = 1, label = "Stutzerimonas",offset.text = 15, align = TRUE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT")+
  geom_strip(205, 252, barsize=1, offset = 2, label = "Pseudomonas_E",offset.text = 15, align = TRUE, fontsize = 4, angle = 0, geom = "text", hjust = 0, fill = NA, family = "Futura Bk BT",color = "grey50")


###### F2B - Visual comparison of thermal traits across Orders ######
#
#
# In this code section we combine the TPC fits and parameters at the Order level
# colored by climate profile for a subset of relatively common clades (Genera >1%)
#
#
#### F2C AUC and S6 ####

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

F2B.auc <- ave_preds_gcplyr_combined_major_taxa_auc %>% 
  ggplot( aes(x = temp, y = .fitted)) +
  geom_line(lwd = 0.75,alpha = 0.8,aes( group = seqID,color = Site_initials)) +
  geom_ribbon(data = ave_preds_gcplyr_combined_major_taxa_auc %>% 
                group_by(Order,temp) %>%
                dplyr::summarize(sd_fitted = sd(.fitted),
                          .fitted = mean(.fitted)),
              aes(x = temp, ymin = .fitted - sd_fitted, 
                  ymax = .fitted + sd_fitted ),lwd = 2,fill = "grey40",alpha = 0.4)+
  geom_line(data = ave_preds_gcplyr_combined_major_taxa_auc %>% 
              group_by(Order,temp) %>%
              dplyr::summarize(.fitted = mean(.fitted)),
            aes(x = temp, y = .fitted),lwd = 1.3,color = "grey95")+
  facet_wrap(~factor(Order),ncol=1)+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)+
  xlim(13,45)+
  ylim(0,7)+
  ylab(expression(AUC))+
  #scale_color_viridis_d()+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none"
  )

F2C1.auc <- ave_params_gcplyr_combined_major_taxa_auc %>% 
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

F2C2.auc <- ave_params_gcplyr_combined_major_taxa_auc %>%
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


F2C3.auc <- ave_params_gcplyr_combined_major_taxa_auc %>%
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

F2C4.auc <- ave_params_gcplyr_combined_major_taxa_auc %>%
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
F2C.auc <- plot_grid(F2C1.auc,F2C2.auc,F2C3.auc,F2C4.auc, ncol = 1, align = "v", rel_heights = c(1.3,1.3,1.3,1.4),axis = "l")
# Show combined plot
F2C.auc


# Combine plots
F2BC.auc <- plot_grid(F2B.auc,F2C.auc, ncol = 2, align = "h", rel_widths = c(1.3,1.3),axis = "l")
# Show combined plot
F2BC.auc


##### S22B.left #####

# This is a figure displaying the TPCs of the two main Orders,
# colored by site

F2B.rate <- ave_preds_gcplyr_combined_major_taxa_rate %>% 
  ggplot( aes(x = temp, y = .fitted)) +
  geom_line(lwd = 0.75,alpha = 0.8,aes( group = seqID,color = Site_initials)) +
  geom_ribbon(data = ave_preds_gcplyr_combined_major_taxa_rate %>% 
                group_by(Order,temp) %>%
                dplyr::summarize(sd_fitted = sd(.fitted),
                                 .fitted = mean(.fitted)),
              aes(x = temp, ymin = .fitted - sd_fitted, 
                  ymax = .fitted + sd_fitted ),lwd = 2,fill = "grey40",alpha = 0.4)+
  geom_line(data = ave_preds_gcplyr_combined_major_taxa_rate %>% 
              group_by(Order,temp) %>%
              dplyr::summarize(.fitted = mean(.fitted)),
            aes(x = temp, y = .fitted),lwd = 1.3,color = "grey95")+
  facet_wrap(~factor(Order),ncol=1)+
  #scale_color_viridis_d(option = "inferno")+
  scale_color_manual(values = full_sites_palette)+
  xlim(13,45)+
  ylim(0,1)+
  ylab(expression(rate))+
  #scale_color_viridis_d()+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none"
  )

F2C1.rate <- ave_params_gcplyr_combined_major_taxa_rate %>% 
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

F2C2.rate <- ave_params_gcplyr_combined_major_taxa_rate %>%
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


F2C3.rate <- ave_params_gcplyr_combined_major_taxa_rate %>%
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

F2C4.rate <- ave_params_gcplyr_combined_major_taxa_rate %>%
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
F2C.rate <- plot_grid(F2C1.rate,F2C2.rate,F2C3.rate,F2C4.rate, ncol = 1, align = "v", rel_heights = c(1.3,1.3,1.3,1.4),axis = "l")
# Show combined plot
F2C.rate


# Combine plots
S22BC <- plot_grid(F2B.rate,F2C.rate, ncol = 2, align = "h", rel_widths = c(1.3,1.3),axis = "l")
# Show combined plot
S22BC

#### Basic statistics for F2C ####
#
#
# Here we run basic one-way ANOVAs to compare thermal traits between the plotted Orders,
# Climate profiles, and post-hoc pairwise comparisons using Tukey 
#
#

# Stats for F2C.auc by Order

anova_model_Order_aucmax <- aov(rmax ~ Order, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_Order_aucmax)
posthoc_Order_aucmax <- TukeyHSD(anova_model_Order_aucmax)
posthoc_Order_aucmax

emm_aucmax <- emmeans(anova_model_Order_aucmax, ~ Order)

d_aucmax <- eff_size(
  emm,
  sigma = sigma(anova_model_Order_aucmax),
  edf = df.residual(anova_model_Order_aucmax)
)


anova_model_Order_topt_auc <- aov(topt ~ Order, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_Order_topt_auc)
posthoc_Order_topt_auc <- TukeyHSD(anova_model_Order_topt_auc)
posthoc_Order_topt_auc

d_auc_topt <- eff_size(
  emm,
  sigma = sigma(anova_model_Order_topt_auc),
  edf = df.residual(anova_model_Order_topt_auc)
)

anova_model_Order_ctmax_auc <- aov(ctmax ~ Order, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_Order_ctmax_auc)
posthoc_Order_ctmax_auc <- TukeyHSD(anova_model_Order_ctmax_auc)
posthoc_Order_ctmax_auc

d_auc_ctmax <- eff_size(
  emm,
  sigma = sigma(anova_model_Order_ctmax_auc),
  edf = df.residual(anova_model_Order_ctmax_auc)
)

anova_model_Order_ctmin_auc <- aov(ctmin ~ Order, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_Order_ctmin_auc)
posthoc_Order_ctmin_auc <- TukeyHSD(anova_model_Order_ctmin_auc)
posthoc_Order_ctmin_auc

d_auc_ctmin <- eff_size(
  emm,
  sigma = sigma(anova_model_Order_ctmin_auc),
  edf = df.residual(anova_model_Order_ctmin_auc)
)

# Stats for F2C.auc by Climate Profile

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


# Stats for S22C by Order

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

# Stats for S22C by Climate Profile

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

###### F2C.right.auc - Calculating Phylogenetic Signal for major taxa (>1% of collection) AUC ######

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
  formula <- as.formula(
    paste(
      trait_name,
      "~ PC1 + PC2 + temperature_of_isolation + (1|Order/Genus) + (1|Site)"
    )
  )
  
  model <- lmer(formula, data = data)
  
  # Global model to extract Marginal and Conditional R2 values
  
  r2_values <- MuMIn::r.squaredGLMM(model)
  
  model_fit <- data.frame(
    Trait = trait_name,
    metric = "auc",
    Variable = c("Marginal_R2", "Conditional_R2"),
    Effect = "model_fit",
    Variance_observed = c(
      r2_values[1, "R2m"],
      r2_values[1, "R2c"]
    )
  )
  
  ## Partitioning the Marginal P-values and R2
  
  part_model <- partR2(
    model,
    partvars = c("PC1", "PC2", "temperature_of_isolation"),
    data = data
  )
  
  ## Extract fixed effects table
  fixed_summary <- as.data.frame(summary(model)$coefficients)
  fixed_summary$grp <- rownames(fixed_summary)
  rownames(fixed_summary) <- NULL
  
  ## Keep only fixed effects included in partR2 output
  fixed_var <- fixed_summary[
    fixed_summary$grp %in% part_model$R2$term,
  ]
  
  ## Match partial R² values
  fixed_var$share <- part_model$R2$estimate[
    match(fixed_var$grp, part_model$R2$term)
  ]
  
  ## Partitioning the random effects
  observed_var <- as.data.frame(VarCorr(model))
  
  observed_total <- sum(observed_var$vcov)
  
  observed_var$proportion <- observed_var$vcov / observed_total
  
  # Running a permutation test to calculate the robustness of our phylogenetic signal
  permuted_proportions <- matrix(
    NA,
    nrow = n_permutations,
    ncol = nrow(observed_var)
  )
  
  colnames(permuted_proportions) <- observed_var$grp
  
  for (i in seq_len(n_permutations)) {
    
    perm_data <- data
    perm_data[[trait_name]] <- sample(perm_data[[trait_name]])
    
    perm_model <- lmer(formula, data = perm_data)
    
    perm_var <- as.data.frame(VarCorr(perm_model))
    
    permuted_proportions[i, ] <-
      perm_var$vcov / sum(perm_var$vcov)
  }
  
  ci_lower <- apply(
    permuted_proportions,
    2,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  )
  
  ci_upper <- apply(
    permuted_proportions,
    2,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  )
  
  confint_df <- data.frame(
    observed = observed_var$proportion,
    lower_ci = ci_lower,
    upper_ci = ci_upper,
    Factor = observed_var$grp,
    Trait = trait_name,
    significant =
      observed_var$proportion < ci_lower |
      observed_var$proportion > ci_upper
  )
  
  confint_df$Factor <- dplyr::recode(
    confint_df$Factor,
    "Order:Genus" = "Genus",
    .default = confint_df$Factor
  )
  
  list(
    fixed_var = fixed_var,
    phylo_signal = confint_df,
    model_fit = model_fit
  )
}

phylogenetic_signal_pipeline <- function(df) {
  
  ## Run model for each trait
  results_aucmax <- phylogenetic_signal_mixed_model_f("rmax", df)
  results_topt   <- phylogenetic_signal_mixed_model_f("topt", df)
  results_ctmax  <- phylogenetic_signal_mixed_model_f("ctmax", df)
  results_ctmin  <- phylogenetic_signal_mixed_model_f("ctmin", df)
  
  ## Combine all data into a single dataframe
  phylogenetic_signal_thermal_traits_auc <- bind_rows(
    results_aucmax$phylo_signal,
    results_topt$phylo_signal,
    results_ctmax$phylo_signal,
    results_ctmin$phylo_signal
  )
  
  ## Combine fixed effects into a single dataframe
  fixed_effects_all_auc <- bind_rows(
    cbind(Trait = "rmax",  results_aucmax$fixed_var),
    cbind(Trait = "topt",  results_topt$fixed_var),
    cbind(Trait = "ctmax", results_ctmax$fixed_var),
    cbind(Trait = "ctmin", results_ctmin$fixed_var)
  )
  
  ## Combine global R2 values
  model_fit_all <- bind_rows(
    results_aucmax$model_fit,
    results_topt$model_fit,
    results_ctmax$model_fit,
    results_ctmin$model_fit
  )
  
  ## Combine random effects into a single dataframe
  random_effects_clean_auc <- phylogenetic_signal_thermal_traits_auc %>%
    mutate(
      metric = "auc",
      Variable = Factor,
      Effect = "random"
    ) %>%
    rename(
      Variance_observed = observed,
      Lower_CI = lower_ci,
      Upper_CI = upper_ci,
      Significant = significant
    ) %>%
    dplyr::select(
      Trait,
      metric,
      Variable,
      Effect,
      Variance_observed,
      Lower_CI,
      Upper_CI,
      Significant
    )
  
  ## Combine fixed effects into a single df too
  fixed_effects_clean_auc <- fixed_effects_all_auc %>%
    mutate(
      metric = "auc",
      Variable = grp,
      Effect = "fixed",
      p_value = `Pr(>|t|)`
    ) %>%
    rename(
      Estimate = Estimate,
      Std_Error = `Std. Error`,
      t_value = `t value`,
      Partial_R2 = share
    ) %>%
    dplyr::select(
      Trait,
      metric,
      Variable,
      Effect,
      Estimate,
      Std_Error,
      t_value,
      p_value,
      Partial_R2
    )
  
  ## Select relevant columns from the global fit df
  model_fit_clean_auc <- model_fit_all %>%
    mutate(
      Lower_CI = NA_real_,
      Upper_CI = NA_real_,
      Significant = NA
    ) %>%
    dplyr::select(
      Trait,
      metric,
      Variable,
      Effect,
      Variance_observed,
      Lower_CI,
      Upper_CI,
      Significant
    )
  
  ## Put everything together into Supplementary Table 2
  combined_effects_lmer_auc <- bind_rows(
    fixed_effects_clean_auc,
    random_effects_clean_auc,
    model_fit_clean_auc
  ) %>%
    mutate(
      across(where(is.numeric), ~ round(.x, 3))
    )
  
  return(
    list(
      combined_effects = combined_effects_lmer_auc,
      phylo_signal = phylogenetic_signal_thermal_traits_auc,
      fixed_effects_all_auc = fixed_effects_all_auc,
      model_fit = model_fit_all
    )
  )
}

# Running pipeline 
combined_effects_lmer_auc <-
  phylogenetic_signal_pipeline(
    ave_params_gcplyr_combined_major_taxa_auc
  )


write.csv(combined_effects_lmer_auc$combined_effects, "06_TCS_figures/03_Output/combined_effects_lmer_auc.csv", row.names = FALSE)
write.csv(combined_effects_lmer_auc$phylo_signal,"06_TCS_figures/03_Output/phylogenetic_signal_thermal_traits_auc.csv",row.names = FALSE)
write.csv(combined_effects_lmer_auc$fixed_effects_all_auc,"06_TCS_figures/03_Output/fixed_effects_lmer_auc.csv",row.names = FALSE)


#combined_effects_lmer_auc <- fread("06_TCS_figures/03_Output/combined_effects_lmer_auc_full.CSV")
#phylogenetic_signal_thermal_traits_auc <- fread("06_TCS_figures/03_Output/phylogenetic_signal_thermal_traits_auc.csv")


phylogenetic_signal_thermal_traits_auc <- phylogenetic_signal_thermal_traits_auc %>%
  mutate(Factor = dplyr::recode(Factor,
                           "Genus:Order" = "Genus"))

phylogenetic_signal_thermal_traits_auc$Factor <- factor(phylogenetic_signal_thermal_traits_auc$Factor, 
                                                                 levels = c("Genus","Order","Residual","Site"))


F2C.right.auc <- phylogenetic_signal_thermal_traits_auc %>% 
  #filter(!`Taxonomic Level` == "Family") %>% 
  mutate(Trait = fct_relevel(Trait,"ctmin","ctmax","topt","rmax")) %>%
  #mutate(Trait = fct_relevel(Trait,"ctmin_auc","ctmax_auc","topt_auc","aucmax" )) %>%
  ggplot(aes(x = Trait, y = observed)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), linewidth = 0.5,width = 0.2, color = "gray50") +
  geom_point(aes(color = significant), size = 6) +
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  ylim(0,1)+
  scale_y_continuous(breaks = seq(0, 1, 0.5))+
  coord_flip()+
  facet_wrap( ~Factor,ncol = 4)+
  labs(
    y = "Proportion of Variance Explained",
    color = "Significant"
  ) +
  scale_color_manual(values = c("grey60","black"))+
  theme(strip.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "top",
        legend.direction = "horizontal"
  )+
  guides(color = "none")

# Assembling the bottom panel of figure F2

F2BCD.auc <- plot_grid(F2BC.auc,F2C.right.auc, ncol = 2, align = "h", rel_widths = c(1.3,1.3),axis = "l")

F2BCD.auc

# Assembling the bottom panel of figure S22

S22BC <- plot_grid(S22BC,S22B.right, ncol = 2, align = "h", rel_widths = c(1.3,1.3),axis = "l")

S22BC


###### LMERS for non-composite environmental variables ######

#
#
# Here we apply the methods detailed in Lennon et al (2016) to calculate the phylogenetic
# signal of bacterial traits to our thermal traits, by fitting hierarchical linear models
# for each trait and metric.
#
#


n_permutations <- 10000

phylogenetic_signal_mixed_model_simple_vars_f <- function(trait_name, data) {
  
  # Fit hierarchical linear model with Satterthwaite df
  formula <- as.formula(paste(trait_name, "~ MAT + MATR + temperature_of_isolation + (1|Order/Genus) + (1|Site)"))
  model <- lmer(formula, data = data)
  
  # Partial marginal R² partitioning
  part_model <- partR2(model, partvars = c("MAT","MATR","temperature_of_isolation"), data = data)
  
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
    perm_data <- data
    perm_data[[trait_name]] <- sample(perm_data[[trait_name]])
    
    perm_model <- lmer(formula, data = perm_data)         
    perm_var <- as.data.frame(VarCorr(perm_model))   
    permuted_proportions[i, ] <- perm_var$vcov / sum(perm_var$vcov)
  }
  
  ci_lower <- apply(permuted_proportions, 2, quantile, probs = 0.025)
  ci_upper <- apply(permuted_proportions, 2, quantile, probs = 0.975)
  
  confint_df <- data.frame(
    observed = observed_var$proportion,
    lower_ci = ci_lower,
    upper_ci = ci_upper,
    Factor = observed_var$grp,
    Trait = trait_name,
    significant = observed_var$proportion < ci_lower | observed_var$proportion > ci_upper
  )
  
  confint_df$Factor <- recode(
    confint_df$Factor,
    "Order:Genus" = "Genus",
    .default = confint_df$Factor
  )
  
  list(fixed_var = fixed_var, phylo_signal = confint_df)
}


phylogenetic_signal_pipeline_simple <- function(df) {
  
  # Apply to all traits
  results_aucmax <- phylogenetic_signal_mixed_model_simple_vars_f("rmax", df)
  results_topt   <- phylogenetic_signal_mixed_model_simple_vars_f("topt", df)
  results_ctmax  <- phylogenetic_signal_mixed_model_simple_vars_f("ctmax", df)
  results_ctmin  <- phylogenetic_signal_mixed_model_simple_vars_f("ctmin", df)
  
  # Combine phylogenetic signal data
  phylogenetic_signal_thermal_traits_auc_simple <- bind_rows(
    results_aucmax$phylo_signal,
    results_topt$phylo_signal,
    results_ctmax$phylo_signal,
    results_ctmin$phylo_signal
  )
  
  # Combine fixed effect contributions
  fixed_effects_all_auc_simple <- bind_rows(
    cbind(Trait="rmax", results_aucmax$fixed_var),
    cbind(Trait="topt", results_topt$fixed_var),
    cbind(Trait="ctmax", results_ctmax$fixed_var),
    cbind(Trait="ctmin", results_ctmin$fixed_var)
  )
  
  # Clean random effects
  random_effects_clean_auc_simple <- phylogenetic_signal_thermal_traits_auc_simple %>%
    mutate(
      metric = "auc",
      Variable = `Factor`,
      Effect = "random"
    ) %>%
    rename(
      Variance_observed = observed,
      Lower_CI = lower_ci,
      Upper_CI = upper_ci,
      Significant = significant
    ) %>%
    dplyr::select(Trait, metric, Variable, Effect, Variance_observed, Lower_CI, Upper_CI, Significant)
  
  # Clean fixed effects
  fixed_effects_clean_auc_simple <- fixed_effects_all_auc_simple %>%
    mutate(
      metric = "auc",
      Variable = grp,
      Effect = "fixed",
      p_value = `Pr(>|t|)`
    ) %>%
    rename(
      Estimate = Estimate,
      Std_Error = `Std. Error`,
      t_value = `t value`,
      Partial_R2 = share
    ) %>%
    dplyr::select(Trait, metric, Variable, Effect, Estimate, Std_Error, t_value, p_value, Partial_R2)
  
  # Combine
  combined_effects_lmer_auc_simple <- bind_rows(
    fixed_effects_clean_auc_simple,
    random_effects_clean_auc_simple
  ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  return(list(
    combined_effects = combined_effects_lmer_auc_simple,
    phylo_signal = phylogenetic_signal_thermal_traits_auc_simple,
    fixed_effects_all_auc = fixed_effects_all_auc_simple
  ))
}


combined_effects_lmer_auc_simple <- phylogenetic_signal_pipeline_simple(ave_params_gcplyr_combined_major_taxa_auc)


write.csv(combined_effects_lmer_auc_simple$combined_effects, "06_TCS_figures/03_Output/combined_effects_lmer_auc_simple.csv", row.names = FALSE)
write.csv(combined_effects_lmer_auc_simple$phylo_signal,"06_TCS_figures/03_Output/phylogenetic_signal_thermal_traits_auc_simple.csv",row.names = FALSE)
write.csv(combined_effects_lmer_auc_simple$fixed_effects_all_auc,"06_TCS_figures/03_Output/fixed_effects_lmer_auc_simple.csv",row.names = FALSE)

#### Merging the table of models ####


combined_effects_lmer_auc_simple_combined_effects <- fread("06_TCS_figures/03_Output/combined_effects_lmer_auc_simple.csv")

combined_effects_lmer_auc$model <- "lmerTest::lmer(Trait ~ PC1 + PC2 + temperature_of_isolation + (1|Order/Genus) + (1|Site))"
combined_effects_lmer_auc_simple_combined_effects$model <- "lmerTest::lmer(Trait ~ MAT + MATR + temperature_of_isolation + (1|Order/Genus) + (1|Site))"

combined_effects_lmer_auc_full <- bind_rows(combined_effects_lmer_auc,
                                            combined_effects_lmer_auc_simple_combined_effects %>% filter(Effect == "fixed"))


#### F3B ####

fixed_effects_all_auc <- fread("06_TCS_figures/03_Output/fixed_effects_lmer_auc.csv")

fixed_effects_all_auc <- fixed_effects_all_auc %>%
  mutate(Factor = dplyr::recode(grp,
                         "temperature_of_isolation" = "TOI"))

F3B.auc <- fixed_effects_all_auc %>%
  mutate(Factor = fct_relevel(Factor,"TOI","PC2","PC1")) %>%
  mutate(Trait = fct_relevel(Trait,"rmax","topt","ctmax","ctmin")) %>%
  ggplot(aes(x = `t value`, y = Factor,size = share*100))+
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

###### S22B.right - Calculating Phylogenetic Signal for major taxa (>1% of collection) Rate ######

#
#
# Here we apply the methods detailed in Lennon et al (2016) to calculate the phylogenetic
# signal of bacterial traits to our thermal traits, by fitting hierarchical linear models
# for each trait and metric.
#
#


n_permutations <- 10000

phylogenetic_signal_mixed_model_f_rate <- function(trait_name, data) {
  
  # Fit hierarchical linear model with Satterthwaite df
  formula <- as.formula(
    paste(
      trait_name,
      "~ PC1 + PC2 + temperature_of_isolation + (1|Order/Genus) + (1|Site)"
    )
  )
  
  model <- lmer(formula, data = data)
  
  # Global model to extract Marginal and Conditional R2 values
  
  r2_values <- MuMIn::r.squaredGLMM(model)
  
  model_fit <- data.frame(
    Trait = trait_name,
    metric = "rate",
    Variable = c("Marginal_R2", "Conditional_R2"),
    Effect = "model_fit",
    Variance_observed = c(
      r2_values[1, "R2m"],
      r2_values[1, "R2c"]
    )
  )
  
  ## Partitioning the Marginal P-values and R2
  
  part_model <- partR2(
    model,
    partvars = c("PC1", "PC2", "temperature_of_isolation"),
    data = data
  )
  
  ## Extract fixed effects table
  fixed_summary <- as.data.frame(summary(model)$coefficients)
  fixed_summary$grp <- rownames(fixed_summary)
  rownames(fixed_summary) <- NULL
  
  ## Keep only fixed effects included in partR2 output
  fixed_var <- fixed_summary[
    fixed_summary$grp %in% part_model$R2$term,
  ]
  
  ## Match partial R² values
  fixed_var$share <- part_model$R2$estimate[
    match(fixed_var$grp, part_model$R2$term)
  ]
  
  ## Partitioning the random effects
  observed_var <- as.data.frame(VarCorr(model))
  
  observed_total <- sum(observed_var$vcov)
  
  observed_var$proportion <- observed_var$vcov / observed_total
  
  # Running a permutation test to calculate the robustness of our phylogenetic signal
  permuted_proportions <- matrix(
    NA,
    nrow = n_permutations,
    ncol = nrow(observed_var)
  )
  
  colnames(permuted_proportions) <- observed_var$grp
  
  for (i in seq_len(n_permutations)) {
    
    perm_data <- data
    perm_data[[trait_name]] <- sample(perm_data[[trait_name]])
    
    perm_model <- lmer(formula, data = perm_data)
    
    perm_var <- as.data.frame(VarCorr(perm_model))
    
    permuted_proportions[i, ] <-
      perm_var$vcov / sum(perm_var$vcov)
  }
  
  ci_lower <- apply(
    permuted_proportions,
    2,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  )
  
  ci_upper <- apply(
    permuted_proportions,
    2,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  )
  
  confint_df <- data.frame(
    observed = observed_var$proportion,
    lower_ci = ci_lower,
    upper_ci = ci_upper,
    Factor = observed_var$grp,
    Trait = trait_name,
    significant =
      observed_var$proportion < ci_lower |
      observed_var$proportion > ci_upper
  )
  
  confint_df$Factor <- dplyr::recode(
    confint_df$Factor,
    "Order:Genus" = "Genus",
    .default = confint_df$Factor
  )
  
  list(
    fixed_var = fixed_var,
    phylo_signal = confint_df,
    model_fit = model_fit
  )
}

phylogenetic_signal_pipeline_rate <- function(df) {
  
  ## Run model for each trait
  results_ratemax <- phylogenetic_signal_mixed_model_f_rate("rmax", df)
  results_topt   <- phylogenetic_signal_mixed_model_f_rate("topt", df)
  results_ctmax  <- phylogenetic_signal_mixed_model_f_rate("ctmax", df)
  results_ctmin  <- phylogenetic_signal_mixed_model_f_rate("ctmin", df)
  
  ## Combine all data into a single dataframe
  phylogenetic_signal_thermal_traits_rate <- bind_rows(
    results_ratemax$phylo_signal,
    results_topt$phylo_signal,
    results_ctmax$phylo_signal,
    results_ctmin$phylo_signal
  )
  
  ## Combine fixed effects into a single dataframe
  fixed_effects_all_rate <- bind_rows(
    cbind(Trait = "rmax",  results_ratemax$fixed_var),
    cbind(Trait = "topt",  results_topt$fixed_var),
    cbind(Trait = "ctmax", results_ctmax$fixed_var),
    cbind(Trait = "ctmin", results_ctmin$fixed_var)
  )
  
  ## Combine global R2 values
  model_fit_all <- bind_rows(
    results_ratemax$model_fit,
    results_topt$model_fit,
    results_ctmax$model_fit,
    results_ctmin$model_fit
  )
  
  ## Combine random effects into a single dataframe
  random_effects_clean_rate <- phylogenetic_signal_thermal_traits_rate %>%
    mutate(
      metric = "rate",
      Variable = Factor,
      Effect = "random"
    ) %>%
    rename(
      Variance_observed = observed,
      Lower_CI = lower_ci,
      Upper_CI = upper_ci,
      Significant = significant
    ) %>%
    dplyr::select(
      Trait,
      metric,
      Variable,
      Effect,
      Variance_observed,
      Lower_CI,
      Upper_CI,
      Significant
    )
  
  ## Combine fixed effects into a single df too
  fixed_effects_clean_rate <- fixed_effects_all_rate %>%
    mutate(
      metric = "rate",
      Variable = grp,
      Effect = "fixed",
      p_value = `Pr(>|t|)`
    ) %>%
    rename(
      Estimate = Estimate,
      Std_Error = `Std. Error`,
      t_value = `t value`,
      Partial_R2 = share
    ) %>%
    dplyr::select(
      Trait,
      metric,
      Variable,
      Effect,
      Estimate,
      Std_Error,
      t_value,
      p_value,
      Partial_R2
    )
  
  ## Select relevant columns from the global fit df
  model_fit_clean_rate <- model_fit_all %>%
    mutate(
      Lower_CI = NA_real_,
      Upper_CI = NA_real_,
      Significant = NA
    ) %>%
    dplyr::select(
      Trait,
      metric,
      Variable,
      Effect,
      Variance_observed,
      Lower_CI,
      Upper_CI,
      Significant
    )
  
  ## Put everything together into Supplementary Table 2
  combined_effects_lmer_rate <- bind_rows(
    fixed_effects_clean_rate,
    random_effects_clean_rate,
    model_fit_clean_rate
  ) %>%
    mutate(
      across(where(is.numeric), ~ round(.x, 3))
    )
  
  return(
    list(
      combined_effects = combined_effects_lmer_rate,
      phylo_signal = phylogenetic_signal_thermal_traits_rate,
      fixed_effects_all_rate = fixed_effects_all_rate,
      model_fit = model_fit_all
    )
  )
}

# Running pipeline
combined_effects_lmer_rate <-
  phylogenetic_signal_pipeline_rate(
    ave_params_gcplyr_combined_major_taxa_rate
  )

write.csv(combined_effects_lmer_rate$combined_effects, "06_TCS_figures/03_Output/combined_effects_lmer_rate.csv", row.names = FALSE)
write.csv(combined_effects_lmer_rate$phylo_signal,"06_TCS_figures/03_Output/phylogenetic_signal_thermal_traits_rate.csv",row.names = FALSE)
write.csv(combined_effects_lmer_rate$fixed_effects_all_rate,"06_TCS_figures/03_Output/fixed_effects_lmer_rate.csv",row.names = FALSE)


phylogenetic_signal_thermal_traits_rate <- phylogenetic_signal_thermal_traits_rate %>%
  dplyr::mutate(Factor = dplyr::recode(Factor,
                         "Genus:Order" = "Genus"))

phylogenetic_signal_thermal_traits_rate$Factor <- factor(phylogenetic_signal_thermal_traits_rate$Factor, 
                                                        levels = c("Genus","Order","Residual","Site"))

fixed_effects_all_rate <- fread("06_TCS_figures/03_Output/fixed_effects_lmer_rate.csv")

#### S23.right ####

S23.right <- fixed_effects_all_rate %>%
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




phylogenetic_signal_thermal_traits_rate <- fread("06_TCS_figures/03_Output/phylogenetic_signal_thermal_traits_rate.csv")



phylogenetic_signal_thermal_traits_rate$Trait <- factor(phylogenetic_signal_thermal_traits_rate$Trait, 
                                                       levels = c("ctmin","ctmax","topt","rmax"))


phylogenetic_signal_thermal_traits_rate$Factor <- factor(phylogenetic_signal_thermal_traits_rate$Factor, 
                                                                 levels = c("Genus","Order","Residual","Site"))


S23B.right <- phylogenetic_signal_thermal_traits_rate %>% 
  #filter(!`Taxonomic Level` == "Family") %>% 
  #mutate(Trait = fct_relevel(Trait,"ctmin_rate","ctmax_rate","topt_rate","rmax_rate")) %>%
  #mutate(Trait = fct_relevel(Trait,"ctmin_rate","ctmax_rate","topt_rate","ratemax" )) %>%
  ggplot(aes(x = Trait, y = observed)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), linewidth = 0.5,width = 0.2, color = "gray50") +
  geom_point(aes(color = significant), size = 6) +
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  ylim(0,1)+
  scale_y_continuous(breaks = seq(0, 1, 0.5))+
  coord_flip()+
  facet_wrap( ~Factor,ncol = 4)+
  labs(
    y = "Proportion of Variance Explained",
    color = "Significant"
  )  +
  scale_color_manual(values = c("grey60","black"))+
  theme(strip.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "top",
        legend.direction = "horizontal"
  )+
  guides(color = "none")


S23B.right

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


#### S7.top.left - Correlation of thermal traits and bioclimatic variables for auc ####
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

f_auc_rcorr <- rcorr(
  as.matrix(correlation_matrix_auc_3),
  type = "spearman"
)


S7.top.left <- ggcorrplot(f_auc_rcorr$r, hc.order = FALSE, type = "upper",
                                                     outline.col = "white",
                                                     ggtheme = ggplot2::theme_classic(
                                                       base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                                                     colors = c("#6D9EC1", "white", "#E46726"),
                                                     lab = TRUE,
                                                     insig = "pch",
                           p.mat = f_auc_rcorr$P,
                           sig.level = 0.05,
                           pch = 0,
                           pch.col = "red3",
                           pch.cex = 12.5,
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

S7.top.left


#### S7.top.center and right - Correlation of thermal traits and bioclimatic variables for auc by Order ####

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


library("Hmisc")


f_auc_Enterobacterales_rcorr <- rcorr(
  as.matrix(correlation_matrix_auc_3_Enterobacterales),
  type = "spearman"
)


# Get the upper triangle
S7.top.center <- ggcorrplot(f_auc_Enterobacterales_rcorr$r, hc.order = FALSE, 
                             type = "upper",
                                                     outline.col = "white",
                                                     ggtheme = ggplot2::theme_classic(
                                                       base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                                                     colors = c("#6D9EC1", "white", "#E46726"),
                                                     lab = TRUE,
                                                     insig = "pch",
                             p.mat = f_auc_Enterobacterales_rcorr$P,
                             sig.level = 0.05,
                             pch = 0,
                             pch.col = "red3",
                             pch.cex = 12.5,
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

f_auc_Pseudomonadales_rcorr <- rcorr(
  as.matrix(correlation_matrix_auc_3_Pseudomonadales),
  type = "spearman"
)


# Get the upper triangle
S7.top.right <- ggcorrplot(f_auc_Pseudomonadales_rcorr$r, hc.order = FALSE, type = "upper",
                                                                      outline.col = "white",
                                                                      ggtheme = ggplot2::theme_classic(
                                                                        base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                                                                      colors = c("#6D9EC1", "white", "#E46726"),
                                                                      lab = TRUE,
                                                                      insig = "pch",
                            p.mat = f_auc_Pseudomonadales_rcorr$P,
                            sig.level = 0.05,
                            pch = 0,
                            pch.col = "red3",
                            pch.cex = 12.5,
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

S7.top.left
S7.top.center
S7.top.right

#### Correlation matrix for AUC ####
correlation_matrix_auc <- ave_params_gcplyr_combined_major_taxa_auc %>% 
  dplyr::select(c("seqID","Order","Family","Genus","temperature_of_isolation","Site_initials","rmax",     
                  "topt" ,"ctmin","ctmax","thermal_safety_margin","thermal_tolerance","PC1","PC2"))



correlation_matrix_auc_2 <- correlation_matrix_auc %>% dplyr::select(c("temperature_of_isolation","Order","Genus","rmax", "topt" ,"ctmin","ctmax","thermal_safety_margin","thermal_tolerance",
                                                                       "PC1","PC2","Site_initials"))


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

model_lmer_auc_rmax <- lmer(log10(rmax) ~ PC1 + PC2 + Order + temperature_of_isolation + (1 | Genus/Order) + (1|Site_initials), data = lmer_auc_rmax)
plot(model_lmer_auc_rmax)              # residuals vs fitted
qqnorm(resid(model_lmer_auc_rmax))     # normality of residuals
qqline(resid(model_lmer_auc_rmax))
summary(model_lmer_auc_rmax)
                     

lmer_auc_topt <- correlation_matrix_auc_2 %>% dplyr::select(!c(rmax, thermal_tolerance,thermal_safety_margin, ctmin, ctmax))

model_lmer_auc_topt <- lmer(topt ~ PC1 + PC2 + Order  + temperature_of_isolation + (1 | Genus/Order) + (1|Site_initials), data = lmer_auc_topt)
plot(model_lmer_auc_topt)              # residuals vs fitted
qqnorm(resid(model_lmer_auc_topt))     # normality of residuals
qqline(resid(model_lmer_auc_topt))
summary(model_lmer_auc_topt)


lmer_auc_ctmax <- correlation_matrix_auc_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmin, rmax))

model_lmer_auc_ctmax <- lmer(ctmax ~ PC1 + PC2 + Order + temperature_of_isolation + (1 | Genus/Order) + (1|Site_initials), data = lmer_auc_ctmax)
plot(model_lmer_auc_ctmax)              # residuals vs fitted
qqnorm(resid(model_lmer_auc_ctmax))     # normality of residuals
qqline(resid(model_lmer_auc_ctmax))
summary(model_lmer_auc_ctmax)


lmer_auc_ctmin <- correlation_matrix_auc_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, rmax, ctmax))

model_lmer_auc_ctmin <- lmer(ctmin ~ PC1 + PC2 + Order + temperature_of_isolation + (1 | Genus/Order) + (1|Site_initials), data = lmer_auc_ctmin)
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


partR2(model_lmer_auc_rmax, partvars = c("PC1", "PC2", "Order", "temperature_of_isolation"))
partR2(model_lmer_auc_topt, partvars = c("PC1", "PC2", "Order", "temperature_of_isolation"))
partR2(model_lmer_auc_ctmax, partvars = c("PC1", "PC2", "Order", "temperature_of_isolation"))
partR2(model_lmer_auc_ctmin, partvars = c("PC1", "PC2", "Order", "temperature_of_isolation"))


#### S7.bot.left - Correlation of thermal traits and bioclimatic variables for rate ####
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

f_rate_rcorr <- rcorr(
  as.matrix(correlation_matrix_rate_3),
  type = "spearman"
)


library(ggcorrplot)

S7.bot.left <- ggcorrplot(f_rate_rcorr$r, hc.order = FALSE, type = "upper",
                           outline.col = "white",
                           ggtheme = ggplot2::theme_classic(
                             base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                           colors = c("#6D9EC1", "white", "#E46726"),
                           lab = TRUE,
                           insig = "pch",
                           p.mat = f_rate_rcorr$P,
                           sig.level = 0.05,
                           pch = 0,
                           pch.col = "red3",
                           pch.cex = 12.5,
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

S7.bot.left


#### S7.bot.center and right - Correlation of thermal traits and bioclimatic variables for rate by Order ####

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


library("Hmisc")


f_rate_Enterobacterales_rcorr <- rcorr(
  as.matrix(correlation_matrix_rate_3_Enterobacterales),
  type = "spearman"
)


# Get the upper triangle
S7.bot.center <- ggcorrplot(f_rate_Enterobacterales_rcorr$r, hc.order = FALSE, 
                             type = "upper",
                             outline.col = "white",
                             ggtheme = ggplot2::theme_classic(
                               base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                             colors = c("#6D9EC1", "white", "#E46726"),
                             lab = TRUE,
                             insig = "pch",
                             p.mat = f_rate_Enterobacterales_rcorr$P,
                             sig.level = 0.05,
                             pch = 0,
                             pch.col = "red3",
                             pch.cex = 12.5,
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

f_rate_Pseudomonadales_rcorr <- rcorr(
  as.matrix(correlation_matrix_rate_3_Pseudomonadales),
  type = "spearman"
)


# Get the upper triangle
S7.bot.right <- ggcorrplot(f_rate_Pseudomonadales_rcorr$r, hc.order = FALSE, type = "upper",
                            outline.col = "white",
                            ggtheme = ggplot2::theme_classic(
                              base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light"),
                            colors = c("#6D9EC1", "white", "#E46726"),
                            lab = TRUE,
                            insig = "pch",
                            p.mat = f_rate_Pseudomonadales_rcorr$P,
                            sig.level = 0.05,
                            pch = 0,
                            pch.col = "red3",
                            pch.cex = 12.5,
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

S7.bot.left
S7.bot.center
S7.bot.right


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

lmer_rate_output <- data.frame(
  Trait = character(),
  PC1_mean = numeric(),
  PC1_sd = numeric(),
  PC2_mean = numeric(),
  PC2_sd = numeric(),
  stringsAsFactors = FALSE
)


# Then we run a model for each trait. For rmax we first transform its distribution to a more normalized shape using log10. 
# This transformation gets us closer to a normal shape, but the Shapiro-Wilk test is very sensitive, so it still fails

lmer_rate_rmax <- correlation_matrix_rate_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmin, ctmax))

model_lmer_rate_rmax <- lmer(log10(rmax) ~ PC1 + PC2 + (1 | Order/Genus), data = lmer_rate_rmax)
plot(model_lmer_rate_rmax)              # residuals vs fitted
qqnorm(resid(model_lmer_rate_rmax))     # normality of residuals
qqline(resid(model_lmer_rate_rmax))
summary(model_lmer_rate_rmax)


lmer_rate_topt <- correlation_matrix_rate_2 %>% dplyr::select(!c(rmax, thermal_tolerance,thermal_safety_margin, ctmin, ctmax))

model_lmer_rate_topt <- lmer(topt ~ PC1 + PC2 + (1 | Order/Genus), data = lmer_rate_topt)
plot(model_lmer_rate_topt)              # residuals vs fitted
qqnorm(resid(model_lmer_rate_topt))     # normality of residuals
qqline(resid(model_lmer_rate_topt))
summary(model_lmer_rate_topt)


lmer_rate_ctmax <- correlation_matrix_rate_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, ctmin, rmax))

model_lmer_rate_ctmax <- lmer(ctmax ~ PC1 + PC2 + (1 | Order/Genus), data = lmer_rate_ctmax)
plot(model_lmer_rate_ctmax)              # residuals vs fitted
qqnorm(resid(model_lmer_rate_ctmax))     # normality of residuals
qqline(resid(model_lmer_rate_ctmax))
summary(model_lmer_rate_ctmax)


lmer_rate_ctmin <- correlation_matrix_rate_2 %>% dplyr::select(!c(topt, thermal_tolerance,thermal_safety_margin, rmax, ctmax))

model_lmer_rate_ctmin <- lmer(ctmin ~ PC1 + PC2 + (1 | Order/Genus), data = lmer_rate_ctmin)
plot(model_lmer_rate_ctmin)              # residuals vs fitted
qqnorm(resid(model_lmer_rate_ctmin))     # normality of residuals
qqline(resid(model_lmer_rate_ctmin))
summary(model_lmer_rate_ctmin)

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
lmer_rate_ctmin_coefficients <- extract_coefficients_lmer(model_lmer_rate_ctmin, "CTmin")
lmer_rate_ctmax_coefficients <- extract_coefficients_lmer(model_lmer_rate_ctmax, "CTmax")
lmer_rate_topt_coefficients  <- extract_coefficients_lmer(model_lmer_rate_topt,  "Topt")
lmer_rate_rmax_coefficients  <- extract_coefficients_lmer(model_lmer_rate_rmax,  "rmax")


# Combine into a data frame
lmer_rate_output <- data.frame(bind_rows(lmer_rate_rmax_coefficients, 
                                       lmer_rate_topt_coefficients, 
                                       lmer_rate_ctmax_coefficients, 
                                       lmer_rate_ctmin_coefficients), 
                             stringsAsFactors = FALSE)

# Convert numeric columns from character (after cbind) back to numeric
lmer_rate_output[, 2:5] <- lapply(lmer_rate_output[, 2:5], as.numeric)






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


#### S9 - Plotting the TPC triangles by Climate Profile ####

# Here we need to obtain the triangle shapes by mapping the coordinates 
# of each vertex to x and y columns of a new dataframe 

full_triangle_data <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,Family,Genus,`Climate profile`,topt,rmax,ctmin,ctmax,metric) %>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(seqID,metric) %>%
  dplyr::summarize(
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
  dplyr::summarize(
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
  dplyr::summarize(
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
  dplyr::summarize(sd_topt = sd(topt),
            sd_max = sd(rmax),
            sd_ctmax = sd(ctmax),
            sd_ctmin = sd(ctmin),
            topt = mean(topt),
            rmax = mean(rmax),
            ctmax = mean(ctmax),
            ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(`Climate profile`,Order,metric) %>%
  dplyr::summarize(
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
  dplyr::summarize(sd_topt = sd(topt),
            sd_max = sd(rmax),
            sd_ctmax = sd(ctmax),
            sd_ctmin = sd(ctmin),
            topt = mean(topt),
            rmax = mean(rmax),
            ctmax = mean(ctmax),
            ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(`Climate profile`,Order,metric) %>%
  dplyr::summarize(
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
  dplyr::summarize(sd_topt = sd(topt),
            sd_max = sd(rmax),
            sd_ctmax = sd(ctmax),
            sd_ctmin = sd(ctmin),
            topt = mean(topt),
            rmax = mean(rmax),
            ctmax = mean(ctmax),
            ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(`Climate profile`,Order,metric) %>%
  dplyr::summarize(
    x = list(c(ctmin,topt,topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order),
    `Climate profile` = first(`Climate profile`)
  ) %>%
  unnest(cols = c(x,y)) %>%  # Unnest to get a tidy dataframe for plotting
  mutate(group = interaction(Order,metric))  


full_triangle_data$`Climate profile` <- factor(full_triangle_data$`Climate profile`)


S9.left <- full_triangle_data %>% 
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
  geom_errorbar(data = full_triangle_data_Order %>%
                  filter(metric == "auc")%>%
                  filter(y == 0),
                width = 0.7,
                color = "grey40",
                lwd = 1,
                aes(x = x, y = y,
                    xmin = x - z, xmax = x + z))+
  geom_errorbar(data = full_triangle_data_Order %>%
                  filter(metric == "auc")%>%
                  filter(y != 0),
                width = 3,
                lwd = 1,
                color = "grey40",
                aes(x = x, y = y,
                    ymin = y - z, ymax = y + z)
  )+
  geom_point(data = full_triangle_data_Order %>% 
               filter(metric == "auc"),
             shape = 21,color = "grey40",aes(x = x, y = y),
             size = 2,
             stroke = 2.25)+
  geom_point(data = full_triangle_data_Order %>% 
               filter(metric == "auc"),
             shape = 21,color = "white",aes(x = x, y = y),size = 2)+
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

S9.right <- full_triangle_data %>% 
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
  geom_errorbar(data = full_triangle_data_Order %>%
                  filter(metric == "rate")%>%
                  filter(y == 0),
                width = 0.1,
                color = "grey40",
                lwd = 1,
                aes(x = x, y = y,
                    xmin = x - z, xmax = x + z))+
  geom_errorbar(data = full_triangle_data_Order %>%
                  filter(metric == "rate")%>%
                  filter(y != 0),
                width = 3,
                lwd = 1,
                color = "grey40",
                aes(x = x, y = y,
                    ymin = y - z, ymax = y + z)
  )+
  geom_point(data = full_triangle_data_Order %>% 
               filter(metric == "rate"),
             shape = 21,color = "grey40",aes(x = x, y = y),
             size = 2,
             stroke = 2.25)+
  geom_point(data = full_triangle_data_Order %>% 
               filter(metric == "rate"),
             shape = 21,color = "white",aes(x = x, y = y),size = 2)+
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


#### F4A - Plotting the TPC triangles ####

# Here we need to obtain the triangle shapes by mapping the coordinates 
# of each vertex to x and y columns of a new dataframe 

full_triangle_data_noenv <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order,Family,Genus,topt,rmax,ctmin,ctmax,metric) %>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(seqID,metric) %>%
  dplyr::summarize(
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
  dplyr::summarize(
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
  dplyr::summarize(
    x = list(c(ctmin, topt, topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus)
  ) %>%
  unnest(cols = c(x, y)) %>%  
  mutate(group = interaction(seqID,metric))  

variance_thermal_traits_by_Order <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order, rmax, topt,ctmax,ctmin,metric) %>%
  group_by(Order,metric) %>%
  dplyr::summarize(sd_topt = sd(topt),
                   sd_max = sd(rmax),
                   sd_ctmax = sd(ctmax),
                   sd_ctmin = sd(ctmin))

full_triangle_data_noenv_Order <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order, rmax, topt,ctmax,ctmin,metric) %>%
  group_by(Order,metric) %>%
  dplyr::summarize(sd_topt = sd(topt)/mean(topt),
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
  dplyr::summarize(
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
  dplyr::summarize(sd_topt = sd(topt)/mean(topt),
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
  dplyr::summarize(
    x = list(c(ctmax,topt,topt)),
    y = list(c(0,0,rmax)),
    Order = first(Order)
  ) %>%
  unnest(cols = c(x,y)) %>%  # Unnest to get a tidy dataframe for plotting
  mutate(group = interaction(Order,metric))  

lower_triangle_data_noenv_Order <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order, rmax, topt, ctmax,ctmin,metric) %>%
  group_by(Order,metric) %>%
  dplyr::summarize(sd_topt = sd(topt)/mean(topt),
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
  dplyr::summarize(
    x = list(c(ctmin,topt,topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order)
  ) %>%
  unnest(cols = c(x,y)) %>%  # Unnest to get a tidy dataframe for plotting
  mutate(group = interaction(Order,metric))  




### Using absolute variances

full_triangle_data_noenv_Order <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order, rmax, topt,ctmax,ctmin,metric) %>%
  group_by(Order,metric) %>%
  dplyr::summarize(sd_topt = sd(topt),
                   sd_max = sd(rmax),
                   sd_ctmax = sd(ctmax),
                   sd_ctmin = sd(ctmin),
                   topt = mean(topt),
                   rmax = mean(rmax),
                   ctmax = mean(ctmax),
                   ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(Order,metric) %>%
  dplyr::summarize(
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
  dplyr::summarize(sd_topt = sd(topt),
                   sd_max = sd(rmax),
                   sd_ctmax = sd(ctmax),
                   sd_ctmin = sd(ctmin),
                   topt = mean(topt),
                   rmax = mean(rmax),
                   ctmax = mean(ctmax),
                   ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(Order,metric) %>%
  dplyr::summarize(
    x = list(c(ctmax,topt,topt)),
    y = list(c(0,0,rmax)),
    Order = first(Order)
  ) %>%
  unnest(cols = c(x,y)) %>%  # Unnest to get a tidy dataframe for plotting
  mutate(group = interaction(Order,metric))  

lower_triangle_data_noenv_Order <- ave_params_gcplyr_combined_major_taxa %>%
  dplyr::select(seqID, Order, rmax, topt, ctmax,ctmin,metric) %>%
  group_by(Order,metric) %>%
  dplyr::summarize(sd_topt = sd(topt),
                   sd_max = sd(rmax),
                   sd_ctmax = sd(ctmax),
                   sd_ctmin = sd(ctmin),
                   topt = mean(topt),
                   rmax = mean(rmax),
                   ctmax = mean(ctmax),
                   ctmin = mean(ctmin))%>%
  ungroup()%>%
  distinct() %>%  # Remove duplicates if necessary
  group_by(Order,metric) %>%
  dplyr::summarize(
    x = list(c(ctmin,topt,topt)),
    y = list(c(0,0, rmax)),
    Order = first(Order)
  ) %>%
  unnest(cols = c(x,y)) %>%  # Unnest to get a tidy dataframe for plotting
  mutate(group = interaction(Order,metric))

F4A.auc_noenv <- full_triangle_data_noenv %>% 
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

F4A.rate_noenv <- full_triangle_data_noenv %>% 
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
  geom_errorbar(data = full_triangle_data_noenv_Order %>%
                  filter(metric == "rate")%>%
                  filter(y == 0),
                width = 0.03,
                color = "grey40",
                lwd = 1.25,
                aes(x = x, y = y,
                    xmin = x - z, xmax = x + z))+
  geom_errorbar(data = full_triangle_data_noenv_Order %>%
                  filter(metric == "rate")%>%
                  filter(y != 0),
                width = 5,
                lwd = 1.25,
                color = "grey40",
                aes(x = x, y = y,
                    ymin = y - z, ymax = y + z)
  )+
  geom_point(data = full_triangle_data_noenv_Order %>% 
               filter(metric == "rate"),
             shape = 21,color = "grey40",aes(x = x, y = y),
             size = 4,
             stroke = 2.25)+
  geom_point(data = full_triangle_data_noenv_Order %>% 
               filter(metric == "rate"),
             shape = 21,color = "white",aes(x = x, y = y),size = 4)+
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = "rate"
  ) +
  facet_wrap(~Order,nrow= 1)+
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


#### F4C (Topt - MAT) ####

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

ave_params_gcplyr_combined_major_taxa_auc$Climate_profile <- ave_params_gcplyr_combined_major_taxa_auc$`Climate profile`

anova_model_topt_mat <- aov((topt - MAT) ~ Climate_profile, data = ave_params_gcplyr_combined_major_taxa_auc)
summary(anova_model_topt_mat)
posthoc_topt_mat <- TukeyHSD(anova_model_topt_mat, which = )
posthoc_topt_mat

F4C.auc <- ave_params_gcplyr_combined_major_taxa %>%
  filter(metric == "auc")%>%
  ggplot(aes(x=`Climate profile`, y=(topt - MAT)))+
  #scale_color_viridis_b()+
  geom_hline(yintercept = 0,linetype = "dashed",color = "grey60",linewidth = 0.75)+
  geom_jitter(alpha = 0.8,width = 0.2,size= 5,aes(color = Order))+
  geom_boxplot(alpha = 0.8,outliers = FALSE)+
  #facet_wrap(~Order,ncol = 1)+
  ylim(-10,30)+
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

S24C <- ave_params_gcplyr_combined_major_taxa %>%
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

#### F4B top - Randomizing thermal traits for AUC and plotting slope above Topt ####

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
  dplyr::summarize(superopt_slope = max(superopt_slope))

model_max_simulated_slopes <- lm(superopt_slope ~ rmax,data = simulated_slopes_auc_df_max)
summary(model_max_simulated_slopes)

model_simulated_slopes <- lm(superopt_slope ~ rmax,data = simulated_slopes_auc_df)
summary(model_simulated_slopes)

F4Btop.auc <- ave_params_gcplyr_combined_major_taxa_auc %>%
  ggplot(aes(x=rmax, y= superopt_slope))+
  geom_point(data = simulated_slopes_auc_df, aes(x = rmax, y = superopt_slope),color = "grey85")+
  #scale_color_viridis_b()+
  geom_smooth(method = "lm", se = TRUE, fill = "grey60", color="grey90") +
  geom_point(alpha = 1,aes(color = Order, size = ctmax))+
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
  #scale_fill_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)
  )

#### F4B bot - Comparing slope ranks AUC ####

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


F4Bbot.auc <- slope_data_ranks_full_2 %>% 
  filter(superopt_slope < 0) %>%
  filter(sim_id == 0) %>%
  ggplot(aes(x = rmax, y = slope_ratio_mean,color = Order)) +
  #geom_smooth()+
  ylim(0,2.5)+
  geom_hline(yintercept = 1,linetype = "dashed",linewidth = 0.5)+
  #geom_point(size = 4,alpha = 1)+
  geom_point(size = 4)+
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

F4B.auc <- plot_grid(F4Btop.auc + theme(legend.position = "none"), 
                     F4Bbot.auc, ncol = 1, align = "v", rel_heights = c(1,1),axis = "l")


#### S24B top - Randomizing thermal traits for rate and plotting slope above Topt ####

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

S24B.top <- ave_params_gcplyr_combined_major_taxa_rate %>%
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

#### S24B bot - Comparing slope ranks rate ####

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


S24B.bot <- slope_data_ranks_full_2 %>% 
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

S24B <- plot_grid(S24B.top + theme(legend.position = "none"), 
                  S24B.bot, ncol = 1, align = "v", rel_heights = c(1,1),
                  axis = "l")
S24B

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
p_value_Welch_auc <- 2 * pt(-abs(t_stat_Welch_auc), Welch_auc)

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
p_value_Welch_rate <- 2 * pt(-abs(t_stat_Welch_rate), Welch_rate)

# Display the test statistic and p-value
t_stat_Welch_rate
p_value_Welch_rate


#### Processing data to plot F5 by Climate Profile AUC ####
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


# Scale the datasets to a (0,1) range
continuous_proportion_data_growing_season_site <- continuous_proportion_data_growing_season_site %>%
  mutate(Proportion_scaled = (Proportion_continuous - min(Proportion_continuous)) /
           (max(Proportion_continuous) - min(Proportion_continuous)))

temperature_probability_df <- continuous_proportion_data_growing_season_site %>%
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
  dplyr::summarize(
    Fitness = mean(Fitness),
    Fitted_scaled = mean(Fitted_scaled)
  )%>%
  filter(Original_climate == New_climate)


# We integrate it by curve and Climate profile

Fitness_by_curve_climate_auc_summarized <- Fitness_by_curve_climate_auc %>%
  group_by(seqID,Order,Original_climate,New_climate,Original_site,New_site) %>%
  dplyr::summarize("Integrated Fitness" = gcplyr::auc(Temperature,Fitness)) 


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
  dplyr::summarize(sd_Native_fitness = sd(Native_fitness),
            Native_fitness = mean(Native_fitness))%>% 
  dplyr::select(New_site,Native_fitness) 

Fitness_by_curve_climate_auc_summarized_Entero <- Fitness_by_curve_climate_auc_summarized_Entero %>%
  left_join(native_fitness, by = "New_site")

Fitness_by_curve_climate_auc_summarized_Entero <- Fitness_by_curve_climate_auc_summarized_Entero %>%
  group_by(seqID,Original_site,New_site,Original_climate,New_climate)%>%
  dplyr::summarize(delta_A_H = `Integrated Fitness` - Original_fitness,
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
  dplyr::summarize(sd_Native_fitness = sd(Native_fitness),
            Native_fitness = mean(Native_fitness))%>% 
  dplyr::select(New_site,Native_fitness) 

Fitness_by_curve_climate_auc_summarized_Pseudo <- Fitness_by_curve_climate_auc_summarized_Pseudo %>%
  left_join(native_fitness, by = "New_site")

Fitness_by_curve_climate_auc_summarized_Pseudo <- Fitness_by_curve_climate_auc_summarized_Pseudo %>%
  group_by(seqID,Original_site,New_site,Original_climate,New_climate)%>%
  dplyr::summarize(delta_A_H = `Integrated Fitness` - Original_fitness,
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
  dplyr::summarize(sd_delta_A_H = sd(delta_A_H),
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
  dplyr::summarize(sd_delta_A_H = sd(delta_A_H),
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

F5A1A.auc <- ggplot() +
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


F5A1B.auc <- ggplot() +
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


F5A1C.auc <- ggplot() +
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

F5A1D.auc <- ggplot() +
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



F5A2B.auc <- Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Cold variable")) %>%
  filter(Original_climate == New_climate)%>%
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

F5A2D.auc <- Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Hot variable")) %>%
  filter(Original_climate == New_climate)%>%
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

F5A2A.auc <- Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Cold stable")) %>%
  filter(Original_climate == New_climate)%>%
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

F5A2C.auc <- Fitness_by_curve_climate_auc_summarized %>% 
  filter(Original_climate %in% c("Hot stable")) %>%
  filter(Original_climate == New_climate)%>%
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

F5AA.auc<- F5A1A.auc+ 
  annotation_custom(ggplotGrob(F5A2A.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5AB.auc <- F5A1B.auc + 
  annotation_custom(ggplotGrob(F5A2B.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5AC.auc <- F5A1C.auc + 
  annotation_custom(ggplotGrob(F5A2C.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5AD.auc <- F5A1D.auc + 
  annotation_custom(ggplotGrob(F5A2D.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5A.auc.site <- plot_grid(F5AB.auc,F5AD.auc,F5AA.auc,F5AC.auc, 
                     ncol = 2,
                     rel_widths = c(1,1,1,1),
                     align = "hv",
                     axis = "l")
F5A.auc.site



#### F5 - Plotting the temperature distribution by Climate profile, with the TPCs overlaid AUC ####

# Scale the datasets to a (0,1) range
continuous_proportion_data_growing_season <- continuous_proportion_data_growing_season %>%
  mutate(Proportion_scaled = (Proportion_continuous - min(Proportion_continuous)) /
           (max(Proportion_continuous) - min(Proportion_continuous)))

F5A1A.auc <- ggplot() +
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


F5A1B.auc <- ggplot() +
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


F5A1C.auc <- ggplot() +
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

F5A1D.auc <- ggplot() +
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

F5AA.auc<- F5A1A.auc+ 
  annotation_custom(ggplotGrob(F5A2A.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5AB.auc <- F5A1B.auc + 
  annotation_custom(ggplotGrob(F5A2B.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5AC.auc <- F5A1C.auc + 
  annotation_custom(ggplotGrob(F5A2C.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5AD.auc <- F5A1D.auc + 
  annotation_custom(ggplotGrob(F5A2D.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5A.auc <- plot_grid(F5AB.auc,F5AD.auc,F5AA.auc,F5AC.auc, 
                     ncol = 2,
                     rel_widths = c(1,1,1,1),
                     align = "hv",
                     axis = "l")
F5A.auc


#### S14 - Plotting the yearly temperature distribution by Climate profile, with the TPCs overlaid AUC ####
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
  dplyr::summarize(
    Fitness = mean(Fitness),
    Fitted_scaled = mean(Fitted_scaled)
  )%>%
  filter(Original_climate == New_climate)


# We integrate it by curve and Climate profile

yearly_Fitness_by_curve_climate_auc_summarized <- yearly_Fitness_by_curve_climate_auc %>%
  group_by(seqID,Order,Original_climate,New_climate) %>%
  dplyr::summarize("Integrated Fitness" = gcplyr::auc(Temperature,Fitness)) 


S14A2B.auc <- yearly_Fitness_by_curve_climate_auc_summarized %>% 
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

S14A2D.auc <- yearly_Fitness_by_curve_climate_auc_summarized %>% 
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

S14A2A.auc <- yearly_Fitness_by_curve_climate_auc_summarized %>% 
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

S14A2C.auc <- yearly_Fitness_by_curve_climate_auc_summarized %>% 
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


S14A1A.auc <- ggplot() +
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
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


S14A1B.auc <- ggplot() +
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
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 


S14A1C.auc <- ggplot() +
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
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

S14A1D.auc <- ggplot() +
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
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  ) 

S14AA.auc<- S14A1A.auc+ 
  annotation_custom(ggplotGrob(S14A2A.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

S14AB.auc <- S14A1B.auc + 
  annotation_custom(ggplotGrob(S14A2B.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

S14AC.auc <- S14A1C.auc + 
  annotation_custom(ggplotGrob(S14A2C.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

S14AD.auc <- S14A1D.auc + 
  annotation_custom(ggplotGrob(S14A2D.auc),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

S14 <- plot_grid(S14AB.auc,S14AD.auc,S14AA.auc,S14AC.auc, 
                ncol = 2,
                rel_widths = c(1,1,1,1),
                align = "hv",
                axis = "l")
S14


#### Calculating the fraction of time bacteria are exposed to temperatures above Topt AUC ####

Fitness_by_order_climate_auc_topt <- Fitness_by_order_climate_auc %>%
  group_by(Original_climate,Order) %>%
  dplyr::summarize(Topt = Temperature[Fitted_scaled == max(Fitted_scaled)],
            CTmax = Temperature[Fitted_scaled == min(Fitted_scaled [Temperature > Topt])])

growth_season_continuous_proportion_data_auc_topt <- continuous_proportion_data_growing_season %>%
  mutate("Original_climate" = continuous_proportion_data_growing_season$`Climate profile`)%>%
  left_join(Fitness_by_order_climate_auc_topt)

growth_season_supraoptimal_time_auc <- growth_season_continuous_proportion_data_auc_topt %>%
  group_by(Order, `Climate profile`) %>%
  dplyr::summarize(
    Hard_times = (sum(Proportion_continuous[Temp_continuous > Topt])/sum(Proportion_continuous))*100,
    Harder_times = (sum(Proportion_continuous[Temp_continuous > 40])/sum(Proportion_continuous))*100,
    CTmax = first(CTmax),
    Topt = first(Topt)
  )

yearly_continuous_proportion_data_auc_topt <- continuous_proportion_data %>%
  mutate("Original_climate" = continuous_proportion_data$`Climate profile`)%>%
  left_join(Fitness_by_order_climate_auc_topt)

yearly_supraoptimal_time <- yearly_continuous_proportion_data_auc_topt %>%
  group_by(Order, `Climate profile`) %>%
  dplyr::summarize(
    Hard_times = (sum(Proportion_continuous[Temp_continuous > Topt])/sum(Proportion_continuous))*100,
    Harder_times = (sum(Proportion_continuous[Temp_continuous > 40])/sum(Proportion_continuous))*100,
    CTmax = first(CTmax),
    Topt = first(Topt)
  )

#### Processing data to plot S25 ####
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
continuous_proportion_data_growing_season <- continuous_proportion_data_growing_season %>%
  mutate(Proportion_scaled = (Proportion_continuous - min(Proportion_continuous)) /
           (max(Proportion_continuous) - min(Proportion_continuous)))

# Scale the datasets to a (0,1) range
continuous_proportion_data_growing_season_site <- continuous_proportion_data_growing_season_site %>%
  mutate(Proportion_scaled = (Proportion_continuous - min(Proportion_continuous)) /
           (max(Proportion_continuous) - min(Proportion_continuous)))

temperature_probability_df <- continuous_proportion_data_growing_season_site %>%
  rename("Temperature" = "Temp_continuous")

#### Summarizing data by taxon, climate profile, and curve rate ####

# Summarize by curve
means_by_curve_climate_rate <- ave_preds_gcplyr_combined_major_taxa %>% filter(metric == "rate")

means_by_curve_climate_rate <- means_by_curve_climate_rate %>%
  #group_by(`Climate profile`) %>%
  mutate(Fitted_scaled = (.fitted - min(.fitted)) / (max(.fitted) - min(.fitted)))%>%
  ungroup()

means_by_curve_climate_rate <- means_by_curve_climate_rate %>%
  rename("Temperature" = "temp")

#### Unifying temperature scales with curve df rate ######
# The temperature scales between trait data and temperature time series differ, and 
# thus we need to harmonize them

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
  dplyr::summarize(
    Fitness = mean(Fitness),
    Fitted_scaled = mean(Fitted_scaled)
  )%>%
  filter(Original_climate == New_climate)


# We integrate it by curve and Climate profile

Fitness_by_curve_climate_rate_summarized <- Fitness_by_curve_climate_rate %>%
  group_by(seqID,Order,Original_climate,New_climate,Original_site,New_site) %>%
  dplyr::summarize("Integrated Fitness" = gcplyr::auc(Temperature,Fitness)) 


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
  dplyr::summarize(sd_Native_fitness = sd(Native_fitness),
            Native_fitness = mean(Native_fitness))%>% 
  dplyr::select(New_site,Native_fitness) 

Fitness_by_curve_climate_rate_summarized_Entero <- Fitness_by_curve_climate_rate_summarized_Entero %>%
  left_join(native_fitness, by = "New_site")

Fitness_by_curve_climate_rate_summarized_Entero <- Fitness_by_curve_climate_rate_summarized_Entero %>%
  group_by(seqID,Original_site,New_site,Original_climate,New_climate)%>%
  dplyr::summarize(delta_A_H = `Integrated Fitness` - Original_fitness,
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
  dplyr::summarize(sd_Native_fitness = sd(Native_fitness),
            Native_fitness = mean(Native_fitness))%>% 
  dplyr::select(New_site,Native_fitness) 

Fitness_by_curve_climate_rate_summarized_Pseudo <- Fitness_by_curve_climate_rate_summarized_Pseudo %>%
  left_join(native_fitness, by = "New_site")

Fitness_by_curve_climate_rate_summarized_Pseudo <- Fitness_by_curve_climate_rate_summarized_Pseudo %>%
  group_by(seqID,Original_site,New_site,Original_climate,New_climate)%>%
  dplyr::summarize(delta_A_H = `Integrated Fitness` - Original_fitness,
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
  dplyr::summarize(sd_delta_A_H = sd(delta_A_H),
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
  dplyr::summarize(sd_delta_A_H = sd(delta_A_H),
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


F5_by_site_rate <- Fitness_by_curve_climate_rate %>%
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

F5A1A.rate <- ggplot() +
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


F5A1B.rate <- ggplot() +
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


F5A1C.rate <- ggplot() +
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

F5A1D.rate <- ggplot() +
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



F5A2B.rate <- Fitness_by_curve_climate_rate_summarized %>% 
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

F5A2D.rate <- Fitness_by_curve_climate_rate_summarized %>% 
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

F5A2A.rate <- Fitness_by_curve_climate_rate_summarized %>% 
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

F5A2C.rate <- Fitness_by_curve_climate_rate_summarized %>% 
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

F5AA.rate<- F5A1A.rate+ 
  annotation_custom(ggplotGrob(F5A2A.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5AB.rate <- F5A1B.rate + 
  annotation_custom(ggplotGrob(F5A2B.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5AC.rate <- F5A1C.rate + 
  annotation_custom(ggplotGrob(F5A2C.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5AD.rate <- F5A1D.rate + 
  annotation_custom(ggplotGrob(F5A2D.rate),xmin = 37, xmax = 57, 
                    ymin = 0.5, ymax = 1.1)

F5A.rate.site <- plot_grid(F5AB.rate,F5AD.rate,F5AA.rate,F5AC.rate, 
                     ncol = 2,
                     rel_widths = c(1,1,1,1),
                     align = "hv",
                     axis = "l")
F5A.rate.site



#### S25 - Plotting the temperature distribution by Climate profile, with the TPCs overlaid rate ####

S25A1A <- ggplot() +
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


S25A1B <- ggplot() +
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


S25A1C <- ggplot() +
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

S25A1D <- ggplot() +
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

S25 <- plot_grid(S25A1B,S25A1D,S25A1A,S25A1C, 
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

#### S24A - Correlations between TPC traits rate ####
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

S24A <- full_triangle_data_noenv %>% 
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




#### S12 - Coloring slope above Topt by climate profile ####

S12 <- ave_params_gcplyr_combined_major_taxa %>%
  filter(metric == "auc")%>%
  ggplot(aes(x=rmax, y=superopt_slope))+
  #scale_color_viridis_b()+
  #geom_smooth(method = "lm", se = TRUE, fill = "grey40", color="grey90") +
  geom_abline(slope = -0.061272,intercept = c(0,0))+
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

summary(lm(data = ave_params_gcplyr_combined_major_taxa 
           %>% filter(metric == "auc"), 
           superopt_slope ~ rmax))

slope_increase_cs <- lm(data = ave_params_gcplyr_combined_major_taxa 
                                %>% filter(metric == "auc" & 
                                             `Climate profile` == "Cold stable"), 
                                superopt_slope ~ rmax)

linearHypothesis(slope_increase_cs, "rmax = -0.061272")

slope_increase_cv <- lm(data = ave_params_gcplyr_combined_major_taxa 
                        %>% filter(metric == "auc" & 
                                     `Climate profile` == "Cold variable"), 
                        superopt_slope ~ rmax)

linearHypothesis(slope_increase_cv, "rmax = -0.061272")


slope_increase_hs <- lm(data = ave_params_gcplyr_combined_major_taxa 
                        %>% filter(metric == "auc" & 
                                     `Climate profile` == "Hot stable"), 
                        superopt_slope ~ rmax)

linearHypothesis(slope_increase_hs, "rmax = -0.061272")

slope_increase_hv <- lm(data = ave_params_gcplyr_combined_major_taxa 
                        %>% filter(metric == "auc" & 
                                     `Climate profile` == "Hot variable"), 
                        superopt_slope ~ rmax)

linearHypothesis(slope_increase_hv, "rmax = -0.061272")


#### S13 Calculating the variance of Topt and CTmax for bins of AUC ####

regression_ctmax_topt_rmax_Entero <- lm(
  data = ave_params_gcplyr_combined_major_taxa_auc %>%
    filter(Order == "Enterobacterales"), ctmax - topt ~ rmax)
summary(regression_ctmax_topt_rmax_Entero)

regression_ctmax_topt_rmax_Pseudo <- lm(
  data = ave_params_gcplyr_combined_major_taxa_auc %>%
    filter(Order == "Pseudomonadales"), ctmax - topt ~ rmax)
summary(regression_ctmax_topt_rmax_Pseudo)

regression_ctmax_topt_Entero <- lm(
  data = ave_params_gcplyr_combined_major_taxa_auc %>%
    filter(Order == "Enterobacterales"), ctmax - topt ~ ctmax)
summary(regression_ctmax_topt_Entero)

regression_ctmax_topt_Pseudo <- lm(
  data = ave_params_gcplyr_combined_major_taxa_auc %>%
    filter(Order == "Pseudomonadales"), ctmax - topt ~ ctmax)
summary(regression_ctmax_topt_Pseudo)

ave_params_gcplyr_combined_major_taxa_auc_Entero <- 
  ave_params_gcplyr_combined_major_taxa_auc %>%
  filter(Order == "Enterobacterales") %>%
  mutate("Growth quartile" = case_when(
    rmax <= quantile(rmax)[2] ~ "q1",
    rmax <= quantile(rmax)[3] & rmax > quantile(rmax)[2] ~ "q2",
    rmax <= quantile(rmax)[4] & rmax > quantile(rmax)[3] ~ "q3",
    rmax <= quantile(rmax)[5] & rmax > quantile(rmax)[4] ~ "q4",
  )) %>%
    mutate("CTmax quartile" = case_when(
      ctmax <= quantile(ctmax)[2] ~ "q1",
      ctmax <= quantile(ctmax)[3] & ctmax > quantile(ctmax)[2] ~ "q2",
      ctmax <= quantile(ctmax)[4] & ctmax > quantile(ctmax)[3] ~ "q3",
      ctmax <= quantile(ctmax)[5] & ctmax > quantile(ctmax)[4] ~ "q4",
    ))

ave_params_gcplyr_combined_major_taxa_auc_Pseudo <- 
  ave_params_gcplyr_combined_major_taxa_auc %>%
  filter(Order == "Pseudomonadales")%>%
  mutate("Growth quartile" = case_when(
    rmax <= quantile(rmax)[2] ~ "q1",
    rmax <= quantile(rmax)[3] & rmax > quantile(rmax)[2] ~ "q2",
    rmax <= quantile(rmax)[4] & rmax > quantile(rmax)[3] ~ "q3",
    rmax <= quantile(rmax)[5] & rmax > quantile(rmax)[4] ~ "q4",
  ))%>%
    mutate("CTmax quartile" = case_when(
      ctmax <= quantile(ctmax)[2] ~ "q1",
      ctmax <= quantile(ctmax)[3] & ctmax > quantile(ctmax)[2] ~ "q2",
      ctmax <= quantile(ctmax)[4] & ctmax > quantile(ctmax)[3] ~ "q3",
      ctmax <= quantile(ctmax)[5] & ctmax > quantile(ctmax)[4] ~ "q4",
    ))

ave_params_gcplyr_combined_major_taxa_binned <- 
  bind_rows(ave_params_gcplyr_combined_major_taxa_auc_Entero,
            ave_params_gcplyr_combined_major_taxa_auc_Pseudo)


S13A <- ave_params_gcplyr_combined_major_taxa_binned %>%
  ggplot(aes(x = `Growth quartile`, y = ctmax - topt,color = Order))+
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

S13B <- ave_params_gcplyr_combined_major_taxa_binned %>%
  ggplot(aes(x = `CTmax quartile`, y = ctmax - topt,color = Order))+
  geom_boxplot(alpha = 0.6, outliers = FALSE, aes(group = `CTmax quartile`))+
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


ave_params_gcplyr_combined_major_taxa_binned_variance_growth <- ave_params_gcplyr_combined_major_taxa_binned %>%
  group_by(Order, `Growth quartile`) %>%
  dplyr::summarize(sd_ctmax_topt = sd(ctmax - topt))

ave_params_gcplyr_combined_major_taxa_binned_variance_ctmax <- ave_params_gcplyr_combined_major_taxa_binned %>%
  group_by(Order, `CTmax quartile`) %>%
  dplyr::summarize(sd_ctmax_topt = sd(ctmax - topt))

S13C <- ave_params_gcplyr_combined_major_taxa_binned_variance_growth %>%
  ggplot(aes(x = `Growth quartile`, y = sd_ctmax_topt,color = Order))+
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

S13D <- ave_params_gcplyr_combined_major_taxa_binned_variance_ctmax %>%
  ggplot(aes(x = `CTmax quartile`, y = sd_ctmax_topt,color = Order))+
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


S13 <- plot_grid(S13A,S13B,S13C,S13D, ncol = 2,rel_widths = c(1,1,1,1),align = "hv")

S13

#### S3 - Additional assay to evaluate the influence of carbon concentration on TPC shape ####
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
# We have two clean datasets with TPC fits and parameters from another study currently in progress.

ave_params_CCA_gcplyr_combined <- fread("06_TCS_figures/02_Data/250915_ave_params_CCA_gcplyr_combined.csv")
ave_preds_CCA_gcplyr_combined <- fread("06_TCS_figures/02_Data/250915_ave_preds_CCA_gcplyr_combined.csv")

summary(aov(data = ave_params_CCA_gcplyr_combined %>% dplyr::filter(metric == "auc"), topt ~ Glucose_stock_X))
summary(aov(data = ave_params_CCA_gcplyr_combined %>% dplyr::filter(metric == "auc"), ctmax ~ Glucose_stock_X))
summary(aov(data = ave_params_CCA_gcplyr_combined %>% dplyr::filter(metric == "auc"), ctmin ~ Glucose_stock_X))



S3A <- ave_preds_CCA_gcplyr_combined %>%
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

lm_concentration_topt <- lm(topt_100 ~ topt_10, data = ave_params_CCA_gcplyr_combined_wide %>% dplyr::filter(metric == "auc"))
ci_lm_concentration_topt <- confint(lm_concentration_topt, "topt_10")
b_lm_concentration_topt  <- coef(summary(lm_concentration_topt))["topt_10", "Estimate"]
se_lm_concentration_topt <- coef(summary(lm_concentration_topt))["topt_10", "Std. Error"]
t_lm_concentration_topt <- (b_lm_concentration_topt - 1) / se_lm_concentration_topt
p_lm_concentration_topt <- 2 * pt(-abs(t_lm_concentration_topt), df.residual(lm_concentration_topt))

lm_concentration_ctmax<- lm(ctmax_100 ~ ctmax_10, data = ave_params_CCA_gcplyr_combined_wide %>% dplyr::filter(metric == "auc"))
ci_lm_concentration_ctmax <- confint(lm_concentration_ctmax, "ctmax_10")
b_lm_concentration_ctmax  <- coef(summary(lm_concentration_ctmax))["ctmax_10", "Estimate"]
se_lm_concentration_ctmax <- coef(summary(lm_concentration_ctmax))["ctmax_10", "Std. Error"]
t_lm_concentration_ctmax <- (b_lm_concentration_ctmax - 1) / se_lm_concentration_ctmax
p_lm_concentration_ctmax <- 2 * pt(-abs(t_lm_concentration_ctmax), df.residual(lm_concentration_ctmax))

lm_concentration_ctmin <- lm(ctmin_100 ~ ctmin_10, data = ave_params_CCA_gcplyr_combined_wide %>% dplyr::filter(metric == "auc"))
ci_lm_concentration_ctmin <- confint(lm_concentration_ctmin, "ctmin_10")
b_lm_concentration_ctmin  <- coef(summary(lm_concentration_ctmin))["ctmin_10", "Estimate"]
se_lm_concentration_ctmin <- coef(summary(lm_concentration_ctmin))["ctmin_10", "Std. Error"]
t_lm_concentration_ctmin <- (b_lm_concentration_ctmin - 1) / se_lm_concentration_ctmin
p_lm_concentration_ctmin <- 2 * pt(-abs(t_lm_concentration_ctmin), df.residual(lm_concentration_ctmin))

S3B <- ave_params_CCA_gcplyr_combined_wide %>%
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

S3C <- ave_params_CCA_gcplyr_combined_wide %>%
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

S3D <- ave_params_CCA_gcplyr_combined_wide %>%
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

S3E <- ave_params_CCA_gcplyr_combined_wide %>%
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


S3F <- ave_params_CCA_gcplyr_combined_wide %>%
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

S3A.1 <- S3A

S3A.2 <- plot_grid(S3C,S3D,S3E, 
                   ncol = 3,
                   rel_widths = c(1,1,1.2),
                   rel_heights = c(1,1,1),
                   align = "v",
                   axis = "l")
S3 <- plot_grid(S3A.1,S3A.2, ncol = 2, rel_widths = c(0.5,1))

S3

#### Checking Topt and Isolation temperature ####

ave_params_gcplyr_combined_major_taxa$temperature_of_isolation <- factor(ave_params_gcplyr_combined_major_taxa$temperature_of_isolation)
ave_params_gcplyr_combined_major_taxa$climate_profile <- factor(ave_params_gcplyr_combined_major_taxa$`Climate profile`)

toi_topt_model <- (aov(data = ave_params_gcplyr_combined_major_taxa, topt ~ temperature_of_isolation))
summary(toi_topt_model)

TukeyHSD()
posthoc_toi_topt_model <- TukeyHSD(toi_topt_model)
posthoc_toi_topt_model

toi_topt_cp_model <- (aov(data = ave_params_gcplyr_combined_major_taxa, topt ~ temperature_of_isolation + climate_profile))
summary(toi_topt_cp_model)
posthoc_toi_topt_cp_model <- TukeyHSD(toi_topt_cp_model)
posthoc_toi_topt_cp_model


toi_topt_cp_model_hv <- (aov(data = ave_params_gcplyr_combined_major_taxa %>%
                            dplyr::filter(climate_profile == "Hot variable"), 
                            topt ~ temperature_of_isolation))
summary(toi_topt_cp_model_hv)
posthoc_toi_topt_cp_model_hv <- TukeyHSD(toi_topt_cp_model_hv)
posthoc_toi_topt_cp_model_hv

toi_topt_cp_model_hs <- (aov(data = ave_params_gcplyr_combined_major_taxa %>%
                               dplyr::filter(climate_profile == "Hot stable"), 
                             topt ~ temperature_of_isolation))
summary(toi_topt_cp_model_hs)
posthoc_toi_topt_cp_model_hs <- TukeyHSD(toi_topt_cp_model_hs)
posthoc_toi_topt_cp_model_hs


toi_topt_cp_model_cv <- (aov(data = ave_params_gcplyr_combined_major_taxa %>%
                               dplyr::filter(climate_profile == "Cold variable"), 
                             topt ~ temperature_of_isolation))
summary(toi_topt_cp_model_cv)
posthoc_toi_topt_cp_model_cv <- TukeyHSD(toi_topt_cp_model_cv)
posthoc_toi_topt_cp_model_cv

toi_topt_cp_model_cs <- (aov(data = ave_params_gcplyr_combined_major_taxa %>%
                               dplyr::filter(climate_profile == "Cold stable"), 
                             topt ~ temperature_of_isolation))
summary(toi_topt_cp_model_cs)
posthoc_toi_topt_cp_model_cs <- TukeyHSD(toi_topt_cp_model_cs)
posthoc_toi_topt_cp_model_cs



S8A <- ave_params_gcplyr_combined_major_taxa %>%
  filter(metric == "auc")%>%
  ggplot(aes(x=temperature_of_isolation, y= topt))+
  #scale_color_viridis_b()+
  geom_jitter(alpha = 0.8,size= 5,aes(color = Order))+
  geom_boxplot(alpha = 0.8, aes(group = temperature_of_isolation),outliers = FALSE)+
  #facet_wrap(~Order,ncol = 1)+
  #ylim(16,40)+
  #ylim(0.5,3.5)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression(Topt[rate]~" - MAT"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "turbo")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)
  )

S8B <- ave_params_gcplyr_combined_major_taxa %>%
  filter(metric == "auc")%>%
  ggplot(aes(x=temperature_of_isolation, y= topt))+
  #scale_color_viridis_b()+
  geom_jitter(alpha = 0.8,size= 5,aes(color = Order))+
  geom_boxplot(alpha = 0.8, aes(group = temperature_of_isolation),outliers = FALSE)+
  facet_wrap(~`Climate profile`,ncol = 4)+
  #ylim(16,40)+
  #ylim(0.5,3.5)+
  #scale_fill_manual(values = order_color)+
  #ylab(expression(Topt[rate]~" - MAT"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  #scale_color_viridis_d(option = "turbo")+
  scale_color_manual(values = c("#DE5925","#0090B5"))+
  scale_size_continuous(guide = "none")+
  theme(strip.background = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)
  )

S8 <- plot_grid(S8A,S8B, align = "h", 
                 axis = "l", 
                 ncol = 2, 
                 labels = c("A","B"),
                 label_fontfamily = "Atkinson Hyperlegible Next VF Light")


S8

#### Checking if the global patterns in the data stem from the overrepresentation of Pseudomonadales ####

## Splitting the ave_params_gcplyr_combined_major_taxa_auc dataframe by Order
ave_params_gcplyr_combined_major_taxa_auc_Pseudomonadales <- ave_params_gcplyr_combined_major_taxa_auc %>% filter(Order == "Pseudomonadales")
ave_params_gcplyr_combined_major_taxa_auc_Enterobacterales <- ave_params_gcplyr_combined_major_taxa_auc %>% filter(Order == "Enterobacterales")

## Making an empty list to hold randomized dataframes from the Pseudomonadales df


# number of dataframes we want to generate
n_strata <- 50

set.seed(123)  # For reproducibility

stratified_Pseudomonadales <- replicate(
  n_strata,
  ave_params_gcplyr_combined_major_taxa_auc_Pseudomonadales %>%
    group_by(Site, temperature_of_isolation) %>% # Stratify the dataframe by Site and Isolation temperature
    slice_sample(n = 1) %>%   # Randomly sample one element per group
    ungroup(),
  simplify = FALSE
)

stratified_datasets <- lapply(
  stratified_Pseudomonadales,
  function(df) {
    bind_rows(df, ave_params_gcplyr_combined_major_taxa_auc_Enterobacterales)
  }
)
# Make function to calculate phylosig 

n_permutations_2 <- 1000

phylogenetic_signal_mixed_model_f <- function(trait_name, data) {
  
  # Fit hierarchical linear model with Satterthwaite df
  formula <- as.formula(paste(trait_name, "~ PC1 + PC2 + temperature_of_isolation + (1|Order/Genus) + (1|Site)"))
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
  permuted_proportions <- matrix(NA, nrow = n_permutations_2, ncol = nrow(observed_var))
  colnames(permuted_proportions) <- observed_var$grp
  
  for(i in 1:n_permutations_2){
    perm_data <- data
    perm_data[[trait_name]] <- sample(perm_data[[trait_name]])
    
    perm_model <- lmer(formula, data = perm_data)         
    perm_var <- as.data.frame(VarCorr(perm_model))   
    permuted_proportions[i, ] <- perm_var$vcov / sum(perm_var$vcov)
  }
  
  ci_lower <- apply(permuted_proportions, 2, quantile, probs = 0.025)
  ci_upper <- apply(permuted_proportions, 2, quantile, probs = 0.975)
  
  confint_df <- data.frame(
    observed = observed_var$proportion,
    lower_ci = ci_lower,
    upper_ci = ci_upper,
    Factor = observed_var$grp,
    Trait = trait_name,
    significant = observed_var$proportion < ci_lower | observed_var$proportion > ci_upper
  )
  
  confint_df$Factor <- recode(
    confint_df$Factor,
    "Order:Genus" = "Genus",
    .default = confint_df$Factor
  )
  
  list(fixed_var = fixed_var, phylo_signal = confint_df)
}

# Run mixed linear regression models

phylogenetic_signal_pipeline <- function(df) {
  
  # Apply to all traits
  results_aucmax <- phylogenetic_signal_mixed_model_f("rmax", df)
  results_topt   <- phylogenetic_signal_mixed_model_f("topt", df)
  results_ctmax  <- phylogenetic_signal_mixed_model_f("ctmax", df)
  results_ctmin  <- phylogenetic_signal_mixed_model_f("ctmin", df)
  
  # Combine phylogenetic signal data
  phylogenetic_signal_thermal_traits_auc <- bind_rows(
    results_aucmax$phylo_signal,
    results_topt$phylo_signal,
    results_ctmax$phylo_signal,
    results_ctmin$phylo_signal
  )
  
  # Combine fixed effect contributions
  fixed_effects_all_auc <- bind_rows(
    cbind(Trait="rmax", results_aucmax$fixed_var),
    cbind(Trait="topt", results_topt$fixed_var),
    cbind(Trait="ctmax", results_ctmax$fixed_var),
    cbind(Trait="ctmin", results_ctmin$fixed_var)
  )
  
  # Clean random effects
  random_effects_clean_auc <- phylogenetic_signal_thermal_traits_auc %>%
    mutate(
      metric = "auc",
      Variable = `Factor`,
      Effect = "random"
    ) %>%
    rename(
      Variance_observed = observed,
      Lower_CI = lower_ci,
      Upper_CI = upper_ci,
      Significant = significant
    ) %>%
    dplyr::select(Trait, metric, Variable, Effect, Variance_observed, Lower_CI, Upper_CI, Significant)
  
  # Clean fixed effects
  fixed_effects_clean_auc <- fixed_effects_all_auc %>%
    mutate(
      metric = "auc",
      Variable = grp,
      Effect = "fixed",
      p_value = `Pr(>|t|)`
    ) %>%
    rename(
      Estimate = Estimate,
      Std_Error = `Std. Error`,
      t_value = `t value`,
      Partial_R2 = share
    ) %>%
    dplyr::select(Trait, metric, Variable, Effect, Estimate, Std_Error, t_value, p_value, Partial_R2)
  
  # Combine
  combined_effects_lmer_auc <- bind_rows(
    fixed_effects_clean_auc,
    random_effects_clean_auc
  ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  return(combined_effects_lmer_auc)
}

set.seed(123)

combined_results_phylosig <- lapply(
  seq_along(stratified_datasets),
  function(i) {
    phylogenetic_signal_pipeline(stratified_datasets[[i]]) %>%
      mutate(iteration = i)
  }
)

combined_results_phylosig_df <- bind_rows(combined_results_phylosig)


combined_results_phylosig_df <- combined_results_phylosig_df %>%
  mutate(Variable = dplyr::recode(Variable,
                           "temperature_of_isolation" = "TOI"))


combined_results_phylosig_df %>%
  filter(Effect == "random")%>%
  group_by(Trait,Variable)%>%
  dplyr::summarize(Variance_observed = mean(Variance_observed),
            Lower_CI = mean(Lower_CI),
            Upper_CI = mean(Upper_CI))%>%
  ggplot(aes(x = Variance_observed, y = Variable))+
  geom_point(data = combined_results_phylosig_df  %>% 
               filter(Effect == "random"), 
             aes(x = Variance_observed, y = Variable), 
             alpha = 0.3,
             size = 3)+
  geom_point(size = 6)+
  geom_errorbar(aes(x = Variance_observed, y = Variable, xmin = Lower_CI, xmax = Upper_CI))+
  facet_grid(~Trait)+
  #scale_fill_manual(values = c("#DE5925","#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )

combined_results_phylosig_df <- combined_results_phylosig_df %>%
  mutate(Variable = dplyr::recode(Variable,
                           "Genus:Order" = "Genus"))

S5B <- combined_results_phylosig_df %>%
  filter(Effect == "random")%>%
  group_by(Trait,Variable)%>%
  dplyr::summarize(Variance_observed = mean(Variance_observed),
            Lower_CI = mean(Lower_CI),
            Upper_CI = mean(Upper_CI),
            Significant_Proportion = mean(Significant, na.rm = TRUE))%>%
  ggplot(aes(x = Variance_observed, y = Variable))+
  geom_point(data = combined_results_phylosig_df  %>% 
               filter(Effect == "random"), 
             aes(x = Variance_observed, y = Variable), 
             alpha = 0.3,
             size = 3)+
  geom_errorbar(aes(x = Variance_observed, y = Variable, xmin = Lower_CI, xmax = Upper_CI))+
  geom_point(size = 7, color = "white")+
  geom_point(size = 6, aes(color = Significant_Proportion))+
  facet_grid(~Trait)+
  scale_color_scico(palette = "vik")+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )


combined_effects_lmer_auc <- fread("PATH_TO_YOUR_DIRECTORY\\05) TCS_figures\\03_Output\\combined_effects_lmer_auc.csv")


combined_effects_lmer_auc$iteration <- 0

combined_effects_lmer_auc <- combined_effects_lmer_auc %>%
  mutate(Variable = dplyr::recode(Variable,
                         "temperature_of_isolation" = "TOI"))

comparison_stratified_full <- bind_rows(combined_results_phylosig_df,combined_effects_lmer_auc)

comparison_stratified_full %>%
  filter(Effect == "fixed")%>%
  filter(iteration > 0)%>%
  ggplot(aes(x = p_value))+
  geom_density(aes(y = after_stat(scaled)),
               fill = "grey80")+
  geom_vline(data = comparison_stratified_full %>%
               filter(Effect == "fixed")%>%
               filter(iteration == 0), aes(xintercept = p_value), 
                                           color = "red3",
             lwd = 1)+
  geom_vline(aes(xintercept = 0.05), 
             color = "black",
             lwd = 0.75,
             linetype = "dashed")+
  facet_grid(Variable~Trait)+
  #scale_fill_manual(values = c("#DE5925","#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )

S5A <- comparison_stratified_full %>%
  filter(Effect == "fixed")%>%
  filter(iteration > 0)%>%
  ggplot(aes(x = Partial_R2))+
  geom_density(aes(y = after_stat(scaled)),
               fill = "grey80")+
  geom_vline(data = comparison_stratified_full %>%
               filter(Effect == "fixed")%>%
               filter(iteration == 0), aes(xintercept = Partial_R2), 
             color = "red3",
             lwd = 1)+
  facet_grid(Variable~Trait)+
  #scale_fill_manual(values = c("#DE5925","#0090B5"))+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        #axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white")
  )


fixed_effects_p_value_distribution <- combined_results_phylosig_df %>%
  filter(Effect == "fixed")%>%
  dplyr::select(c(Partial_R2,p_value,Trait,Variable))

library(poolr)

significance_fixed_effects <- fixed_effects_p_value_distribution %>%
  group_by(Trait,Variable)%>%
  dplyr::summarize(
    "mean_Partial_R2" = mean(Partial_R2),
    "Significance" = poolr::fisher(p_value)[["p"]],
    "Significant" = case_when(Significance <= 0.05 ~ "Yes",
                              TRUE ~ "No")
  )


combined_effects_lmer_auc_fixed_effects <- combined_effects_lmer_auc %>%
  filter(Effect == "fixed")%>%
  dplyr::select(c(p_value,Trait,Variable))%>%
  mutate(
    "Significant_full_dataset" = case_when(p_value <= 0.05 ~ "Yes",
                              TRUE ~ "No")
  )

significance_fixed_effects_comparison <- left_join(significance_fixed_effects, 
                                                   combined_effects_lmer_auc_fixed_effects, 
                                                   by = c("Trait","Variable"))



S5 <- plot_grid(S5A, S5B, ncol = 1, align = "hv",axis = "x", labels = c("A","B"))
S5

###### Evaluating the phylogenetic signal of the multivariate association between thermal traits ######

library(factoextra)
ave_params_gcplyr_combined_major_taxa_auc_Enterobacterales <- ave_params_gcplyr_combined_major_taxa_auc %>%
  filter(Order == "Enterobacterales")

ave_params_gcplyr_combined_major_taxa_auc_Pseudomonadales <- ave_params_gcplyr_combined_major_taxa_auc %>%
  filter(Order == "Pseudomonadales")

pc_Enterobacterales <- prcomp(ave_params_gcplyr_combined_major_taxa_auc_Enterobacterales %>%
                                dplyr::select(c(rmax,topt,ctmin,ctmax)),
                              center = TRUE,
                              scale. = TRUE)

pc_Pseudomonadales <- prcomp(ave_params_gcplyr_combined_major_taxa_auc_Pseudomonadales %>%
                               dplyr::select(c(rmax,topt,ctmin,ctmax)),
                             center = TRUE,
                             scale. = TRUE)

pc_major_taxa_auc <- prcomp(ave_params_gcplyr_combined_major_taxa_auc %>%
                              dplyr::select(c(rmax,topt,ctmin,ctmax)),
                            center = TRUE,
                            scale. = TRUE)

S15A <- fviz_pca_biplot(
  label = "var",
  pc_Enterobacterales,
  habillage   = ave_params_gcplyr_combined_major_taxa_auc_Enterobacterales$Genus %>%
    droplevels(),
  addEllipses = TRUE,
  repel       = TRUE,
  palette = full_sites_palette
)+
  ylab("PC2 (29.7%)")+
  xlab("PC1 (43.1%)")+
  theme_classic(base_size = 14,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(
    plot.title = element_blank(),
    legend.position = "bottom",
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

S15B <- fviz_pca_biplot(
  label = "var",
  pc_Pseudomonadales,
  habillage   = ave_params_gcplyr_combined_major_taxa_auc_Pseudomonadales$Genus %>%
    droplevels(),
  addEllipses = TRUE,
  repel       = TRUE,
  palette = full_sites_palette
)+
  ylab("PC2 (33.6%)")+
  xlab("PC1 (35.9%)")+
  theme_classic(base_size = 14,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(
    plot.title = element_blank(),
    legend.position = "bottom",
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

F2D <- fviz_pca_biplot(
  title = NULL,
  label = "var",
  pc_major_taxa_auc,
  habillage   = ave_params_gcplyr_combined_major_taxa_auc$Order %>%
    droplevels(),
  addEllipses = TRUE,
  repel       = TRUE,
  palette = c("#DE5925","#0090B5")
)+
  ylab("PC2 (33.2%)")+
  xlab("PC1 (40.2%)")+
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  theme(
    legend.position = "top",
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


pc_major_taxa_auc_df <- data.frame(
  "PC1" = pc_major_taxa_auc$x[,1],
  "PC2" = pc_major_taxa_auc$x[,2],
  "seqID" = ave_params_gcplyr_combined_major_taxa_auc$seqID,
  "Order" = ave_params_gcplyr_combined_major_taxa_auc$Order)



#
#
# Here we apply the methods detailed in Lennon et al (2016) to calculate the phylogenetic
# signal of bacterial traits to our thermal traits, by fitting hierarchical linear models
# for each trait and metric.
#
#

pc_major_taxa_auc_df$pc1_thermal_traits <- pc_major_taxa_auc_df$PC1

ave_params_gcplyr_combined_major_taxa_auc_pctt <- ave_params_gcplyr_combined_major_taxa_auc %>% 
  left_join(pc_major_taxa_auc_df %>% dplyr::select(seqID,pc1_thermal_traits),by = "seqID")

n_permutations <- 10000

phylogenetic_signal_mixed_model_multivar <- function(trait_name, data) {
  
  # Fit hierarchical linear model with Satterthwaite df
  formula <- as.formula(
    paste(
      trait_name,
      "~ PC1 + PC2 + temperature_of_isolation + (1|Order/Genus) + (1|Site)"
    )
  )
  
  model <- lmer(formula, data = data)
  
  # Global model to extract Marginal and Conditional R2 values
  
  r2_values <- MuMIn::r.squaredGLMM(model)
  
  model_fit <- data.frame(
    Trait = trait_name,
    metric = "auc",
    Variable = c("Marginal_R2", "Conditional_R2"),
    Effect = "model_fit",
    Variance_observed = c(
      r2_values[1, "R2m"],
      r2_values[1, "R2c"]
    )
  )
  
  ## Partitioning the Marginal P-values and R2
  
  part_model <- partR2(
    model,
    partvars = c("PC1", "PC2", "temperature_of_isolation"),
    data = data
  )
  
  ## Extract fixed effects table
  fixed_summary <- as.data.frame(summary(model)$coefficients)
  fixed_summary$grp <- rownames(fixed_summary)
  rownames(fixed_summary) <- NULL
  
  ## Keep only fixed effects included in partR2 output
  fixed_var <- fixed_summary[
    fixed_summary$grp %in% part_model$R2$term,
  ]
  
  ## Match partial R² values
  fixed_var$share <- part_model$R2$estimate[
    match(fixed_var$grp, part_model$R2$term)
  ]
  
  ## Partitioning the random effects
  observed_var <- as.data.frame(VarCorr(model))
  
  observed_total <- sum(observed_var$vcov)
  
  observed_var$proportion <- observed_var$vcov / observed_total
  
  # Running a permutation test to calculate the robustness of our phylogenetic signal
  permuted_proportions <- matrix(
    NA,
    nrow = n_permutations,
    ncol = nrow(observed_var)
  )
  
  colnames(permuted_proportions) <- observed_var$grp
  
  for (i in seq_len(n_permutations)) {
    
    perm_data <- data
    perm_data[[trait_name]] <- sample(perm_data[[trait_name]])
    
    perm_model <- lmer(formula, data = perm_data)
    
    perm_var <- as.data.frame(VarCorr(perm_model))
    
    permuted_proportions[i, ] <-
      perm_var$vcov / sum(perm_var$vcov)
  }
  
  ci_lower <- apply(
    permuted_proportions,
    2,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  )
  
  ci_upper <- apply(
    permuted_proportions,
    2,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  )
  
  confint_df <- data.frame(
    observed = observed_var$proportion,
    lower_ci = ci_lower,
    upper_ci = ci_upper,
    Factor = observed_var$grp,
    Trait = trait_name,
    significant =
      observed_var$proportion < ci_lower |
      observed_var$proportion > ci_upper
  )
  
  confint_df$Factor <- dplyr::recode(
    confint_df$Factor,
    "Order:Genus" = "Genus",
    .default = confint_df$Factor
  )
  
  list(
    fixed_var = fixed_var,
    phylo_signal = confint_df,
    model_fit = model_fit
  )
}

phylogenetic_signal_pipeline_multivar <- function(df) {
  
  ## Run model for each trait
  results_pc1_thermal_traits <- phylogenetic_signal_mixed_model_multivar("pc1_thermal_traits", df)
  
  ## Combine all data into a single dataframe
  phylogenetic_signal_pc1_thermal_traits <- results_pc1_thermal_traits$phylo_signal
  
  ## Combine fixed effects into a single dataframe
  fixed_effects_pc1_thermal_traits <- results_pc1_thermal_traits$fixed_var
  
  ## Combine global R2 values
  model_fit_pc1_thermal_traits <- results_pc1_thermal_traits$model_fit
  
  ## Combine random effects into a single dataframe
  random_effects_pc1_thermal_traits_clean<- phylogenetic_signal_pc1_thermal_traits %>%
    mutate(
      metric = "auc",
      Variable = Factor,
      Effect = "random"
    ) %>%
    rename(
      Variance_observed = observed,
      Lower_CI = lower_ci,
      Upper_CI = upper_ci,
      Significant = significant
    ) %>%
    dplyr::select(
      metric,
      Variable,
      Effect,
      Variance_observed,
      Lower_CI,
      Upper_CI,
      Significant
    )
  
  ## Combine fixed effects into a single df too
  fixed_effects_pc1_thermal_traits_clean <- fixed_effects_pc1_thermal_traits %>%
    mutate(
      metric = "auc",
      Variable = grp,
      Effect = "fixed",
      p_value = `Pr(>|t|)`
    ) %>%
    rename(
      Estimate = Estimate,
      Std_Error = `Std. Error`,
      t_value = `t value`,
      Partial_R2 = share
    ) %>%
    dplyr::select(
      metric,
      Variable,
      Effect,
      Estimate,
      Std_Error,
      t_value,
      p_value,
      Partial_R2
    )
  
  ## Select relevant columns from the global fit df
  model_fit_pc1_thermal_traits_clean <- model_fit_pc1_thermal_traits %>%
    mutate(
      Lower_CI = NA_real_,
      Upper_CI = NA_real_,
      Significant = NA
    ) %>%
    dplyr::select(
      metric,
      Variable,
      Effect,
      Variance_observed,
      Lower_CI,
      Upper_CI,
      Significant
    )
  
  ## Put everything together into Supplementary Table 2
  combined_effects_pc1_thermal_traits <- bind_rows(
    fixed_effects_pc1_thermal_traits_clean,
    random_effects_pc1_thermal_traits_clean,
    model_fit_pc1_thermal_traits_clean
  ) %>%
    mutate(
      across(where(is.numeric), ~ round(.x, 3))
    )
  
  return(
    list(
      combined_effects = combined_effects_pc1_thermal_traits,
      phylo_signal = phylogenetic_signal_pc1_thermal_traits,
      fixed_effects_all_auc = fixed_effects_pc1_thermal_traits,
      model_fit = model_fit_pc1_thermal_traits
    )
  )
}

# Running pipeline 
combined_effects_pc1_thermal_traits <-
  phylogenetic_signal_pipeline_multivar(
    ave_params_gcplyr_combined_major_taxa_auc_pctt
  )


phylogenetic_signal_multivar <- combined_effects_pc1_thermal_traits$phylo_signal

phylogenetic_signal_multivar <- phylogenetic_signal_multivar %>%
  dplyr::mutate(Factor = dplyr::recode(Factor,
                                       "Genus:Order" = "Genus"))

phylogenetic_signal_multivar$Factor <- factor(phylogenetic_signal_multivar$Factor, 
                                              levels = c("Genus","Order","Residual","Site"))


F2E <- phylogenetic_signal_multivar %>% 
  #filter(!`Taxonomic Level` == "Family") %>% 
  #mutate(Trait = fct_relevel(Trait,"ctmin_rate","ctmax_rate","topt_rate","rmax_rate")) %>%
  #mutate(Trait = fct_relevel(Trait,"ctmin_rate","ctmax_rate","topt_rate","ratemax" )) %>%
  ggplot(aes(x = Factor, y = observed)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), linewidth = 0.5,width = 0.2, color = "gray50") +
  geom_point(aes(color = significant), size = 6) +
  theme_classic(base_size = 18,base_family = "Atkinson Hyperlegible Next VF Light")+
  ylim(0,1)+
  scale_y_continuous(breaks = seq(0, 1, 0.5))+
  coord_flip()+
  #facet_wrap( ~Factor,ncol = 4)+
  labs(
    y = "Proportion of Variance Explained",
    color = "Significant"
  )  +
  scale_color_manual(values = c("grey60","black"))+
  theme(strip.background = element_blank(),
        axis.line.y = element_blank(),
        #     axis.text.y = element_blank(),
        #      axis.text.x = element_blank(),
        #     axis.ticks.y = element_blank(),
        #    axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        # strip.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "top",
        legend.direction = "horizontal"
  )



F2DE<- plot_grid(F2D, F2E, 
                 ncol = 2, 
                 align = "hv",axis = "x", 
                 labels = c("A","B"),
                 label_fontfamily = "Atkinson Hyperlegible Next VF Light")

S15<- plot_grid(S15A, S15B, 
                ncol = 2, 
                align = "hv",axis = "xy", 
                labels = c("C","D"),
                label_fontfamily = "Atkinson Hyperlegible Next VF Light")

