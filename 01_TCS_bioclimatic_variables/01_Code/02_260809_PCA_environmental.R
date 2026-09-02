######################################################
# Paper: Conserved upper thermal limits and small   #
# safety margins in soil copiotrophic bacteria.     #

#######################################################
####### Calculating PCA of bioclimatic variables ######
#######################################################
## Code: Ariel Favier (afavier@uci.edu)
## Last updated: Aug 19 2026
#######################################################

# Load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(viridis)
library(viridisLite)
linrar(PCAtools)
font_import("Futura")
loadfonts("win")
fonts()

# Setting the directory 
# parent <- ("YOUR PATH/Paper_PhylogeneticThermalConstraints/") set the path to the parent directory
paste(parent,"01_TCS_bioclimatic_variables/", sep="") %>% setwd

environment_database <- fread("03_Output/environment_database_tcs_250512.csv.gz") 



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

write.csv(PCs_env_full, "03_Outpu/PCs_env_full.csv", row.names = FALSE)


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

write.csv(PCs_env_noprec, "03_Output/PCs_env_noprec.csv", row.names = FALSE)

View(pca_result_env_variables_noprec$loadings)



pc_scores_noprec <- as.data.frame(pca_result_env_variables_noprec$rotated)
pc_scores_noprec$Site_initials <- rownames(pc_scores_noprec)



environment_database_noprec <- environment_database %>% left_join(pc_scores_noprec, by = "Site_initials")



write.csv(environment_database_noprec,"03_Output/environment_database_tcs_noprec_250512.csv",row.names = FALSE)

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

