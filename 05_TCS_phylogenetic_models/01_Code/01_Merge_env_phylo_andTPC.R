######################################################
# Paper: Conserved upper thermal limits and small   #
# safety margins in soil copiotrophic bacteria.     #

#######################################################
########## Phylogenetic correlations  #################
#######################################################
## Code: Ariel Favier (afavier@uci.edu)
## Edited by; Maria Rebolleda-Gomez (mreboll1@uci.edu)
## Last updated: Sept 3 2026
#######################################################

### Load libraries
library(data.table)
library(ape)
library(phytools)
library(tidyverse)

### Set parent directory 
setwd("PATH TO DIRECTORY/Paper_PhylogeneticThermalConstraints")

### Load and format TPC data
# TPC parameter output file of the 250324_TCS_TPC.R script
ave_params_gcplyr_combined <- fread("04_TCS_thermal_performance_curves/02_Output/251029_tcs_params_gcplyr.csv.gz")

# TPC fits output file of the 250324_TCS_TPC.R script (we replace negative values with 0)
ave_preds_gcplyr_combined <- fread("04_TCS_thermal_performance_curves/02_Output/251029_tcs_fits_gcplyr.csv.gz")


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
taxonomy <-  fread("02_TCS_genome_collection/03_Output/Kbase_taxonomy_quality_filtered.csv")

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


write.files <- function(data_sets){
    for (i in 1:length(data_sets)){
        data <- get(data_sets[i])
        file <- paste("05_TCS_phylogenetic_models/03_Output/", data_sets[i], ".csv.gz", sep = "")
        fwrite(data, file)
    }
}

write.files(c("ave_params_gcplyr_combined_common_taxa_auc",
    "ave_params_gcplyr_combined_common_taxa_rate",
    "ave_preds_gcplyr_combined_common_taxa_auc",
    "ave_preds_gcplyr_combined_common_taxa_rate"))

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

write.files(c("ave_params_gcplyr_combined_major_taxa_auc",
    "ave_params_gcplyr_combined_major_taxa_rate",
    "ave_preds_gcplyr_combined_major_taxa_auc",
    "ave_preds_gcplyr_combined_major_taxa_rate"))


### Open and format tree and taxonomy ####
tree_kbase <-  read.tree(file="05_TCS_phylogenetic_models/02_Data/TCS_genomes_annotated_tree_without_others-labels.newick")

# Extract tip labels from the tree
tree_kbase$tip.label <- gsub("output.*DRAM", "", tree_kbase$tip.label)

# Now we organize the tree, by selecting the taxa we want to keep in the tree as those that have trait data
tipstokeep <- intersect(tree_kbase$tip.label,unique(as.character(ave_params_gcplyr_combined_common_taxa_auc$seqID)))
tpc_tree <- keep.tip(tree_kbase,tipstokeep) 

# We save the trimmed tree
write.tree(tpc_tree,file="05_TCS_phylogenetic_models/03_Output/tpc_tree_kbase_relaxed.tre")

# We then subset the trait metadata to map on our phylogenetic tree tips by metric 
tpc_tree_metadata_auc <- ave_params_gcplyr_combined_common_taxa_auc %>% filter(as.character(seqID) %in% tpc_tree$tip.label & metric == "auc") %>% droplevels()
tpc_tree_metadata_rate <- ave_params_gcplyr_combined_common_taxa_rate %>% filter(as.character(seqID) %in% tpc_tree$tip.label & metric == "rate") %>% droplevels()

tpc_tree_metadata_auc$iso1 <- tpc_tree_metadata_auc$seqID
tpc_tree_metadata_auc$iso2 <- tpc_tree_metadata_auc$seqID

# We save the tree metadata
write.files(c("tpc_tree_metadata_auc",
    "tpc_tree_metadata_rate"))
