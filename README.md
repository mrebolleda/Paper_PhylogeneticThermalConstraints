# Repository for paper "Conserved upper thermal limits and small safety margins in soil copiotrophic bacteria"

### Authors
Corresponding: María Rebolleda-Gómez (mreboll1@uci.edu)  
Ariel Favier (afavier@uci.edu) (data and code curation)
Alejandra Hernández-Terán (ipomeap@gmail.com) 
Celia Symons (csymons@uci.edu) 

### Abstract
One of the key uncertainties in climate change models is how microbes will adapt to rising temperatures. Large-scale comparisons of bacterial thermal performances show a clear boundary between mesophiles and thermophiles. Here, we investigated whether phylogenetic constraints limit the adaptive potential of bacteria to warming soils. Focusing on copiotrophs within Gammaproteobacteria, we found that both thermal optima and upper thermal limits are constrained; variation in these traits decreases above 42°C across the phylogeny, with minimal influence of present-day bioclimatic variables. This, along with the reduced thermal safety margins found in fluctuating hot climates, suggests that many isolates may already be maladapted to local temperature variability. Our findings indicate that these constraints interact with the geometry of thermal performance curves, imposing a trade-off between high-temperature performance and the risk of substantial fitness losses above Topt. Overall, this work underscores potential limits to thermal adaptation and their implications for bacterial fitness.

## Overview
The overall aim of this work was to evaluate the role of phylogenetic constraints in shaping the thermal performance curves (TPCs) of bacteria and their adaptability to warming soils. For that purpose, we obtained bacterial isolates from soil samples spanning a wide range of climates across California and measured their TPCs in minimal media. We then contrasted their shape against hypotheses of local thermal adaptation and phylogenetic conservatism of thermal traits. 

This repository contains data files and R scripts used to analyze the data for the paper. Each folder contains the necessary data for the script (s) to run independently and their outputs. Most data files already contain the necessary subsets of metadata for plotting and wrangling purposes.

- **01_TCS_bioclimatic_variables** - Contains two R scripts used to analyze climate time series. 
- **02_TCS_Genome_collection** - Contains a short R script to quality filter the assembled genomes from our KBase narrative, along with the taxonomy assignment of the genome collection.
- **03_TCS_Growth_curves** - Contains the aggregated raw optical density data for our growth kinetics, the R script used to process and curate the data, intermediate steps of data filtering, the final output of the script (used as input for the TPC fitting script in the next folder), and a PDF file with plotted individual growth curves.
- **04_TCS_Thermal_performance_curves** - Contains the main script to fit the TPCs from the growth curve parameter file, along with several datasets summarized in different ways for plotting or analytic purposes, and PDF files with plotted TPCs for different growth metrics.
- **05_TCS_Figures** - Contains a very large R script and input data used to make the main and supplementary figures, along with all statistical analyses. 

## General instructions on how to use this repository
Within each major directories there are three directories: 
- **01_Code** - Contains the scripts necessary for the analyses. 
- **02_Data** - Contains raw data or aggregated files with the raw data from (originally) disparate files. 
- **03_Output** - Contains the files produced by the scripts in `01_Code`. 
*Some major directories do not have a data folder and instead use output from previous steps (numbers). We provide the ouptup files in this repository and therefore be able to skip steps.*

At the beginning of each script there is a code to define the parent directory (where is this directory stored in your computer) and set the working directory. After that, all code to open files is set with a relative path based on the structure of this repository. 

## Supplementary data not contained in this repository
- Raw sequencing reads are available in the NCBI Sequence Read Archive (SRA) under BioProject accession PRJNA1414602. 
- We partially processed the genomic data using the KBase Narrative (https://doi.org/10.25982/200311.171/3014786)

## R code and package versions used
|     Package     |     Version   |    In-text citation                                 |
|-----------------|---------------|-----------------------------------------------------|
|     brms        |     2.23.0    |    Main Manuscript                                  |
|     lmerTest    |     3.1.3     |    Main Manuscript                                  |
|     partR2      |     0.9.2     |    Main Manuscript                                  |
|     stats       |     4.3.2     |    Main Manuscript and Supplementary information    |
|     ggplot2     |     4.01      |    Main Manuscript                                  |
|     corrplot    |     0.92      |    Main Manuscript                                  |
|     ggtree      |     3.10.1    |    Main Manuscript                                  |
|     gcplyr      |     1.10.0    |    Main Manuscript and Supplementary information    |
|     rTPC        |     1.0.5     |    Main Manuscript and Supplementary information    |

## Detailed repository layout
### 01_TCS_bioclimatic_variables
- **01_Code:** Has all the R code necessary for characterizing sites and climate profiles based on soil bioclimatic variables.
    - `01_251029_Environmental_data_processing.R`: This script allows us to compile and summarize the climate data from soil time series and precipitation data across sampling sites. A big portion of the script is dedicated to harmonizing measurements from different sites, determining a shared time interval among all sites, and calculating the bioclimatic variables across timescales. We gathered environmental data from weather stations at all the sampled sites. For all sites except the Merced Vernal Pools and Grassland Reserve (MVPGR) and Pinyon Flats/Pinyon Crest (PF), we obtained climate data from the Dendra database (https://dendra.science/). We obtained We obtained MVPGR data from the San Joaquin River weather station (https://cdec.water.ca.gov/) and we obtained PF data from the Boyd Deep Canyon weather station at Pinyon Crest (https://deepcanyon.ucnrs.org/weather-data/) *All the data was downloaded July 22nd 2024 and it is stored in the /01_TCS_bioclimatic_variables/02_Data/Raw_stable_files directory. Details of how the data was downloaded are described in detail in the paper's supplemental materials.*
    - `02_260809_PCA_environmental.R`: Using the environmental data formated in the previous script, we compress the multiple variables into two main Principal Components, that explain most of the differences between sites, and use them to generate broader climate categories based on temperature mean and variability. 
    - `03_260809_densityprobabilities_temp.R`: 
    - `251001_precipitation_db_comparison.R`: This script assesses the interchangeability of the precipitation data obtained from the UCNRS Dendra database and PRISM. Its main output is the supplemental figure S1
      
- **02_Data:** Contains relevant data (mostly metadata) for the analyses described above.
  - `UCM_temperature_data.csv.gz`:
  *We obtained these data through the CDEC portal (https://cdec.water.ca.gov/), selecting the tab “Historical Data”, with filters: Station ID = UCM; Sensor Number = 194-(hourly) - SOIL TEMP, DEPTH 1, Start Date = 2017-05-16, End Date = 2018-07-10; then downloaded as a .csv file*
    - Soil 100mm Temperature degC: Hourly recording of temperature values in degrees Celsius at a soil depth of 100 mm
    - STATION_ID: Weather station identifier. In this case, the station “UCM” is the one closest to the UC Merced campus, where our sampling site is located.
  - temperature-data-ucnrs-240502.csv.gz: 
  *We obtained these data manually on the Dendra website on July 22 2024, by accessing the UCNRS time series. This website has been updated as of August 8 2026, and the original datastreams can be accessed by clicking on the “Explore” and then “UC Nature Research” tabs.*
      - time: timestamp of temperature collection. Soil temperature values are measured every 10 minutes.
      - [Site name] Soil Temp 50 mm Avg degC: Average temperature values in degrees Celsius at a soil depth of 50 mm.
  - `TCS-reserves-MOD13A1-061-results.csv.gz` and `TCS-reserves-2-MOD13A1-061-results.csv.gz`: *Satellite data for plant growth season calculation by site based on NDVI. Downloaded from the NASA AppEARS (https://appeears.earthdatacloud.nasa.gov) website. We manually obtained these data by selecting the tab “Extract”, “Point”, “Start a new request”, “Upload coordinates from a file”, using the files “TCS_reserves.csv” and “TCS_reserves_2.csv”, within this folder. We obtained the 500 m point sample NDVI measurements at a 16-day resolution for all our sampling sites for the period 2013-2022.*
      - ID: Sampling site initials, equivalent to Site_initials in other datasets
      - MOD13A1_061__500m_16_days_NDVI: Normalized Difference Vegetation Index (NDVI) for a pixel size of 500 m/side over 16 days using MODIS VI products.
      - MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description: Quality assessment of the NDVI value obtained. We use this column to filter out unreliable values due to factors such as cloud cover.
  - `rainfall_data_ucnrs_240722.csv.gz`: *Precipitation data obtained from the Dendra database for sampling sites except the Merced Vernal Pools and Grassland Reserve and Pinyon Crest/Flats site. We obtained these data manually on the Dendra website on July 22 2024, by accessing the UCNRS time series. This website has been updated as of August 8 2026, and the original datastreams can be accessed by clicking on the “Explore” and then “UC Nature Research” tabs.*
      - time: timestamp of precipitation collection. Precipitation values are measured every 10 minutes.
      - [Site name]-mm: Yearly cumulative precipitation values in mm per site. Measurements start on 10/01 and finish on 09/30 of the next calendar year.
  - `rainfall_data_PRISM_240725.csv.gz`:
  *Precipitation data for all sites obtained from the PRISM Climate Group Time Series database for the period 2012-2021. File condensing all data downloaded from PRISM*
      - [Site name]: Monthly cumulative precipitation for a specific sampling site in mm
  - `PRISM_raw_files.tar.gz`: *Directory containing the raw files from PRISM. Raw files from PRISM (Site_initials corresponds to the variable of the same name across our dataset). The data obtained from each file were condensed into the rainfall_data_PRISM_240725.csv. We obtained each file manually from the PRISM website by selecting the nearest coordinate pair to each of our sampling sites, with a resolution of 4 km from January 2012 through December 2021. Each file contains a header with additional information on the download date and search values.*
  - `Pinyon_crest_temperature_data.csv.gz`
      - Soil Temp 50 mm Avg degC: Average temperature values in degrees Celsius at a soil depth of 50 mm, after applying offset calculations with Air temperature.
  - `air-soil-offset-bdc.csv`: This dataset allows us to estimate soil temperatures at the PF site from the difference between air and soil temperature time series of the nearby DC (Boyd Deep Canyon) site.
    - Deep Canyon Air Temp Avg degC: Average air temperature values in degrees Celsius, at 2 m aboveground.
    - Deep Canyon Soil Temp 50 mm Avg degC: Average temperature values in degrees Celsius at a soil depth of 50 mm.
  

      
- **03_Output:** Mainly contains intermediate datasets generated after time-intensive steps, and final output files to be used in other scripts
  - `climate_profile_temperature_probability.csv`: This dataframe assigns a probability for each temperature within a range of temperatures from -5 to 62.5°C, at a 0.1°C resolution, for each Climate Profile
  - `climate_profile_temperature_probability_growing_season.csv`: This dataframe assigns a probability for each temperature within a range of temperatures from -5 to 62.5°C, at a 0.1°C resolution, for each Climate Profile during the plant growth season
  - `dendra_data_long.csv.gz`: Long-format soil temperature data obtained from the Dendra database for sampling sites except the Merced Vernal Pools and Grassland Reserve and Pinyon Crest/Flats site. 
      - Soil Temp 50 mm Avg degC: Average temperature values in degrees Celsius at a soil depth of 50 mm.
      - Site: sampling site identifier, equivalent to Site_initials.
  - `environment_database_tcs_noprec_250512.csv.gz`: Contains the condensed temperature bioclimatic variables and PC loadings per site
  - `environment_database_tcs_full_250512.csv.gz`: Contains the condensed temperature and precipitation bioclimatic variables, and PC loadings per site
  - `PCs_env_noprec.csv`: This dataframe contains the values of principal component (PC) loadings (PC1-PC10) obtained while excluding precipitation variables, per site
  - `PCs_env_full.csv`: This dataframe contains the values of principal component (PC) loadings (PC1-PC10) obtained while including precipitation variables, per site
  - `precipitation_database_tcs_240730.csv.gz`: This dataframe contains all the precipitation-related bioclimatic variables across sites
  - `site_temperature_probability_growing_season.csv`: This dataframe assigns a probability for each temperature within a range of temperatures from -5 to 62.5°C, at a 0.1°C resolution, for each sampling site during the plant growth season
  - `summarized_climate_data_growing_season.csv`: This dataframe contains the average soil temperature time series by site, year, and month, but keeping only dates within the estimated plant growth season of each site
  - `summarized_climate_data.csv`
  - `temperature_database_tcs_240505.csv`: This dataframe contains all the soil temperature bioclimatic variables across sites

 
### 02_TCS_genome_collection
- **01_Code:** Has a short R script to parse out contaminated genomes
  - `TCS_genome_quality_assessment.R`: This script filters out contaminated genomes using the CheckM output from the KBase pipeline, and also removes them from the quast output
- **02_Data:** Contains relevant data for the analysis described above. Data here is all output from the KBase pipeline (KBase Narrative https://doi.org/10.25982/200311.171/3014786)
    - `TCS_genome_quality.csv`
      - Bin name: isolate genome identifier. The number at the beginning of the string is equivalent to the value of “seqID” or “curve_id” in other datasets. The substring “TCS_genomes_[number]_DRAM” refers to the batch of DRAM-annotated genomes to which the genome belonged within the KBase pipeline. 
      - Marker Lineage: CheckM infers the position of a query genome within a reference phylogenetic tree. This “Marker Lineage” is then used to find the fraction of marker genes that should be present in the genome (and with it its “completeness”), and to estimate the level of “contamination” with marker sequences from other clades.
      - A gene is considered “marker” when it is present as a single-copy gene in >97% of genomes within a lineage. Columns 0-5 represent the number of copies/versions associated with each Marker Gene from the inferred Maker Lineage gene set found in the query genome. These are then used to calculate the columns “Completeness” and “Contamination”.
      - Completeness: The values of column “0” is contrasted against the “Marker sets” column to obtain the completeness percentage, as in (“Marker sets” – “0”)/(“Marker sets”). We use this column to only keep highly complete genomes (most of them, because we sequenced isolate populations).
      - Contamination: Genome contamination is estimated from the number of multicopy marker genes identified in each marker set (Parks et al. 2015). We kept genomes with contamination levels below 5% (“Low” and “no detectable contamination”). 
    - `Taxonomy_file_kbase.csv`: Taxonomy assignment output from GTDB-Tk in the KBase pipeline
      - seqID: Isolate genome identifier
      - Domain, Phylum, Class, Order, Family, Genus, Species: Taxonomic assignment
      - Warnings: this column provides explanations for isolates with incomplete taxonomic assignments (e.g.: lacking Species name)
    - `quast_stats_unfiltered.csv`
      - curve_id: Isolate genome identifier
      - Other columns: refer to the “Metrics description” in the Quast (v5.3.0) manual
        
- **03_Output:** Contains output files to be used in other scripts
    - `Kbase_taxonomy_quality_unfiltered.csv`
    - `Kbase_taxonomy_quality_filtered.csv`
    - `quast_stats_filtered.csv`

### 03_TCS_growth_curves
- **01_Code:** Has all the R code necessary for fitting growth curves from the raw optical density data
    - `260327_TCS_growthcurves.R`
- **02_Data:** Contains relevant data (mostly metadata) for the analyses described above.
    - `raw_curves_df.csv.gz`: Aggregated optical density file over time for each kinetics measurement
      - seqID: Isolate genome identifier
      - abs: optical density measurement (OD620)
      - t: time since start of kinetics measurement (hs)
      - tmp: experimental growth temperature
      - replicate: experimental replicate
      - Isolation.temperature: temperature at which the organism was experimentally isolated (equivalent to “temperature_of_isolation” in other datasets. 
    - `CCA_gc_full.csv.gz`: Aggregated optical density file over time for each kinetics measurement. Very high resolution data from the carbon concentration analysis test (~12 temperatures). Most columns are the same as in the previous file. 
      - Glucose_stock_X: concentration of the glucose stock (10 = 0.007 M-C; 100 = 0.07 M-C)
    - `RC_gc_full.csv.gz`: Aggregated optical density file over time for each kinetics measurement. Very high resolution data from respirator isolates (Pseudomonadales).
    - `respirator_collection_platemap.csv.gz`: Metadata of respirator isolates used for the thermal performance curves (RC isolates). 
      - seqID: Isolate genome identifier
      - Taxonomic identification columns: Domain, Phylum, Class, Order, Family, Genus, Species. *Taxonomy assigned in the KBase pipeline*
      -  Reserve: Two-letter identifier for the UC reserve where we took the sample (Merced, ME; Anza Borrego, AB; Bodega, B; Boyd Deep Canyon, DC; Elliot Chaparral, EC; James San Jacinto, SJ; Jepson Prairie,J; McLaughlin,M; Point Reyes,P; Quail Ridge,QR; SNARL,S; Valentine,V; White Mountain,W; Yosemite Mariposa Grove,Y; Pinyon crest,PF).
      - Isolation.temperature: temperature at which the organism was experimentally isolated.
      - ID
      - DNAcc: DNA concentration
      - Area: Area of sampling (Southern California, Bay Area, Sierras). 
      - Color: Color of colony in chromogenic agar. 
      - old_well: Well of isolate in the stock plate.
      - Plate: Label of the plate where to find the stock isolate. 
      - well: Well for growth curves. 
      - carbon: In which carbon source was the isolate grown to measure growth. 

- **03_Output:** Mainly contains intermediate datasets generated after time-intensive steps, and final output files to be used in other scripts
  - `y_gcplyr_with_derivatives.csv`: Intermediate file used as input to calculate growth curve parameters. This file contains multiple fits per curve, based on possible combinations of smoothing parameters. In this case, we only provide two combinations, but we screened several others, which extended the processing time considerably.
      - seqID: Isolate genome identifier
      - Domain, Phylum, Class, Order, Family, Genus, Species: Taxonomy assignment, inherited from the KBase pipeline.
      - movavg_[number]: width of the sliding window used to smooth out the original “abs” values. A value of “13” means that the smoothing occurred using abs values along 13 time points (~2.17 hs).
      - derivpercap_[number]_[number]: width of the sliding window used to calculate the maximum derivative along the absorbance, or smoothed absorbance, data. A value of derivpercap_9_13 means that the derivative was calculated using a sliding window of 9 values of movavg_13.
  - `260327_fits_gcplyr.csv`: Predicted growth curves after fitting using gcplyr. Columns are the same as in `y_gcplyr_with_derivatives.csv`
  - `260327_params_gcplyr.csv`: Growth curve parameters obtained from the predicted growth curves
      - seqID: Isolate genome identifier
      - Domain, Phylum, Class, Order, Family, Genus, Species: Taxonomy assignment, inherited from the KBase pipeline.
      - t: time since start of kinetics measurement (hs)
      - tmp: experimental growth temperature
      - replicate: experimental replicate
      - temperature_of_isolation: temperature at which the organism was experimentally isolated (equivalent to “Isolation.temperature” in other datasets).
      - Site_initials: Sampling site from which the isolate was obtained
      - lag_time: lag time calculated by gcplyr. This measurement is extremely sensitive to noise and initial conditions, reason why we did not use it for our downstream analyses
      - r: maximum value of the chosen “derivpercap_[number]_[number]” column from “260327_fits_gcplyr.csv” per growth curve. It represents the maximum growth rate attained by the culture.
      - odmax: maximum value of the chosen “movavg_[number]” column from “260327_fits_gcplyr.csv” per growth curve. It represents the maximum OD620 value attained by the culture, or carrying capacity.
      - auc: Integrated values of movavg_[number] over time per growth curve. Given its nature, it incorporates the values of r and odmax to estimate not only the maximum speed at which a culture grows, but also how quickly that speed is attained and for how long it is sustained.
      - max_percap_time: estimated time from the beginning of kinetic measurements to the moment the maximum growth rate (r) is detected.
      - time_to_stationary: estimated time from the beginning of kinetic measurements to the moment growth rate declines to a value of derivpercap[number]_[number] = 0.5*r, or the onset of the stationary phase.
      - exponential_interval: time from “max_percap_time” to “time_to_stationary”.
   
## 04_TCS_Thermal_performance_curves
- **01_Code:** Has all the R code necessary for fitting thermal performance curves from growth curve parameters.
  - `260327_TCS_TPCs.R`: This file uses the growth curve parameters calculated above as input to fit thermal performance curves for the growth rate (r) and the area under the curve (AUC), using five different models in the package Rtpc. After filtering out unrealistic fits (semi-manual step), it generates a consensus TPC model for each isolate and metric using a weighted average approach, and then extracts the TPC parameters.
    
- **02_Output:** Mainly contains intermediate datasets generated after time-intensive steps, and final output files to be used in other scripts
    - `tpc_input_gcplyr_average.csv`: Intermediate input file after summarizing the growth curve data by isolate and temperature
      - curve_ID: Isolate genome identifier (required by rTPC and equivalent to seqID in other datasets)
      - rate, auc: average value of each growth metric across replicates for each isolate and temperature
      - curve_id: Isolate genome identifier
    - `preds_glu_fits_model_selection_rate_gcplyr.csv.gz` and `preds_glu_fits_model_selection_auc_gcplyr.csv.gz`: Fitted TPCs by model and isolate for the metrics r and AUC
      - curve_id: Isolate genome identifier (equivalent to seqID in other datasets)
      - temp: temperature in Degrees Celsius
      - .fitted: predicted performance at a given temperature
    - `glu_params_rate_gcplyr.csv.gz` and `glu_params_auc_gcplyr.csv.gz`: TPC parameters by model and isolate for the metrics r and AUC
      - curve_id: Isolate genome identifier (equivalent to seqID in other datasets).
      - model_name: name of the model used to predict a given TPC, as provided in the rTPC package.
      - rmax: performance peak for either r or AUC
      - topt: Thermal optimum, or optimal growth temperature estimated by the TPC model for a given isolate and growth metric
      - ctmin: Lower thermal limit. The temperature below which the growth performance is 0 or lower.
      - ctmax: Upper thermal limit. The temperature above which the growth performance is 0 or lower.
    - `ave_preds_glu_rate_gcplyr.csv.gz` and `ave_preds_glu_auc_gcplyr.csv.gz`: Consensus TPC fit by isolate for the metrics r and AUC
      - curve_id: Isolate genome identifier (equivalent to seqID in other datasets)
      - temp: temperature in Degrees Celsius
      - .fitted: predicted performance at a given temperature, averaged across models
      - model_name: this column displays the specific model used to fit the TPC. In this case, it has the value “model average”
    - `ave_params_glu_rate_gcplyr.csv.gz` and `ave_params_glu_auc_gcplyr.csv.gz`: Consensus TPC parameters by isolate for the metrics r and AUC (same columns as the non-averaged file)
    - `251029_tcs_fits_gcplyr.csv.gz`: Condensed TPC fits data frame with values for each metric and isolate
      - seqID: Isolate genome identifier (equivalent to curve_id in other datasets)
      - temp: temperature in Degrees Celsius
      - .fitted: predicted performance at a given temperature, averaged across models
      - model_name: this column displays the specific model used to fit the TPC. In this case, it has the value “model average” 
      - metric: specific growth metric used to fit the TPC. Either r or AUC
    - `251029_tcs_params_gcplyr.csv.gz`
      - Condensed TPC parameter data frame with values for each metric and isolate
      - Relevant column descriptions:
      - seqID: Isolate genome identifier (equivalent to curve_id in other datasets).
      - Glucose_stock_X: glucose stock concentration used for each treatment. A value of 10 results in a 0.007 M-C media, while a value of 100 results in a 0.07 M-C media.
      - model_name: name of the model used to predict a given TPC, as provided in the rTPC package.
      - rmax: performance peak for either r or AUC
      - topt: Thermal optimum, or optimal growth temperature estimated by the TPC model for a given isolate and growth metric
      - ctmin: Lower thermal limit. The temperature below which the growth performance is 0 or lower.
      - ctmax: Upper thermal limit. The temperature above which the growth performance is 0 or lower.
      - metric: specific growth metric used to fit the TPC. Either r or AUC

### 05_TCS_phylogenetic_models
- **01_Code:** Contains two R scripts, the first one to merge taxonomy, thermal traits, and bioclimatic variables with the tree and the second one to perform the bayesian phylogenetic models 
  - `01_Merge_env_phylo_andTPC.R`
  - `02_PhylogeneticModel.R`
- **02_Data:** Contains relevant data (phylogenetic tree). 
  - `TCS_genomes_annotated_tree_without_others-labels.newick`: Tree from the KBase pipeline
- **03_Output:** A series of files integrating thermal performance, taxonomic information, and bioclimatic variables to produce figures in the paper. 
  - `ave_params_gcplyr_combined_common_taxa_auc.csv.gz`, `ave_params_gcplyr_combined_common_taxa_rate.csv.gz`, `ave_preds_gcplyr_combined_common_taxa_auc.csv.gz`, `ave_preds_gcplyr_combined_common_taxa_rate.csv.gz`: Data sets with thermal performance parameters (params) or predictions across temperature values (preds), taxonomic information of isolates, and information about the environment where it came from for the most common groups in our data-set. 
  - `ave_params_gcplyr_combined_major_taxa_auc.csv.gz`, `ave_params_gcplyr_combined_major_taxa_rate.csv.gz`, `ave_preds_gcplyr_combined_major_taxa_auc.csv.gz`, `ave_preds_gcplyr_combined_major_taxa_rate.csv.gz`: Data sets with thermal performance parameters (params) or predictions across temperature values (preds), taxonomic information of isolates, and information about the environment where it came from for Enterobacterales and Pseudomonadales.
  - `tpc_tree_kbase_relaxed.tre`: filtered tree with only the relevant taxa. 
  - `tpc_tree_metadata_auc.csv.gz` and `tpc_tree_metadata_rate.csv.gz`: Associated metadata for the phylogenetic tree, including TPC parameters for rate or auc (depending on the file). 
  

### 06_TCS_Figures
- **01_Code:** Contains a very large R script to generate all the main text figures and 22 out of 25 supplemental figures (the "bioclimatic_variables" scripts generate other two figures, which are still included in the main “Figures” repository), and the "phylogenetic_models" code generates Supplementary Figure 11.
    - `260327_TCS_figures.R`: This script generates most of the output figures by metric, so a big portion of the code, including plots and statistics, is duplicated for AUC and r. Each section of the code is named after the main text or supplemental figure it generates, the statistical analyses performed, and the metric it corresponds to.

- **02_Data:** Contains relevant data (mostly metadata) for the analyses described above.
  - `phylogenetic_distances_rate_kbase.txt`: Matrix of phylogenetic distances between the isolates kept in the phylogenetic tree for the metric r
  - `phylogenetic_distances_auc_kbase.txt`: Matrix of phylogenetic distances between the isolates kept in the phylogenetic tree for the metric AUC
  - `CA_State.shp` and `CA_State.shx`: Shapefile of the state of California
  - `250915_ave_preds_CCA_gcplyr_combined.csv`
    - Predicted TPC fits obtained in a separate experiment to evaluate the influence of carbon concentration on TPC shape
    - Relevant column descriptions:
    - seqID: Isolate genome identifier (equivalent to curve_id in other datasets)
    - Glucose_stock_X: glucose stock concentration used for each treatment. A value of 10 results in a 0.007 M-C media, while a value of 100 results in a 0.07 M-C media.
    - temp: temperature in Degrees Celsius
    - .fitted: predicted performance at a given temperature, averaged across models
    - model_name: this column displays the specific model used to fit the TPC. In this case, it has the value “model average” 
  - `250915_ave_params_CCA_gcplyr_combined.csv`:
    - Predicted TPC parameters obtained in a separate experiment to evaluate the influence of carbon concentration on TPC shape
    - Relevant column descriptions:
    - seqID: Isolate genome identifier (equivalent to curve_id in other datasets).
    - Glucose_stock_X: glucose stock concentration used for each treatment. A value of 10 results in a 0.007 M-C media, while a value of 100 results in a 0.07 M-C media.
    - model_name: name of the model used to predict a given TPC, as provided in the rTPC package.
    - rmax: performance peak for either r or AUC
    - topt: Thermal optimum, or optimal growth temperature estimated by the TPC model for a given isolate and growth metric
    - ctmin: Lower thermal limit. The temperature below which the growth performance is 0 or lower.
    - ctmax: Upper thermal limit. The temperature above which the growth performance is 0 or lower.

- **03_Output**:
  - `combined_effects_lmer_auc_full.csv`: Summary of linear mixed regression models
  - `combined_effects_lmer_auc_simple.csv`: Summary of linear mixed regression models using non-composite bioclimatic variables as explanatory variables
  - `combined_effects_lmer_rate.csv`: Summary of linear mixed regression models for rate
  - `Combined_results_phylosig_df.csv`: 50 iterations of linear mixed regression models for auc, after obtaining stratified samples of Pseudomonadales
  - `fixed_effects_lmer_auc.csv`: Fixed effects of mixed regression models
  - `fixed_effects_lmer_auc_simple.csv`: Fixed effects of linear mixed regression models using non-composite bioclimatic variables as explanatory variables
 - `fixed_effects_lmer_rate.csv`: Fixed of linear mixed regression models for rate
 - `phylogenetic_signal_thermal_traits_auc.csv`: Random effects of mixed regression models
 - `phylogenetic_signal_thermal_traits_auc_simple.csv`: Random effects of linear mixed regression models using non-composite bioclimatic variables as explanatory variables
 - `phylogenetic_signal_thermal_traits_rate.csv`: Random of linear mixed regression models for rate











