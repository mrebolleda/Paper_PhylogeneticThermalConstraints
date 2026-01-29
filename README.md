# Repository for paper "Conserved upper thermal limits and small safety margins in soil copiotrophic bacteria"
### Authors
Corresponding: María Rebolleda-Gómez (mreboll1@uci.edu)  
Ariel Favier (afavier@uci.edu) (data and code curation)
Alejandra Hernández-Terán (ipomeap@gmail.com) 
Celia Symons (csymons@uci.edu) 

### Abstract
One of the key uncertainties in climate change models is how microbes will adapt to rising temperatures. Large-scale comparisons of bacterial thermal performances show a clear boundary between mesophiles and thermophiles. Here, we investigated whether phylogenetic constraints limit the adaptive potential of bacteria to warming soils. Focusing on copiotrophs within Gammaproteobacteria, we found that both thermal optima and upper thermal limits are constrained; variation in these traits decreases above 42°C across the phylogeny, with minimal influence of present-day bioclimatic variables. This, along with the reduced thermal safety margins found in fluctuating hot climates, suggests that many isolates may already be maladapted to local temperature variability. Our findings indicate that these constraints interact with the geometry of thermal performance curves, imposing a trade-off between high-temperature performance and the risk of substantial fitness losses above Topt. Overall, this work underscores potential limits to thermal adaptation and their implications for bacterial fitness.

## Overview
This repository contains data files and R scripts used to analyze the data for the paper. Each folder contains the necessary data for the script (s) to run independently and their outputs. Most data files already contain the necessary subsets of metadata for plotting and wrangling purposes.
Bioclimatic_variables - Contains two R scripts used to analyze climate time series. 

- **Genome_collection** - Contains a short R script to quality filter the assembled genomes from our KBase narrative, along with the taxonomy assignment of the genome collection.
- **Growth_curves** - Contains the aggregated raw optical density data for our growth kinetics, the R script used to process and curate the data, intermediate steps of data filtering, the final output of the script (used as input for the TPC fitting script in the next folder), and a PDF file with plotted individual growth curves.
- **Thermal_performance_curves** - Contains the main script to fit the TPCs from the growth curve parameter file, along with several datasets summarized in different ways for plotting or analytic purposes, and PDF files with plotted TPCs for different growth metrics.
- **Figures** - Contains a very large R script and input data used to make the main and supplementary figures, along with all statistical analyses. 

## Detailed repository layout
### Bioclimatic_variables
- **01_Code:** Has all the R code necessary for characterizing sites and climate profiles based on soil bioclimatic variables.
    - 251029_Environmental_data_processing.R: This script allows us to compile and summarize the climate data from soil time series and precipitation data across sampling sites. A big portion of the script is dedicated to harmonizing measurements from different sites, determining a shared time interval among all sites, and calculating the bioclimatic variables across timescales. We finally compress the multiple variables into two main Principal Components, that explain most of the differences between sites, and use them to generate broader climate categories based on temperature mean and variability.
    - 251001_precipitation_db_comparison.R: This script assesses the interchangeability of the precipitation data obtained from the UCNRS Dendra database and PRISM. Its main output is the supplemental figure S1
      
- **02_Data:** Contains relevant data (mostly metadata) for the analyses described above.
  - UCM_temperature_data.csv
    - Soil 100mm Temperature degC: Hourly recording of temperature values in degrees Celsius at a soil depth of 100 mm
    - STATION_ID: Weather station identifier. In this case, the station “UCM” is the one closest to the UC Merced campus, where our sampling site is located.
  - temperature-data-ucnrs-240502.csv
      - time: timestamp of temperature collection. Soil temperature values are measured every 10 minutes.
      - [Site name] Soil Temp 50 mm Avg degC: Average temperature values in degrees Celsius at a soil depth of 50 mm.
  - TCS-reserves-MOD13A1-061-results.csv and TCS-reserves-2-MOD13A1-061-results.csv: Satellite data for plant growth season calculation by site based on NDVI. Downloaded from NASA AppEARS.
      - ID: Sampling site initials, equivalent to Site_initials in other datasets
      - MOD13A1_061__500m_16_days_NDVI: Normalized Difference Vegetation Index (NDVI) for a pixel size of 500 m/side over 16 days using MODIS VI products.
      - MOD13A1_061__500m_16_days_VI_Quality_MODLAND_Description: Quality assessment of the NDVI value obtained. We use this column to filter out unreliable values due to factors such as cloud cover.
  - rainfall_data_ucnrs_240722.csv
      - time: timestamp of precipitation collection. Precipitation values are measured every 10 minutes.
      - [Site name]-mm: Yearly cumulative precipitation values in mm per site. Measurements start on 10/01 and finish on 09/30 of the next calendar year.
  - rainfall_data_PRISM_240725.csv
      - [Site name]: Monthly cumulative precipitation for a specific sampling site in mm
  - Pinyon_crest_temperature_data.csv
      - Soil Temp 50 mm Avg degC: Average temperature values in degrees Celsius at a soil depth of 50 mm, after applying offset calculations with Air temperature.
  - dendra_data_long.csv: Long-format soil temperature data obtained from the Dendra database for sampling sites except the Merced Vernal Pools and Grassland Reserve and Pinyon Crest/Flats site. *We provide this intermediate dataframe because creating it is a highly time-consuming step.*
      - Soil Temp 50 mm Avg degC: Average temperature values in degrees Celsius at a soil depth of 50 mm.
      - Site: sampling site identifier, equivalent to Site_initials.
  - air-soil-offset-bdc.csv
      - Deep Canyon Air Temp Avg degC: Average air temperature values in degrees Celsius, at 2 m aboveground.
      - Deep Canyon Soil Temp 50 mm Avg degC: Average temperature values in degrees Celsius at a soil depth of 50 mm.
      
- **03_Output:** Mainly contains intermediate datasets generated after time-intensive steps, and final output files to be used in other scripts
    - temperature_database_tcs_240505.csv
    - summarized_climate_data.csv
  - summarized_climate_data_growing_season.csv: This dataframe contains the average soil temperature time series by site, year, and month, but keeping only dates within the estimated plant growth season of each site
  - site_temperature_probability_growing_season.csv: This dataframe assigns a probability for each temperature within a range of temperatures from -5 to 62.5°C, at a 0.1°C resolution, for each sampling site during the plant growth season
  - precipitation_database_tcs_240730.csv: This dataframe contains all the precipitation-related bioclimatic variables across sites
  - PCs_env_noprec.csv: This dataframe contains the values of principal component (PC) loadings (PC1-PC10) obtained while excluding precipitation variables, per site
  - PCs_env_full.csv: This dataframe contains the values of principal component (PC) loadings (PC1-PC10) obtained while including precipitation variables, per site
  - environment_database_tcs_noprec_250512.csv: Contains the condensed temperature bioclimatic variables and PC loadings per site
  - environment_database_tcs_full_250512.csv: Contains the condensed temperature and precipitation bioclimatic variables, and PC loadings per site
  - climate_profile_temperature_probability.csv: This dataframe assigns a probability for each temperature within a range of temperatures from -5 to 62.5°C, at a 0.1°C resolution, for each Climate Profile
  - climate_profile_temperature_probability_growing_season.csv: This dataframe assigns a probability for each temperature within a range of temperatures from -5 to 62.5°C, at a 0.1°C resolution, for each Climate Profile during the plant growth season
 
### Genome_collection
- **01_Code:** Has a short R script to parse out contaminated genomes
  - TCS_genome_quality_assessment.R: This script filters out contaminated genomes using the CheckM output from the KBase pipeline, and also removes them from the quast output
- **02_Data:** Contains relevant data for the analysis described above.
    - TCS_genome_quality.csv
      - Bin name: isolate genome identifier. The number at the beginning of the string is equivalent to the value of “seqID” or “curve_id” in other datasets. The substring “TCS_genomes_[number]_DRAM” refers to the batch of DRAM-annotated genomes to which the genome belonged within the KBase pipeline. 
      - Marker Lineage: CheckM infers the position of a query genome within a reference phylogenetic tree. This “Marker Lineage” is then used to find the fraction of marker genes that should be present in the genome (and with it its “completeness”), and to estimate the level of “contamination” with marker sequences from other clades.
      - A gene is considered “marker” when it is present as a single-copy gene in >97% of genomes within a lineage. Columns 0-5 represent the number of copies/versions associated with each Marker Gene from the inferred Maker Lineage gene set found in the query genome. These are then used to calculate the columns “Completeness” and “Contamination”.
      - Completeness: The values of column “0” is contrasted against the “Marker sets” column to obtain the completeness percentage, as in (“Marker sets” – “0”)/(“Marker sets”). We use this column to only keep highly complete genomes (most of them, because we sequenced isolate populations).
      - Contamination: Genome contamination is estimated from the number of multicopy marker genes identified in each marker set (Parks et al. 2015). We kept genomes with contamination levels below 5% (“Low” and “no detectable contamination”). 
    - Taxonomy_file_kbase.csv: Taxonomy assignment output from GTDB-Tk in the KBase pipeline
      - seqID: Isolate genome identifier
      - Domain, Phylum, Class, Order, Family, Genus, Species: Taxonomic assignment
      - Warnings: this column provides explanations for isolates with incomplete taxonomic assignments (e.g.: lacking Species name)
    - quast_stats_unfiltered.csv
      - curve_id: Isolate genome identifier
      - Other columns: refer to the “Metrics description” in the Quast (v5.3.0) manual
        
- **03_Output:** Contains output files to be used in other scripts
    - Kbase_taxonomy_quality_unfiltered.csv
    - Kbase_taxonomy_quality_filtered.csv
    - quast_stats_filtered.csv

### Growth_curves
- **01_Code:** Has all the R code necessary for fitting growth curves from the raw optical density data
    - 251029_TCS_growthcurves.R
      
- **02_Data:** Contains relevant data (mostly metadata) for the analyses described above.
    - raw_curves_df.csv: Aggregated optical density file over time for each kinetics measurement
      - seqID: Isolate genome identifier
      - abs: optical density measurement (OD620)
      - t: time since start of kinetics measurement (hs)
      - tmp: experimental growth temperature
      - replicate: experimental replicate
      - Isolation.temperature: temperature at which the organism was experimentally isolated (equivalent to “temperature_of_isolation” in other datasets. 
  - y_gcplyr_with_derivatives.csv: Intermediate file used as input to calculate growth curve parameters. This file contains multiple fits per curve, based on possible combinations of smoothing parameters. In this case, we only provide two combinations, but we screened several others, which extended the processing time considerably.
      - seqID: Isolate genome identifier
      - Domain, Phylum, Class, Order, Family, Genus, Species: Taxonomy assignment, inherited from the KBase pipeline.
      - movavg_[number]: width of the sliding window used to smooth out the original “abs” values. A value of “13” means that the smoothing occurred using abs values along 13 time points (~2.17 hs).
      - derivpercap_[number]_[number]: width of the sliding window used to calculate the maximum derivative along the absorbance, or smoothed absorbance, data. A value of derivpercap_9_13 means that the derivative was calculated using a sliding window of 9 values of movavg_13.

- **03_Output and Plots:** Mainly contains intermediate datasets generated after time-intensive steps, and final output files to be used in other scripts
  - 240820_fits_gcplyr.csv
  - 240820_params_gcplyr.csv: Growth curve parameters obtained from the predicted growth curves
      - seqID: Isolate genome identifier
      - Domain, Phylum, Class, Order, Family, Genus, Species: Taxonomy assignment, inherited from the KBase pipeline.
      - t: time since start of kinetics measurement (hs)
      - tmp: experimental growth temperature
      - replicate: experimental replicate
      - temperature_of_isolation: temperature at which the organism was experimentally isolated (equivalent to “Isolation.temperature” in other datasets.
      - Site_initials: Sampling site from which the isolate was obtained
      - lag_time: lag time calculated by gcplyr. This measurement is extremely sensitive to noise and initial conditions, reason why we did not use it for our downstream analyses
      - r: maximum value of the chosen “derivpercap_[number]_[number]” column from “240820_fits_gcplyr.csv” per growth curve. It represents the maximum growth rate attained by the culture.
      - odmax: maximum value of the chosen “movavg_[number]” column from “240820_fits_gcplyr.csv” per growth curve. It represents the maximum OD620 value attained by the culture, or carrying capacity.
      - auc: Integrated values of movavg_[number] over time per growth curve. Given its nature, it incorporates the values of r and odmax to estimate not only the maximum speed at which a culture grows, but also how quickly that speed is attained and for how long it is sustained.
      - max_percap_time: estimated time from the beginning of kinetic measurements to the moment the maximum growth rate (r) is detected.
      - time_to_stationary: estimated time from the beginning of kinetic measurements to the moment growth rate declines to a value of derivpercap[number]_[number] = 0.5*r, or the onset of the stationary phase.
      - exponential_interval: time from “max_percap_time” to “time_to_stationary”.
    - growthcurves_gcplyr_250324.pdf: Growth curve plots, faceted by isolate and temperature, after filtering out contaminations and artifacts
   
## Thermal_performance_curves
- **01_Code:** Has all the R code necessary for fitting thermal performance curves from growth curve parameters.
  - 251029_TCS_TPCs.R: This file uses the growth curve parameters calculated above as input to fit thermal performance curves for the growth rate (r) and the area under the curve (AUC), using five different models in the package Rtpc. After filtering out unrealistic fits (semi-manual step), it generates a consensus TPC model for each isolate and metric using a weighted average approach, and then extracts the TPC parameters.
- **02_Data:** Contains relevant data (mostly metadata) for the analyses described above.
    - tpc_input_gcplyr_average.csv: Intermediate input file after summarizing the growth curve data by isolate and temperature
      - curve_ID: Isolate genome identifier (required by rTPC and equivalent to seqID in other datasets)
      - rate, auc: average value of each growth metric across replicates for each isolate and temperature
      - curve_id: Isolate genome identifier
- **03_Output_and_Plots:** Mainly contains intermediate datasets generated after time-intensive steps, and final output files to be used in other scripts
    - preds_glu_fits_model_selection_rate_gcplyr.csv and preds_glu_fits_model_selection_auc_gcplyr.csv: Fitted TPCs by model and isolate for the metrics r and AUC
      - curve_id: Isolate genome identifier (equivalent to seqID in other datasets)
      - temp: temperature in Degrees Celsius
      - .fitted: predicted performance at a given temperature
    - glu_params_rate_gcplyr.csv and glu_params_auc_gcplyr.csv: TPC parameters by model and isolate for the metrics r and AUC
      - curve_id: Isolate genome identifier (equivalent to seqID in other datasets).
      - model_name: name of the model used to predict a given TPC, as provided in the rTPC package.
      - rmax: performance peak for either r or AUC
      - topt: Thermal optimum, or optimal growth temperature estimated by the TPC model for a given isolate and growth metric
      - ctmin: Lower thermal limit. The temperature below which the growth performance is 0 or lower.
      - ctmax: Upper thermal limit. The temperature above which the growth performance is 0 or lower.
    - ave_preds_glu_rate_gcplyr.csv and ave_preds_glu_auc_gcplyr.csv: Consensus TPC fit by isolate for the metrics r and AUC
      - curve_id: Isolate genome identifier (equivalent to seqID in other datasets)
      - temp: temperature in Degrees Celsius
      - .fitted: predicted performance at a given temperature, averaged across models
      - model_name: this column displays the specific model used to fit the TPC. In this case, it has the value “model average”
    - ave_params_glu_rate_gcplyr.csv and ave_params_glu_auc_gcplyr.csv: Consensus TPC parameters by isolate for the metrics r and AUC (same columns as the non-averaged file)
    - environment_database_tcs_noprec_250512.csv: Contains the condensed temperature bioclimatic variables and PC loadings per site

### Figures
- **01_Code:** Contains a very large R script to generate all the main text figures and 15 out of 17 supplemental figures (the “bioclimatic_variables” scripts generate the other two figures, which are still included in the main “Figures” repository)
    - 251029_TCS_figures.R: This script generates most of the output figures by metric, so a big portion of the code, including plots and statistics, is duplicated for AUC and r. Each section of the code is named after the main text or supplemental figure it generates, the statistical analyses performed, and the metric it corresponds to.

- **02_Data:** Contains relevant data (mostly metadata) for the analyses described above.
  - tpc_tree_kbase_relaxed.tre: Output phylogenetic tree from the KBase pipeline
  - tpc_tree_kbase_metadata_rate.csv
  - tpc_tree_kbase_metadata_auc.csv
  - phylogenetic_distances_rate_kbase.txt
  - CA_State.shp: Shapefile of the state of California
  - ave_preds_CCA_auc_gcplyr.csv: Predicted TPC fits obtained in a separate experiment to evaluate the influence of carbon concentration on TPC shape
    - curve_id: Isolate genome identifier (equivalent to seqID in other datasets)
    - Glucose_stock_X: glucose stock concentration used for each treatment. A value of 10 results in a 0.007 M-C media, while a value of 100 results in a 0.07 M-C media.
    - temp: temperature in Degrees Celsius
    - .fitted: predicted performance at a given temperature, averaged across models
    - model_name: this column displays the specific model used to fit the TPC. In this case, it has the value “model average” 
  - ave_params_CCA_auc_gcplyr.csv: Predicted TPC parameters obtained in a separate experiment to evaluate the influence of carbon concentration on TPC shape
    - rmax: performance peak for either r or AUC
    - topt: Thermal optimum, or optimal growth temperature estimated by the TPC model for a given isolate and growth metric
    - ctmin: Lower thermal limit. The temperature below which the growth performance is 0 or lower.
    - ctmax: Upper thermal limit. The temperature above which the growth performance is 0 or lower.










