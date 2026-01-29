
#
#
# This is a script to process the growth curve data and fit a thermal performance curve (TPC) by isolate
# We will use the package rTPC, which allows for the use of several TPC models per curve
# The output will include 
# 1) A file with rate, odmax, and auc-derived TPC parameters
# 2) A file with rate, odmax, and auc-derived TPC fits

#### Set the directory and packages ####
setwd("PATH_TO_YOUR_DIRECTORY")

library(data.table)
library(tidyr)
library(rTPC)
library(nls.multstart)
library(broom)
library(dplyr)
library(ggforce)
library(grDevices)
library(MuMIn)
library(purrr)

#### Read growth curve data ####

y_params_gcplyr <- fread("240820_params_gcplyr.csv")
y_fits_gcplyr <- fread("240820_fits_gcplyr.csv")

#### OPTIONAL: Manually remove values from isolates that did not have a good TPC fit ####
#
#
# Here we remove isolates that did not yield biologically realistic TPCs after running the script
# Since fitting the TPCs is a time-consuming process, this step can save time. We include this step if you 
# decide to inspect the removed data 
#
#

y_params_gcplyr <- y_params_gcplyr %>% 
  filter(!seqID %in% c("63","66", "73","74","76","144","145","146","150","152","165","170","172","174","185","188","194","218","234",
                       "253","255","257","259","266","267","268","269","270","274","280","281","284","291","298","312",
                       "313","322","337","342","345","346","352","355","356","361","364","365","366","371","372","381","387", "388","391",
                       "392","396","399"))

y_fits_gcplyr <- y_fits_gcplyr %>% 
  filter(!seqID %in% c("63","66", "73","74","76","144","145","146","150","152","165","170","172","174","185","188","194","218","234",
                       "253","255","257","259","266","267","268","269","270","274","280","281","284","291","298","312",
                       "313","322","337","342","345","346","352","355","356","361","364","365","366","371","372","381","387", "388","391",
                       "392","396","399"))

y_params_gcplyr_metadata <- y_params_gcplyr %>% ungroup() %>% dplyr::select(!c(tmp,replicate,r,odmax,auc,
                                                                               max_percap_time,time_to_stationary,
                                                                               lag_time,exponential_interval,derivative,
                                                                               moving_avg)) %>% distinct()




#### Preparing data for TPC fitting by removing replicates ####


tpc_input_gcplyr<- y_params_gcplyr %>% dplyr::select(!c(replicate)) 

tpc_input_gcplyr <- tpc_input_gcplyr %>% rename("temp" = "tmp",
                                  "curve_id" = "seqID",
                                  "rate" = "r",
                                  "od" = "odmax",
                                  "auc" = "auc")

#### S4 - Identifying outliers with incongruent r and od pairs ####
S4 <- tpc_input_gcplyr %>%
  ggplot(aes(x = rate, y = od, color = temp)) +
  geom_hline(yintercept = 0, 
             linetype = "solid",
             lwd = 0.75,
             color = "white")+
  geom_vline(xintercept = 0, 
             linetype = "solid",
             lwd = 0.75,
             color = "white")+
  geom_hline(yintercept = 0.1, 
             linetype = "dashed",
             lwd = 0.75,
             color = "white")+
  geom_vline(xintercept = 0.1, 
             linetype = "dashed",
             lwd = 0.75,
             color = "white")+
  geom_point(size = 2, alpha = 0.85)+
  scale_color_viridis_c(option = "viridis")+
  ylab("Rate")+
  xlab("ODmax")+
  facet_grid(moving_avg~derivative)+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    strip.background = element_blank(),
    #strip.text = element_blank(),
    strip.placement = NULL,
    panel.background = element_rect(color = "black",fill = "black"),
    legend.position = "right",
    #legend.position = "none",
    legend.direction = "vertical")






#### Obtain average values by curve and smoothing combination ####

# Here we calculate average values across replicates for 
# ODmax (average_od), 
# maximum growth rate (average_rate), 
# Area Under the Curve (average_auc),
# duration of the exponential growth phase (average_exponential_interval),
# time until the maximum growth rate (average_max_percap_time),
# time until stationary phase (average_time_to_stationary),
# lag phase (average_lag_time)

average_rates <- tpc_input_gcplyr %>%
  group_by(curve_id,temp,derivative,moving_avg) %>%
  mutate(od = ifelse(od == 0 & ave(od != 0, curve_id, temp, FUN = sum) >= 2, NA, od),
         rate = ifelse(rate == 0 & ave(rate != 0, curve_id, temp, FUN = sum) >= 2, NA, rate),
         auc = ifelse(auc == 0 & ave(auc != 0, curve_id, temp, FUN = sum) >= 2, NA, auc),
         exponential_interval = ifelse(exponential_interval == 0 & ave(exponential_interval != 0, curve_id, temp, FUN = sum) >= 2, NA, exponential_interval),
         max_percap_time = ifelse(max_percap_time == 0 & ave(max_percap_time != 0, curve_id, temp, FUN = sum) >= 2, NA, max_percap_time),
         time_to_stationary = ifelse(time_to_stationary == 0 & ave(time_to_stationary != 0, curve_id, temp, FUN = sum) >= 2, NA, time_to_stationary),
         lag_time = ifelse(lag_time == 0 & ave(lag_time != 0, curve_id, temp, FUN = sum) >= 2, NA, lag_time)) %>%
  ungroup() %>%
  group_by(curve_id, temp,derivative,moving_avg) %>%
  summarise(average_rate = mean(rate, na.rm = TRUE),
            average_od = mean(od, na.rm = TRUE),
            average_auc = mean(auc, na.rm = TRUE),
            average_exponential_interval = mean(exponential_interval, na.rm = TRUE),
            average_max_percap_time = mean(max_percap_time, na.rm = TRUE),
            average_time_to_stationary = mean(time_to_stationary, na.rm = TRUE),
            average_lag_time = mean(lag_time, na.rm = TRUE)
  ) %>%
  ungroup()

#### Subset smoothed data ####
#
# Here we select the specific combination of derivative and moving average we will be working with.
# The criterion for this is to use moving averages that smooth out very noisy measurements at 
# the beginning of the kinetics and/or at very high temperatures, and to use a sliding window
# to calculate the derivative (growth rate) that is wide enough to remove artifacts but narrower
# than the total time of exponential growth, so we do not underestimate the maximum growth rate
#
#

tpc_input_gcplyr <- average_rates %>% 
  filter (moving_avg == 21 & derivative == 29) %>% 
  left_join(tpc_input_gcplyr, by = c("curve_id", "temp","moving_avg","derivative"))

tpc_input_gcplyr <- tpc_input_gcplyr %>% 
  dplyr::select(-c(rate,od,auc,lag_time,time_to_stationary,exponential_interval,max_percap_time)) %>% 
  distinct() %>% 
  rename("rate" = "average_rate",
         "od" = "average_od",
         "auc" = "average_auc",
         "lag_time" = "average_lag_time")


# Then we subset data to only keep isolates with measurements in at least 6 temperatures. 
# TPC fitting requires at least 6 temperatures to be reliable, but not all of our isolates reached that number
# We measured a total of 13 different temperatures

tpc_input_gcplyr$temp <- factor(tpc_input_gcplyr$temp)
complete_levels <- names(table(tpc_input_gcplyr$curve_id))[table(tpc_input_gcplyr$curve_id) >= length(unique(tpc_input_gcplyr$temp))-7] 

tpc_input_gcplyr <- tpc_input_gcplyr %>% filter(curve_id %in% complete_levels)
tpc_input_gcplyr$temp <- as.numeric(levels(tpc_input_gcplyr$temp))[tpc_input_gcplyr$temp] 
tpc_input_gcplyr$curve_id <- as.character(tpc_input_gcplyr$curve_id)


write.csv(tpc_input_gcplyr, "tpc_input_gcplyr_average.csv",row.names = FALSE)


#### S3 - Check the duration of exponential intervals ####

exponential_interval_db <- tpc_input_gcplyr %>%
  filter(rate > 0 & average_exponential_interval > 0 & average_exponential_interval < 24)

median_exponential_interval <- median(exponential_interval_db$average_exponential_interval)
mean_exponential_interval <- mean(exponential_interval_db$average_exponential_interval)

S3 <- ggplot(exponential_interval_db, aes(average_exponential_interval)) +
  geom_histogram(fill = "white",color = "white",bins = 100)+
  geom_vline(xintercept = mean_exponential_interval, 
             linetype = "solid",
             lwd = 1.5,
             color = "#DE5925")+
  geom_vline(xintercept = 4.83, 
             linetype = "dashed",
             lwd = 1.5,
             color = "#0090B5")+
  ylab("Count")+
  xlab("Exponential Phase Duration [hs]")+
  theme_classic(base_size = 18, base_family = "Atkinson Hyperlegible Next VF Light") +
  theme(
    strip.background = element_blank(),
    #strip.text = element_blank(),
    strip.placement = NULL,
    panel.background = element_rect(color = "black",fill = "black"),
    legend.position = "right",
    #legend.position = "none",
    legend.direction = "vertical")









#### Create custom functions for TPC fitting and filtering ####

# Custom function to show the progress of TPC calculations
nls_multstart_progress <- function(formula, data = parent.frame(), iter, start_lower, 
                                   start_upper, supp_errors = c("Y", "N"), convergence_count = 100, 
                                   control, modelweights, pb = NULL, ...){
  if (!is.null(pb)) {
    pb$tick()
  }
  nls_multstart(formula = formula, data = data, iter = iter, start_lower = start_lower, 
                start_upper = start_upper, supp_errors = supp_errors, convergence_count = convergence_count, 
                control = control, modelweights = modelweights, ...)
}

# Custom function to calculate weights by AIC when AICc cannot be calculated
# This function will be useful when comparing the fits of different models to the data of an individual isolate
# It is not realistic to assume that a single model will be the best fit for fit all bacterial TPCs, 
# and thus we will obtain weighted averages of five different models for each isolate

calculate_weights <- function(aicc, aic) {
  if (all(is.infinite(aicc))) {
    return(MuMIn::Weights(aic))
  } else {
    return(MuMIn::Weights(aicc))
  }
}



# Specify number of curves and models to be applied. We tried applying all models available in 
# the rTPC package, but only 5 phenomenological models fitted all the data 

number_of_models <- 5
number_of_curves <- length(unique(tpc_input_gcplyr$curve_id))


###### Sequentially fit models for three different metrics of growth ######
#### RATE - Fit 5 different models (Runtime for 342 curves: 22 min) ####
tpc_input_gcplyr <- fread("tpc_input_gcplyr_average.csv")


# Set up progress bar
pb <- progress::progress_bar$new(total = number_of_curves * number_of_models,
                                 clear = FALSE,
                                 format = "[:bar] :percent :elapsedfull")


# We fit five chosen model formulations in rTPC. 
# The output contains a column with a list of TPC results by model, that we need to then unnest

tpc_model_selection_rate_gcplyr <- nest(
  tpc_input_gcplyr,
  data = c(temp, rate, od, auc, lag_time,average_exponential_interval,average_time_to_stationary,average_max_percap_time)
) %>% mutate(
  briere2 = purrr::map(data, ~{
    tryCatch(
      {
        nls_multstart_progress(
          rate ~ briere2_1999(temp = temp, tmin, tmax, a, b),
          data = .x,
          iter = c(4, 4, 4, 4),
          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 1,
          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 1,
          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
          supp_errors = 'Y',
          convergence_count = FALSE,
          pb = pb
        )
      },
      error = function(e) {
        message("Error occurred in briere2 model fitting: ", e$message)
        NULL
      }
    )
  }),
  gaussian = purrr::map(data, ~{
    tryCatch(
      {
        nls_multstart_progress(
          rate ~ gaussian_1987(temp = temp, rmax, topt, a),
          data = .x,
          iter = c(4, 4, 4),
          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
          supp_errors = 'Y',
          convergence_count = FALSE,
          pb = pb
        )
      },
      error = function(e) {
        message("Error occurred in gaussian model fitting: ", e$message)
        NULL
      }
    )
  }),
  rezende = purrr::map(data, ~{
    tryCatch(
      {
        nls_multstart_progress(
          rate ~ rezende_2019(temp = temp, q10, a, b, c),
          data = .x,
          iter = c(4, 4, 4, 4),
          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') - 10,
          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') + 10,
          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
          supp_errors = 'Y',
          convergence_count = FALSE,
          pb = pb
        )
      },
      error = function(e) {
        message("Error occurred in rezende model fitting: ", e$message)
        NULL
      }
    )
  }),
  spain = purrr::map(data, ~{
    tryCatch(
      {
        nls_multstart_progress(
          rate ~ spain_1982(temp = temp, a, b, c, r0),
          data = .x,
          iter = c(3, 3, 3, 3),
          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 10,
          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 10,
          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
          supp_errors = 'Y',
          convergence_count = FALSE,
          pb = pb
        )
      },
      error = function(e) {
        message("Error occurred in spain model fitting: ", e$message)
        NULL
      }
    )
  }),
  lactin2 = purrr::map(data, ~{
    tryCatch(
      {
        nls_multstart_progress(
          rate ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
          data = .x,
          iter = c(3, 3, 3, 3),
          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') - 10,
          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') + 10,
          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
          supp_errors = 'Y',
          convergence_count = FALSE,
          pb = pb
        )
      },
      error = function(e) {
        message("Error occurred in lactin2 model fitting: ", e$message)
        NULL
      }
    )
  })
)  

#### Stack models ####
d_stack_glu_rate_gcplyr <- dplyr::select(tpc_model_selection_rate_gcplyr, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', briere2:lactin2)

# Save the dataframe
saveRDS(d_stack_glu_rate_gcplyr, "d_stack_glu_rate_gcplyr.rds")
d_stack_glu_rate_gcplyr <- readRDS("d_stack_glu_rate_gcplyr.rds")


# Remove values of fit that are NULL (fit is a vector of lists)

d_stack_glu_rate_gcplyr <- d_stack_glu_rate_gcplyr %>%
  filter(map_lgl(fit, ~ !is.null(.x)))


glu_params_rate_gcplyr <- d_stack_glu_rate_gcplyr %>%
  mutate(., params = purrr::map(fit, calc_params)) %>%
  dplyr::select(-fit) %>%
  unnest(params)

write.csv(glu_params_rate_gcplyr, "glu_params_rate_gcplyr.csv",row.names = FALSE)

#### Calculate estimated parameters ####

glu_params_rate_gcplyr <- fread("glu_params_rate_gcplyr.csv")

# And then filter out fits based on the relationship between TPC breadth and peak. 
# Some models have the tendency to produce artifacts by overestimating the peak values of very low curves, 
# This is evidenced by very high peak and low breadth values

horrible_fits <- glu_params_rate_gcplyr %>%
  filter(rmax > 1 & breadth < 6 | rmax > 2 | breadth < 3)

glu_params_rate_gcplyr <- glu_params_rate_gcplyr %>%
  filter(!(rmax > 1 & breadth < 6 | rmax > 2 | breadth < 3))

d_stack_glu_rate_gcplyr <- d_stack_glu_rate_gcplyr %>%
  anti_join(horrible_fits, by = c("curve_id", "model_name"))

#### Obtain fit predictions for the range of temperatures measured after removing questionable fits ####
# Here we predict the trajectory of the TPCs based on the fit of each model before calculating the weighted averages

newdata <- tibble(temp = seq(min(tpc_input_gcplyr$temp), max(tpc_input_gcplyr$temp), length.out = 100))

d_preds_glu_rate_gcplyr <- d_stack_glu_rate_gcplyr %>%
  mutate(preds = purrr::map(fit, ~{
    tryCatch(augment(.x, newdata = newdata), 
             error = function(e) NULL)
  })) %>%
  filter(!is.null(preds)) %>%
  dplyr::select(-fit) %>%
  unnest(preds)

write.csv(d_preds_glu_rate_gcplyr,"preds_glu_fits_model_selection_rate_gcplyr.csv",row.names = FALSE)

#### Obtain AICc from each model and curve after removing null values and calculate averages ####

d_stack_glu_rate_gcplyr <- d_stack_glu_rate_gcplyr %>%
  mutate(error_encountered = grepl("Improper input parameters.", fit))

# We first flag and remove fitting errors 
d_stack_glu_rate_gcplyr <- d_stack_glu_rate_gcplyr %>% filter (error_encountered != TRUE)
null_rows <- sapply(d_stack_glu_rate_gcplyr$fit, function(x) is.null(x))
d_stack_glu_rate_gcplyr <- d_stack_glu_rate_gcplyr[!null_rows, ]

# We calculate AIC values
d_ic_glu_rate_gcplyr <- d_stack_glu_rate_gcplyr %>%
  mutate(., info = purrr::map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  dplyr::select(-fit) %>%
  unnest(info) %>%
  dplyr::select(curve_id, model_name, sigma, AIC, AICc, BIC, df.residual)

# And get weights by model based on fits
d_ic_glu_rate_gcplyr <- d_ic_glu_rate_gcplyr %>%
  group_by(curve_id) %>%
  mutate(weight = calculate_weights(AICc, AIC))

# We finally obtain average parameters based on model weights
ave_params_glu_rate_gcplyr <- left_join(glu_params_rate_gcplyr, dplyr::select(d_ic_glu_rate_gcplyr, curve_id, model_name, weight)) %>%
  group_by(curve_id)%>%
  summarise(across(rmax:skewness, ~ {
    filtered_values <- .[!is.infinite(.)]
    sum(filtered_values * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE)
  }))%>% 
  filter (!is.na(rmax))%>%
  mutate(model_name = 'model average')

# We finally obtain average fits based on model weights
ave_preds_glu_rate_gcplyr <- left_join(d_preds_glu_rate_gcplyr, dplyr::select(d_ic_glu_rate_gcplyr, curve_id, model_name, weight)) %>%
  mutate(across(c(.fitted, weight), as.numeric)) %>%  # Convert needed columns to numeric formate
  filter(!is.na(.fitted) & !is.na(weight)) %>%  # Remove rows with NA values in .fitted or weight
  group_by(curve_id, temp) %>%
  summarise(
    .fitted = sum(.fitted * weight),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(model_name = 'model average')%>% 
  filter (curve_id %in% ave_params_glu_rate_gcplyr$curve_id)

write.csv(ave_preds_glu_rate_gcplyr,"ave_preds_glu_rate_gcplyr.csv",row.names = FALSE)  
write.csv(ave_params_glu_rate_gcplyr,"ave_params_glu_rate_gcplyr.csv",row.names = FALSE)  



#### AUC - Fit 5 different models (Runtime for 331 curves: 17 min) ####
tpc_input_gcplyr <- fread("tpc_input_gcplyr_average.csv")


# Set up progress bar
pb <- progress::progress_bar$new(total = number_of_curves * number_of_models,
                                 clear = FALSE,
                                 format = "[:bar] :percent :elapsedfull")


# We fit five chosen model formulations in rTPC. 
# The output contains a column with a list of TPC results by model, that we need to then unnest
tpc_model_selection_auc_gcplyr <- nest(
  tpc_input_gcplyr,
  data = c(temp, rate, od, auc, lag_time,average_exponential_interval,average_time_to_stationary,average_max_percap_time)
) %>% mutate(
  briere2 = purrr::map(data, ~{
    tryCatch(
      {
        nls_multstart_progress(
          auc ~ briere2_1999(temp = temp, tmin, tmax, a, b),
          data = .x,
          iter = c(4, 4, 4, 4),
          start_lower = get_start_vals(.x$temp, .x$auc, model_name = 'briere2_1999') - 1,
          start_upper = get_start_vals(.x$temp, .x$auc, model_name = 'briere2_1999') + 1,
          lower = get_lower_lims(.x$temp, .x$auc, model_name = 'briere2_1999'),
          upper = get_upper_lims(.x$temp, .x$auc, model_name = 'briere2_1999'),
          supp_errors = 'Y',
          convergence_count = FALSE,
          pb = pb
        )
      },
      error = function(e) {
        message("Error occurred in briere2 model fitting: ", e$message)
        NULL
      }
    )
  }),
  gaussian = purrr::map(data, ~{
    tryCatch(
      {
        nls_multstart_progress(
          auc ~ gaussian_1987(temp = temp, rmax, topt, a),
          data = .x,
          iter = c(4, 4, 4),
          start_lower = get_start_vals(.x$temp, .x$auc, model_name = 'gaussian_1987') - 10,
          start_upper = get_start_vals(.x$temp, .x$auc, model_name = 'gaussian_1987') + 10,
          lower = get_lower_lims(.x$temp, .x$auc, model_name = 'gaussian_1987'),
          upper = get_upper_lims(.x$temp, .x$auc, model_name = 'gaussian_1987'),
          supp_errors = 'Y',
          convergence_count = FALSE,
          pb = pb
        )
      },
      error = function(e) {
        message("Error occurred in gaussian model fitting: ", e$message)
        NULL
      }
    )
  }),
  rezende = purrr::map(data, ~{
    tryCatch(
      {
        nls_multstart_progress(
          auc ~ rezende_2019(temp = temp, q10, a, b, c),
          data = .x,
          iter = c(4, 4, 4, 4),
          start_lower = get_start_vals(.x$temp, .x$auc, model_name = 'rezende_2019') - 10,
          start_upper = get_start_vals(.x$temp, .x$auc, model_name = 'rezende_2019') + 10,
          lower = get_lower_lims(.x$temp, .x$auc, model_name = 'rezende_2019'),
          upper = get_upper_lims(.x$temp, .x$auc, model_name = 'rezende_2019'),
          supp_errors = 'Y',
          convergence_count = FALSE,
          pb = pb
        )
      },
      error = function(e) {
        message("Error occurred in rezende model fitting: ", e$message)
        NULL
      }
    )
  }),
  spain = purrr::map(data, ~{
    tryCatch(
      {
        nls_multstart_progress(
          auc ~ spain_1982(temp = temp, a, b, c, r0),
          data = .x,
          iter = c(3, 3, 3, 3),
          start_lower = get_start_vals(.x$temp, .x$auc, model_name = 'spain_1982') - 10,
          start_upper = get_start_vals(.x$temp, .x$auc, model_name = 'spain_1982') + 10,
          lower = get_lower_lims(.x$temp, .x$auc, model_name = 'spain_1982'),
          upper = get_upper_lims(.x$temp, .x$auc, model_name = 'spain_1982'),
          supp_errors = 'Y',
          convergence_count = FALSE,
          pb = pb
        )
      },
      error = function(e) {
        message("Error occurred in spain model fitting: ", e$message)
        NULL
      }
    )
  }),
  lactin2 = purrr::map(data, ~{
    tryCatch(
      {
        nls_multstart_progress(
          auc ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
          data = .x,
          iter = c(3, 3, 3, 3),
          start_lower = get_start_vals(.x$temp, .x$auc, model_name = 'lactin2_1995') - 10,
          start_upper = get_start_vals(.x$temp, .x$auc, model_name = 'lactin2_1995') + 10,
          lower = get_lower_lims(.x$temp, .x$auc, model_name = 'lactin2_1995'),
          upper = get_upper_lims(.x$temp, .x$auc, model_name = 'lactin2_1995'),
          supp_errors = 'Y',
          convergence_count = FALSE,
          pb = pb
        )
      },
      error = function(e) {
        message("Error occurred in lactin2 model fitting: ", e$message)
        NULL
      }
    )
  })
)  


#### Stack models ####
d_stack_glu_auc_gcplyr <- dplyr::select(tpc_model_selection_auc_gcplyr, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', briere2:lactin2)

# Save the dataframe
saveRDS(d_stack_glu_auc_gcplyr, "d_stack_glu_auc_gcplyr.rds")
d_stack_glu_auc_gcplyr <- readRDS("d_stack_glu_auc_gcplyr.rds")

# We remove values of fit that are NULL (fit is a vector of lists)
d_stack_glu_auc_gcplyr <- d_stack_glu_auc_gcplyr %>%
  filter(map_lgl(fit, ~ !is.null(.x)))

glu_params_auc_gcplyr <- d_stack_glu_auc_gcplyr %>%
  mutate(., params = purrr::map(fit, calc_params)) %>%
  dplyr::select(-fit) %>%
  unnest(params)

write.csv(glu_params_auc_gcplyr, "glu_params_auc_gcplyr.csv",row.names = FALSE)

#### Calculate estimated parameters ####

glu_params_auc_gcplyr <- fread("glu_params_auc_gcplyr.csv")

hist(glu_params_auc_gcplyr$rmax)
hist(glu_params_auc_gcplyr$breadth)

glu_params_auc_gcplyr %>% 
  ggplot(aes(x = breadth, y = rmax))+
  geom_point()+
  facet_wrap(~model_name)+
  theme_classic()

# And then filter out fits based on the relationship between TPC breadth and peak. 
# Some models have the tendency to produce artifacts by overestimating the peak values of very low curves, 
# This is evidenced by very high peak and low breadth values

horrible_fits <- glu_params_auc_gcplyr %>%
  filter(rmax > 2 & breadth < 6 | rmax > 7.5 | breadth < 3)

glu_params_auc_gcplyr <- glu_params_auc_gcplyr %>%
  filter(!(rmax > 2 & breadth < 6 | rmax > 7.5 | breadth < 3))

d_stack_glu_auc_gcplyr <- d_stack_glu_auc_gcplyr %>%
  anti_join(horrible_fits, by = c("curve_id", "model_name"))

#### Obtain fit predictions for the range of temperatures measured after removing questionable fits ####
# Here we predict the trajectory of the TPCs based on the fit of each model before calculating the weighted averages
newdata <- tibble(temp = seq(min(tpc_input_gcplyr$temp), max(tpc_input_gcplyr$temp), length.out = 100))

d_preds_glu_auc_gcplyr <- d_stack_glu_auc_gcplyr %>%
  mutate(preds = purrr::map(fit, ~{
    tryCatch(augment(.x, newdata = newdata), 
             error = function(e) NULL)
  })) %>%
  filter(!is.null(preds)) %>%
  dplyr::select(-fit) %>%
  unnest(preds)

write.csv(d_preds_glu_auc_gcplyr,"preds_glu_fits_model_selection_auc_gcplyr.csv",row.names = FALSE)

#### Obtain AICc from each model and curve after removing null values and calculate averages ####

d_stack_glu_auc_gcplyr <- d_stack_glu_auc_gcplyr %>%
  mutate(error_encountered = grepl("Improper input parameters.", fit))

# We first flag and remove fitting errors 
d_stack_glu_auc_gcplyr <- d_stack_glu_auc_gcplyr %>% filter (error_encountered != TRUE)
null_rows <- sapply(d_stack_glu_auc_gcplyr$fit, function(x) is.null(x))
d_stack_glu_auc_gcplyr <- d_stack_glu_auc_gcplyr[!null_rows, ]

# We calculate AIC values
d_ic_glu_auc_gcplyr <- d_stack_glu_auc_gcplyr %>%
  mutate(., info = purrr::map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  dplyr::select(-fit) %>%
  unnest(info) %>%
  dplyr::select(curve_id, model_name, sigma, AIC, AICc, BIC, df.residual)

# And get weights by model based on fits
d_ic_glu_auc_gcplyr <- d_ic_glu_auc_gcplyr %>%
  group_by(curve_id) %>%
  mutate(weight = calculate_weights(AICc, AIC))

# We finally obtain average parameters based on model weights
ave_params_glu_auc_gcplyr <- left_join(glu_params_auc_gcplyr, dplyr::select(d_ic_glu_auc_gcplyr, curve_id, model_name, weight)) %>%
  group_by(curve_id)%>%
  summarise(across(rmax:skewness, ~ {
    filtered_values <- .[!is.infinite(.)]
    sum(filtered_values * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE)
  }))%>% 
  filter (!is.na(rmax))%>%
  mutate(model_name = 'model average')

# We finally obtain average fits based on model weights
ave_preds_glu_auc_gcplyr <- left_join(d_preds_glu_auc_gcplyr, dplyr::select(d_ic_glu_auc_gcplyr, curve_id, model_name, weight)) %>%
  mutate(across(c(.fitted, weight), as.numeric)) %>%  # Convert needed columns to numeric formate
  filter(!is.na(.fitted) & !is.na(weight)) %>%  # Remove rows with NA values in .fitted or weight
  group_by(curve_id, temp) %>%
  summarise(
    .fitted = sum(.fitted * weight),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(model_name = 'model average')%>% 
  filter (curve_id %in% ave_params_glu_auc_gcplyr$curve_id)

write.csv(ave_preds_glu_auc_gcplyr,"ave_preds_glu_auc_gcplyr.csv",row.names = FALSE)  
write.csv(ave_params_glu_auc_gcplyr,"ave_params_glu_auc_gcplyr.csv",row.names = FALSE)  


#### Add metadata to TPC fits #### 

gcplyr_results <- list(ave_preds_glu_rate_gcplyr,
                         ave_preds_glu_auc_gcplyr,
                         ave_params_glu_rate_gcplyr,
                         ave_params_glu_auc_gcplyr)

gcplyr_results <- gcplyr_results %>%
  purrr::map(~ .x %>%
               rename(seqID = curve_id) %>%
               left_join(y_params_gcplyr_metadata, by = "seqID")
  )

list2env(
  gcplyr_results,
  envir = .GlobalEnv
)


#### Remove bad fits (Skip if you have done this in the first step) ####

# Bad fits have unrealistic TPC shapes (e.g. they fit a peak performance that greatly exceeds the maximum growth measured)

ave_preds_glu_rate_gcplyr <- ave_preds_glu_rate_gcplyr %>% 
  filter(!seqID %in% c("63","66", "73","76","103","144","145","146","150","152","165","170","172","174","182","185","188","194","218","234",
                       "253","255","257","259","266","267","268","269","270","274","280","281","284","291","298","312",
                       "313","322","337","342","345","346","352","355","356","361","364","365","366","371","372","381","387", "388","391",
                       "392","396","399"))

ave_preds_glu_auc_gcplyr <- ave_preds_glu_auc_gcplyr %>% 
  filter(!seqID %in% c("63","66", "73","76","103","144","145","146","150","152","165","170","172","174","182","185","188","194","218","234",
                       "253","255","257","259","266","267","268","269","270","274","280","281","284","291","298","312",
                       "313","322","337","342","345","346","352","355","356","361","364","365","366","371","372","381","387", "388","391",
                       "392","396","399"))

ave_params_glu_rate_gcplyr <- ave_params_glu_rate_gcplyr %>% 
  filter(!seqID %in% c("63","66", "73","76","103","144","145","146","150","152","165","170","172","174","182","185","188","194","218","234",
                       "253","255","257","259","266","267","268","269","270","274","280","281","284","291","298","312",
                       "313","322","337","342","345","346","352","355","356","361","364","365","366","371","372","381","387", "388","391",
                       "392","396","399"))

ave_params_glu_auc_gcplyr <- ave_params_glu_auc_gcplyr %>% 
  filter(!seqID %in% c("63","66", "73","76","103","144","145","146","150","152","165","170","172","174","182","185","188","194","218","234",
                       "253","255","257","259","266","267","268","269","270","274","280","281","284","291","298","312",
                       "313","322","337","342","345","346","352","355","356","361","364","365","366","371","372","381","387", "388","391",
                       "392","396","399"))


#### Condensing all TPC metrics #### 
# We create a single file containing the TPC parameters obtained from each metric

ave_params_glu_rate_gcplyr$metric <- "rate"
ave_params_glu_auc_gcplyr$metric <- "auc"

ave_params_gcplyr_combined <- bind_rows(ave_params_glu_rate_gcplyr,
                                   ave_params_glu_auc_gcplyr)


ave_preds_glu_rate_gcplyr$metric <- "rate"
ave_preds_glu_auc_gcplyr$metric <- "auc"

ave_preds_gcplyr_combined <- bind_rows(ave_preds_glu_rate_gcplyr,
                                   ave_preds_glu_auc_gcplyr)


#### Adding environmental variables and genomic info ####


environmental_database <- fread("environment_database_tcs_noprec_250512.csv")

# We add environmental data to the resulting dataframes

ave_params_gcplyr_combined <- left_join(ave_params_gcplyr_combined,environmental_database)
ave_preds_gcplyr_combined <- left_join(ave_preds_gcplyr_combined,environmental_database)

ave_preds_glu_rate_gcplyr <- ave_preds_gcplyr_combined %>% filter(metric == "rate")
ave_preds_glu_auc_gcplyr <- ave_preds_gcplyr_combined %>% filter(metric == "auc")
ave_params_glu_rate_gcplyr <- ave_params_gcplyr_combined %>% filter(metric == "rate")
ave_params_glu_auc_gcplyr <- ave_params_gcplyr_combined %>% filter(metric == "auc")

# We add genomic data to the resulting dataframes

gc_perc_data <- fread("quast_stats_filtered.csv") %>%
  rename("seqID" = "curve_id")


ave_params_gcplyr_combined <- left_join(ave_params_gcplyr_combined,gc_perc_data)
ave_preds_gcplyr_combined <- left_join(ave_preds_gcplyr_combined,gc_perc_data)


#### Write output files ####

write.csv(ave_preds_glu_rate_gcplyr, "251029_ave_preds_glu_rate_gcplyr.csv",row.names = FALSE)
write.csv(ave_params_glu_rate_gcplyr, "251029_ave_params_glu_rate_gcplyr.csv",row.names = FALSE)
write.csv(ave_preds_glu_auc_gcplyr, "251029_ave_preds_glu_auc_gcplyr.csv",row.names = FALSE)
write.csv(ave_params_glu_auc_gcplyr, "251029_ave_params_glu_auc_gcplyr.csv",row.names = FALSE)

write.csv(ave_params_gcplyr_combined,"251029_tcs_params_gcplyr.csv",row.names = FALSE)
write.csv(ave_preds_gcplyr_combined,"251029_tcs_fits_gcplyr.csv",row.names = FALSE)


#### Plot TPC fits by metric ####

d_preds_glu_rate_gcplyr <- fread("preds_glu_fits_model_selection_rate_gcplyr.csv")
d_preds_glu_auc_gcplyr <- fread("preds_glu_fits_model_selection_rate_gcplyr.csv")


ave_preds_glu_auc_gcplyr <- fread("251029_ave_preds_glu_auc_gcplyr.csv")
ave_params_glu_auc_gcplyr <- fread("251029_ave_params_glu_auc_gcplyr.csv")
ave_preds_glu_rate_gcplyr <- fread("251029_ave_preds_glu_rate_gcplyr.csv")
ave_params_glu_rate_gcplyr <- fread("251029_ave_params_glu_rate_gcplyr.csv")

tpc_input_gcplyr <- fread("tpc_input_gcplyr_average.csv")



gcplyr_results_2 <- list(ave_preds_glu_rate_gcplyr,
                       ave_preds_glu_auc_gcplyr,
                       ave_params_glu_rate_gcplyr,
                       ave_params_glu_auc_gcplyr)

gcplyr_results_2 <- gcplyr_results_2 %>%
  purrr::map(~ .x %>%
               rename(curve_id = seqID) 
  )

list2env(
  gcplyr_results_2,
  envir = .GlobalEnv
)



pdf("~/average_glu_rate_gcplyr.pdf")
for(i in 1:100){
  print(ggplot(d_preds_glu_rate_gcplyr, aes(temp, .fitted)) +
          geom_line(aes(col = model_name),alpha = 0.3) +
          geom_line(data = ave_preds_glu_rate_gcplyr, col = 'blue') +
          geom_point(aes(temp, rate),tpc_input_gcplyr)+
          facet_wrap_paginate(~curve_id, scales = 'free_y',nrow = 5,ncol=5, page = i)+
          theme(
            legend.position = 'none',
            panel.background = element_rect(fill='transparent'), #transparent panel bg
            plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank(), #remove minor gridlines
            legend.background = element_rect(fill='transparent'), #transparent legend bg
            legend.box.background = element_rect(fill='transparent'))+
          coord_cartesian(ylim = c(0, 2))+
          geom_hline(yintercept = 1.5)
        
  )
}



pdf("~/average_glu_auc_gcplyr.pdf")
for(i in 1:100){
  print(ggplot(d_preds_glu_auc_gcplyr, aes(temp, .fitted)) +
          geom_line(aes(col = model_name),alpha = 0.3) +
          geom_line(data = ave_preds_glu_auc_gcplyr, col = 'blue') +
          geom_point(aes(temp, auc),tpc_input_gcplyr)+
          facet_wrap_paginate(~curve_id, scales = 'free_y',nrow = 5,ncol=5, page = i)+
          theme(
            legend.position = 'none',
            panel.background = element_rect(fill='transparent'), #transparent panel bg
            plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank(), #remove minor gridlines
            legend.background = element_rect(fill='transparent'), #transparent legend bg
            legend.box.background = element_rect(fill='transparent'))+
          coord_cartesian(ylim = c(0, 6))
        
  )
}

graphics.off()


#### S5 - Proof of concept of filtering steps ####

d_stack_glu_auc_gcplyr <- readRDS("d_stack_glu_auc_gcplyr.rds")
glu_params_auc_gcplyr <- fread("glu_params_auc_gcplyr.csv")

# Filter out fits  based on the relationship between TPC breadth and peak

horrible_fits <- glu_params_auc_gcplyr %>%
  filter(rmax > 2 & breadth < 6 | rmax > 7.5 | breadth < 3)

glu_params_auc_gcplyr <- glu_params_auc_gcplyr %>%
  filter(!(rmax > 2 & breadth < 6 | rmax > 7.5 | breadth < 3))

# Remove values of fit that are NULL (fit is a vector of lists)

d_stack_glu_auc_gcplyr <- d_stack_glu_auc_gcplyr %>%
  filter(map_lgl(fit, ~ !is.null(.x)))

# Keep the bad fits
d_stack_glu_auc_gcplyr_bad_fits <- d_stack_glu_auc_gcplyr %>%
  filter(curve_id %in% horrible_fits$curve_id)

# Remove the bad fits
d_stack_glu_auc_gcplyr_bad_fits_corrected <- d_stack_glu_auc_gcplyr_bad_fits %>%
  anti_join(horrible_fits, by = c("curve_id", "model_name"))

# Obtain fit predictions for the range of temperatures measured after removing questionable fits
newdata_bad_fits <- tibble(temp = seq(min(tpc_input_gcplyr$temp), max(tpc_input_gcplyr$temp), length.out = 100))

d_preds_glu_auc_gcplyr_bad_fits <- d_stack_glu_auc_gcplyr_bad_fits %>%
  mutate(preds = purrr::map(fit, ~{
    tryCatch(augment(.x, newdata = newdata_bad_fits), 
             error = function(e) NULL)
  })) %>%
  filter(!is.null(preds)) %>%
  dplyr::select(-fit) %>%
  unnest(preds)

d_preds_glu_auc_gcplyr_bad_fits_corrected <- d_stack_glu_auc_gcplyr_bad_fits_corrected %>%
  mutate(preds = purrr::map(fit, ~{
    tryCatch(augment(.x, newdata = newdata_bad_fits), 
             error = function(e) NULL)
  })) %>%
  filter(!is.null(preds)) %>%
  dplyr::select(-fit) %>%
  unnest(preds)

# Obtain AICc from each model and curve after removing null values and calculate averages

d_stack_glu_auc_gcplyr_bad_fits <- d_stack_glu_auc_gcplyr_bad_fits %>%
  mutate(error_encountered = grepl("Improper input parameters.", fit))

d_stack_glu_auc_gcplyr_bad_fits_corrected <- d_stack_glu_auc_gcplyr_bad_fits_corrected %>%
  mutate(error_encountered = grepl("Improper input parameters.", fit))

# Remove bad fits
d_stack_glu_auc_gcplyr_bad_fits <- d_stack_glu_auc_gcplyr_bad_fits %>% filter (error_encountered != TRUE)
null_rows <- sapply(d_stack_glu_auc_gcplyr_bad_fits$fit, function(x) is.null(x))
d_stack_glu_auc_gcplyr_bad_fits <- d_stack_glu_auc_gcplyr_bad_fits[!null_rows, ]

# Remove bad fits
d_stack_glu_auc_gcplyr_bad_fits_corrected <- d_stack_glu_auc_gcplyr_bad_fits_corrected %>% filter (error_encountered != TRUE)
null_rows <- sapply(d_stack_glu_auc_gcplyr_bad_fits_corrected$fit, function(x) is.null(x))
d_stack_glu_auc_gcplyr_bad_fits_corrected <- d_stack_glu_auc_gcplyr_bad_fits_corrected[!null_rows, ]

# Calculate AIC
d_ic_glu_auc_gcplyr_bad_fits <- d_stack_glu_auc_gcplyr_bad_fits %>%
  mutate(., info = purrr::map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  dplyr::select(-fit) %>%
  unnest(info) %>%
  dplyr::select(curve_id, model_name, sigma, AIC, AICc, BIC, df.residual)

# Calculate AIC
d_ic_glu_auc_gcplyr_bad_fits_corrected <- d_stack_glu_auc_gcplyr_bad_fits_corrected %>%
  mutate(., info = purrr::map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  dplyr::select(-fit) %>%
  unnest(info) %>%
  dplyr::select(curve_id, model_name, sigma, AIC, AICc, BIC, df.residual)

# Get weights by model based on fits
d_ic_glu_auc_gcplyr_bad_fits <- d_ic_glu_auc_gcplyr_bad_fits %>%
  group_by(curve_id) %>%
  mutate(weight = calculate_weights(AICc, AIC))

# Get weights by model based on fits_corrected
d_ic_glu_auc_gcplyr_bad_fits_corrected <- d_ic_glu_auc_gcplyr_bad_fits_corrected %>%
  group_by(curve_id) %>%
  mutate(weight = calculate_weights(AICc, AIC))

# Obtain average parameters based on model weights
ave_params_glu_auc_gcplyr_bad_fits <- left_join(glu_params_auc_gcplyr, dplyr::select(d_ic_glu_auc_gcplyr_bad_fits, curve_id, model_name, weight)) %>%
  group_by(curve_id)%>%
  summarise(across(rmax:skewness, ~ {
    filtered_values <- .[!is.infinite(.)]
    sum(filtered_values * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE)
  }))%>% 
  filter (!is.na(rmax))%>%
  mutate(model_name = 'model average')

# Obtain average parameters based on model weights
ave_params_glu_auc_gcplyr_bad_fits_corrected <- left_join(glu_params_auc_gcplyr, dplyr::select(d_ic_glu_auc_gcplyr_bad_fits_corrected, curve_id, model_name, weight)) %>%
  group_by(curve_id)%>%
  summarise(across(rmax:skewness, ~ {
    filtered_values <- .[!is.infinite(.)]
    sum(filtered_values * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE)
  }))%>% 
  filter (!is.na(rmax))%>%
  mutate(model_name = 'model average')

# Obtain average fits based on model weights
ave_preds_glu_auc_gcplyr_bad_fits <- left_join(d_preds_glu_auc_gcplyr_bad_fits, dplyr::select(d_ic_glu_auc_gcplyr_bad_fits, curve_id, model_name, weight)) %>%
  mutate(across(c(.fitted, weight), as.numeric)) %>%  # Convert needed columns to numeric formate
  filter(!is.na(.fitted) & !is.na(weight)) %>%  # Remove rows with NA values in .fitted or weight
  group_by(curve_id, temp) %>%
  summarise(
    .fitted = sum(.fitted * weight),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(model_name = 'model average')%>% 
  filter (curve_id %in% ave_params_glu_auc_gcplyr_bad_fits$curve_id)

# Obtain average fits_corrected based on model weights
ave_preds_glu_auc_gcplyr_bad_fits_corrected <- left_join(d_preds_glu_auc_gcplyr_bad_fits_corrected, dplyr::select(d_ic_glu_auc_gcplyr_bad_fits_corrected, curve_id, model_name, weight)) %>%
  mutate(across(c(.fitted, weight), as.numeric)) %>%  # Convert needed columns to numeric formate
  filter(!is.na(.fitted) & !is.na(weight)) %>%  # Remove rows with NA values in .fitted or weight
  group_by(curve_id, temp) %>%
  summarise(
    .fitted = sum(.fitted * weight),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(model_name = 'model average')%>% 
  filter (curve_id %in% ave_params_glu_auc_gcplyr_bad_fits$curve_id)


S5A <- glu_params_auc_gcplyr %>% 
  ggplot(aes(x = breadth, y = rmax, color = model_name, group = interaction(curve_id,model_name)))+
  geom_point(size = 3)+
  geom_hline(yintercept = 2,linetype = 'dashed')+
  geom_vline(xintercept = 6,linetype = 'dashed')+
  theme_classic()+
  scale_color_viridis_d(option = "plasma")


S5B.top <- ggplot(data = d_preds_glu_auc_gcplyr_bad_fits 
           %>% filter(curve_id %in% c(219,241,279,284,397,398))
           ,aes(x = temp, y = .fitted,group = curve_id)) +
      geom_line(data = ave_preds_glu_auc_gcplyr_bad_fits
                %>% filter(curve_id %in% c(219,241,279,284,397,398))
                ,col = 'grey50',alpha = 0.6,linewidth = 3) +
      geom_line(aes(group = model_name, col = model_name),alpha = 0.6)+
      geom_point(data = tpc_input_gcplyr 
                 %>% filter(curve_id %in% c(219,241,279,284,397,398))
                 , aes(temp, auc))+
      scale_color_viridis_d(option = "plasma")+
      facet_wrap_paginate(~curve_id, scales = 'free_y',nrow = 1,ncol=6)+
      theme_classic()+
      coord_cartesian(ylim = c(0, 10))
    

S5B.bot <- ggplot(data = d_preds_glu_auc_gcplyr_bad_fits_corrected 
       %>% filter(curve_id %in% c(219,241,279,284,397,398))
       ,aes(x = temp, y = .fitted,group = curve_id)) +
  geom_line(data = ave_preds_glu_auc_gcplyr_bad_fits_corrected 
            %>% filter(curve_id %in% c(219,241,279,284,397,398))
            ,col = 'grey50',alpha = 0.6,linewidth = 3) +
  geom_line(aes(group = model_name, col = model_name),alpha = 0.6)+
  geom_point(data = tpc_input_gcplyr 
             %>% filter(curve_id %in% c(219,241,279,284,397,398))
             , aes(temp, auc))+
  scale_color_viridis_d(option = "plasma")+
  facet_wrap_paginate(~curve_id, scales = 'free_y',nrow = 1,ncol=6)+
  theme_classic()+
  coord_cartesian(ylim = c(0, 10))

S5B <- plot_grid(S5B.top + theme(legend.position = "none"),
                 S5B.bot + theme(legend.position = "none")
                 , ncol = 1, align = "v", rel_heights = c(1,1),axis = "l")

S5 <- plot_grid(S5A,
                S5B + theme(legend.position = "none")
                , ncol = 2, align = "hv", rel_widths = c(1,1),axis = "l")

