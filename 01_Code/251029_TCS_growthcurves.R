# This is a script to process the absorbance data of the TCS collection
# After running it we will obtain 
# 1) growth curve fits
# 2) parameters for each isolate across temperatures
# You can rerun the script starting on section "Process curves", which loads the full absorbance dataset

#### Set directory and get absorbance data ####

setwd("PATH_TO_YOUR_DIRECTORY")


library(data.table)
library(ggplot2)
library(grDevices)
library(stats)
library(stringr)
library(utils)
library(base)
library(dplyr)
library(ggforce)
library(gcplyr)

#### Process curves ####
# Get the data

y_gcplyr <- fread("raw_curves_df.csv")

#### Remove probably contaminated isolates ####


#### Smoothing data to minimize OD outliers ####

y_gcplyr <-
  mutate(group_by(y_gcplyr,seqID,well,tmp,replicate),
         movavg_1 = abs,
         movavg_13 = smooth_data(x = t, y = abs,
                                 sm_method = "moving-average", window_width_n = 13),
         movavg_21 = smooth_data(x = t, y = abs,
                                 sm_method = "moving-average", window_width_n = 23))



#### Calculating derivatives of smoothed data ####
#
# This is a very time-intensive step, so we only provide results for the chosen moving average. 
# In total, we tried combinations of 5 moving averages (abs, 3,7,11,15,21) with their derivative windows
#
#

y_gcplyr_with_derivatives <- 
  mutate(group_by(y_gcplyr,seqID,well,tmp,replicate),
         derivpercap_9_21 = calc_deriv(x = t, y = movavg_21,
                                     percapita = TRUE, blank = 0,
                                     window_width_n = 9),
         derivpercap_13_21 = calc_deriv(x = t, y = movavg_21,
                                     percapita = TRUE, blank = 0,
                                     window_width_n = 13),
         derivpercap_21_21 = calc_deriv(x = t, y = movavg_21,
                                     percapita = TRUE, blank = 0,
                                     window_width_n = 21),
         derivpercap_29_21 = calc_deriv(x = t, y = movavg_21,
                                     percapita = TRUE, blank = 0,
                                     window_width_n = 29)
         )

# We removed measurements of a couple of genomes that could not be assigned taxonomy
y_gcplyr_with_derivatives <- y_gcplyr_with_derivatives %>% 
  filter (!is.na(Domain)) %>% 
  filter(!Domain %in% c("Unclassified Bacteria","Unclassified Archaea"))


#### We then removed outliers with no growth or biologically unrealistic growth ####


# In this section we remove specific values after inspecting the curves across temperatures for each isolate and
# comparing (when available) between replicates

y_gcplyr_with_derivatives <- y_gcplyr_with_derivatives %>%
  dplyr::filter (!(seqID == 10), #
                 !(seqID == 23 & tmp == 33 & replicate == 2), #
                 !(seqID == 24 & tmp == 33 & replicate == 2), #
                 !(seqID == 31 & tmp == 33 & replicate == 2), #
                 !(seqID == 32 & tmp == 23 & replicate == 2), #
                 !(seqID == 34 & tmp == 28 & replicate == 2), #
                 !(seqID == 37 & tmp == 33 & replicate == 2), #
                 !(seqID == 38 & tmp == 33 & replicate == 2), #
                 !(seqID == 42 & tmp == 23 & replicate == 2), #
                 !(seqID == 50 & tmp == 23 & replicate == 2), #
                 !(seqID == 56 & tmp == 33 & replicate == 2),
                 !(seqID == 57 & tmp == 33 & replicate == 2),
                 !(seqID == 59 & tmp == 13 & replicate == 2),
                 !(seqID == 59 & tmp == 19 & replicate == 1 & well == "A09"),
                 !(seqID == 59 & tmp == 21 & replicate == 1 & well == "A10"),
                 !(seqID == 148 & tmp == 15 & replicate == 1),
                 !(seqID == 148 & tmp == 23 & replicate == 1),
                 !(seqID == 150 & tmp == 23 & replicate == 2),
                 !(seqID == 151 & tmp == 23 & replicate == 1),
                 !(seqID == 153 & tmp == 23 & replicate == 1),
                 !(seqID == 154 & tmp == 23 & replicate == 1),
                 !(seqID == 155 & tmp == 23 & replicate == 1),
                 !(seqID == 156 & tmp == 23 & replicate == 1),
                 !(seqID == 160 & tmp == 23 & replicate == 1),
                 !(seqID == 162 & tmp == 23 & replicate == 1),
                 !(seqID == 164 & tmp == 23 & replicate == 1),
                 !(seqID == 166 & tmp == 23 & replicate == 1),
                 !(seqID == 167 & tmp == 23 & replicate == 1),
                 !(seqID == 171 & tmp == 28 & replicate == 2),
                 !(seqID == 171 & tmp == 35),
                 !(seqID == 177 & tmp == 38 & replicate == 2),
                 !(seqID == 181 & tmp == 15 & replicate == 2),
                 !(seqID == 181 & tmp == 35 & replicate == 1),
                 !(seqID == 182 & tmp == 23 & replicate == 2),
                 !(seqID == 182 & tmp == 38 & replicate == 2),
                 !(seqID == 183 & tmp == 23 & replicate == 2),
                 !(seqID == 190 & tmp == 23 & replicate == 2),
                 !(seqID == 195 & tmp == 35),
                 !(seqID == 196 & tmp == 19),
                 !(seqID == 196 & tmp == 21),
                 !(seqID == 198 & tmp == 38 & replicate == 2),
                 !(seqID == 204),
                 !(seqID == 208 & tmp == 43 & replicate == 1),
                 !(seqID == 211 & tmp == 33 & replicate == 3),
                 !(seqID == 213),
                 !(seqID == 214 & tmp == 33 & replicate == 3),
                 !(seqID == 215 & tmp == 33 & replicate == 3),
                 !(seqID == 220 & tmp == 35 & replicate == 1),
                 !(seqID == 221 & tmp == 15 & replicate == 1),
                 !(seqID == 221 & tmp == 28 & replicate == 1),
                 !(seqID == 222 & tmp == 35 & replicate == 1),
                 !(seqID == 223 & tmp == 35 & replicate == 1),
                 !(seqID == 224 & tmp == 33 & replicate == 3),
                 !(seqID == 225 & tmp == 33 & replicate == 3),
                 !(seqID == 226 & tmp == 33 & replicate == 3),
                 !(seqID == 227 & tmp == 33 & replicate == 3),
                 !(seqID == 228 & tmp == 33 & replicate == 3),
                 !(seqID == 229 & tmp == 33 & replicate == 3),
                 !(seqID == 229 & tmp == 35 & replicate == 1),
                 !(seqID == 230 & tmp == 15 & replicate == 2),
                 !(seqID == 231 & tmp == 33 & replicate == 3),
                 !(seqID == 233 & tmp == 33 & replicate == 3),
                 !(seqID == 234 & tmp == 33 & replicate == 3),
                 !(seqID == 236 & tmp == 33 & replicate == 3),
                 !(seqID == 237 & tmp == 33 & replicate == 3),
                 !(seqID == 238 & tmp == 33 & replicate == 3),
                 !(seqID == 239 & tmp == 15 & replicate == 2),
                 !(seqID == 239 & tmp == 33 & replicate == 3),
                 !(seqID == 240 & tmp == 33 & replicate == 3),
                 !(seqID == 241 & tmp == 33 & replicate == 3),
                 !(seqID == 241 & tmp == 43 & replicate == 1),
                 !(seqID == 242 & tmp == 33 & replicate == 3),
                 !(seqID == 243 & tmp == 33 & replicate == 3),
                 !(seqID == 245 & tmp == 15 & replicate == 1),
                 !(seqID == 245 & tmp == 33 & replicate == 1),
                 !(seqID == 246 & tmp == 15 & replicate == 2),
                 !(seqID == 247 & tmp == 15 & replicate == 2),
                 !(seqID == 247 & tmp == 33 & replicate == 3),
                 !(seqID == 247 & tmp == 35 & replicate == 1),
                 !(seqID == 248 & tmp == 15 & replicate == 2),
                 !(seqID == 249 & tmp == 15 & replicate == 2),
                 !(seqID == 251 & tmp == 21 & replicate == 2),
                 !(seqID == 252 & tmp == 21 & replicate == 2),
                 !(seqID == 253 & tmp == 19),
                 !(seqID == 253 & tmp == 21),
                 !(seqID == 253 & tmp == 33 & replicate == 2),
                 !(seqID == 253 & tmp == 35 & replicate == 1),
                 !(seqID == 254 & tmp == 21),
                 !(seqID == 255 & tmp == 15 & replicate == 2),
                 !(seqID == 255 & tmp == 19 & replicate == 2),
                 !(seqID == 255 & tmp == 21 & replicate == 2),
                 !(seqID == 255 & tmp == 33 & replicate == 2),
                 !(seqID == 255 & tmp == 35 & replicate == 2),
                 !(seqID == 256 & tmp == 21 & replicate == 2),
                 !(seqID == 257 & tmp == 19),
                 !(seqID == 257 & tmp == 21),
                 !(seqID == 258 & tmp == 19),
                 !(seqID == 258 & tmp == 21),
                 !(seqID == 259 & tmp == 19 & replicate == 2),
                 !(seqID == 259 & tmp == 21 & replicate == 2),
                 !(seqID == 260 & tmp == 21 & replicate == 2),
                 !(seqID == 260 & tmp == 23 & replicate == 2),
                 !(seqID == 261 & tmp == 21 & replicate == 2),
                 !(seqID == 261 & tmp == 23 & replicate == 2),
                 !(seqID == 262 & tmp == 21),
                 !(seqID == 263 & tmp == 19 & replicate == 1),
                 !(seqID == 263 & tmp == 21),
                 !(seqID == 263 & tmp == 33 & replicate == 1),
                 !(seqID == 264 & tmp == 23 & replicate == 2),
                 !(seqID == 264 & tmp == 35 & replicate == 1),
                 !(seqID == 265 & tmp == 19 & replicate == 2),
                 !(seqID == 265 & tmp == 15 & replicate == 2),
                 !(seqID == 266 & tmp == 19),
                 !(seqID == 266 & tmp == 21),
                 !(seqID == 267 & tmp == 19),
                 !(seqID == 267 & tmp == 21),
                 !(seqID == 266 & tmp == 38 & replicate == 1),
                 !(seqID == 267 & tmp == 15 & replicate == 2),
                 !(seqID == 268 & tmp == 19),
                 !(seqID == 268 & tmp == 21),
                 !(seqID == 268 & tmp == 33 & replicate == 2),
                 !(seqID == 268 & tmp == 35),
                 !(seqID == 269 & tmp == 15 & replicate == 2),
                 !(seqID == 269 & tmp == 19),
                 !(seqID == 269 & tmp == 21),
                 !(seqID == 269 & tmp == 35),
                 !(seqID == 270 & tmp == 19),
                 !(seqID == 270 & tmp == 21),
                 !(seqID == 271 & tmp == 21),
                 !(seqID == 271 & tmp == 33 & replicate == 2),
                 !(seqID == 272 & tmp == 19),
                 !(seqID == 272 & tmp == 21 & replicate == 2),
                 !(seqID == 272 & tmp == 23 & replicate == 2),
                 !(seqID == 272 & tmp == 28 & replicate == 2),
                 !(seqID == 273 & tmp == 19),
                 !(seqID == 273 & tmp == 21 & replicate == 2),
                 !(seqID == 273 & tmp == 23 & replicate == 2),
                 !(seqID == 273 & tmp == 28 & replicate == 2),
                 !(seqID == 273 & tmp == 33 & replicate == 2),
                 !(seqID == 273 & tmp == 35 & replicate == 2),
                 !(seqID == 274 & tmp == 19),
                 !(seqID == 274 & tmp == 21),
                 !(seqID == 274 & tmp == 33 & replicate == 2),
                 !(seqID == 275 & tmp == 19),
                 !(seqID == 275 & tmp == 21),
                 !(seqID == 275 & tmp == 38 & replicate == 1),
                 !(seqID == 276 & tmp == 21),
                 !(seqID == 276 & tmp == 35 & replicate == 2),
                 !(seqID == 278 & tmp == 19),
                 !(seqID == 278 & tmp == 21 & replicate == 2),
                 !(seqID == 279 & tmp == 21),
                 !(seqID == 280 & tmp == 19),
                 !(seqID == 280 & tmp == 21),
                 !(seqID == 281 & tmp == 19),
                 !(seqID == 281 & tmp == 21),
                 !(seqID == 281 & tmp == 33 & replicate == 2),
                 !(seqID == 282 & tmp == 21),
                 !(seqID == 283 & tmp == 21),
                 !(seqID == 284 & tmp == 19),
                 !(seqID == 284 & tmp == 21),
                 !(seqID == 284 & tmp == 28 & replicate == 2),
                 !(seqID == 285 & tmp == 21 & replicate == 2),
                 !(seqID == 286 & tmp == 19),
                 !(seqID == 286 & tmp == 21),
                 !(seqID == 287 & tmp == 21),
                 !(seqID == 287 & tmp == 33 & replicate == 2),
                 !(seqID == 287 & tmp == 35 & replicate == 2),
                 !(seqID == 288 & tmp == 33 & replicate == 2),
                 !(seqID == 288 & tmp == 35 & replicate == 2),
                 !(seqID == 289 & tmp == 19),
                 !(seqID == 289 & tmp == 21),
                 !(seqID == 290 & tmp == 21),
                 !(seqID == 291 & tmp == 15 & replicate == 1),
                 !(seqID == 291 & tmp == 21),
                 !(seqID == 292 & tmp == 23 & replicate == 2),
                 !(seqID == 292 & tmp == 35),
                 !(seqID == 293 & tmp == 38 & replicate == 1),
                 !(seqID == 295 & tmp == 19),
                 !(seqID == 295 & tmp == 21),
                 !(seqID == 296 & tmp == 19),
                 !(seqID == 296 & tmp == 21),
                 !(seqID == 296 & tmp == 33 & replicate == 2),
                 !(seqID == 296 & tmp == 35 & replicate == 1),
                 !(seqID == 298 & tmp == 19),
                 !(seqID == 298 & tmp == 21),
                 !(seqID == 298 & tmp == 33 & replicate == 2),
                 !(seqID == 298 & tmp == 35 & replicate == 2),
                 !(seqID == 299 & tmp == 21),
                 !(seqID == 300 & tmp == 21),
                 !(seqID == 300 & tmp == 23 & replicate == 2),
                 !(seqID == 300 & tmp == 33),
                 !(seqID == 300 & tmp == 40 & replicate == 2),
                 !(seqID == 301 & tmp == 19),
                 !(seqID == 301 & tmp == 21),
                 !(seqID == 302 & tmp == 21),
                 !(seqID == 309 & tmp == 19 & replicate == 1 & well == "C09"),
                 !(seqID == 309 & tmp == 35 & replicate == 2),
                 !(seqID == 309 & tmp == 43 & replicate == 1),
                 !(seqID == 320 & tmp == 19 ),
                 !(seqID == 320 & tmp == 21 ),
                 !(seqID == 321 & tmp == 21),
                 !(seqID == 321 & tmp == 35),
                 !(seqID == 324 & tmp == 13 ),
                 !(seqID == 324 & tmp == 19 & replicate == 1 & well == "F07"),
                 !(seqID == 324 & tmp == 21 & replicate == 1 & well == "F07"),
                 !(seqID == 339 & tmp == 21 & replicate == 1 & well == "G05"),
                 !(seqID == 339 & tmp == 35 ),
                 !(seqID == 341 & tmp == 33 ),
                 !(seqID == 341 & tmp == 35 ),
                 !(seqID == 356 & tmp == 35 & replicate == 1),
                 !(seqID == 359 & tmp == 35 & replicate == 1),
                 !(seqID == 361 & tmp == 35 & replicate == 1),
                 !(seqID == 362 & tmp == 35 & replicate == 1),
                 !(seqID == 363 & tmp == 35 & replicate == 1),
                 !(seqID == 364 & tmp == 35 & replicate == 1),
                 !(seqID == 365 & tmp == 35 & replicate == 1), 
                 !(seqID == 374 & tmp == 35), 
                 !(seqID == 375 & tmp == 35) 
                 ) #

# We then rename a few variables for data compatibility with other databases
y_gcplyr_with_derivatives <- y_gcplyr_with_derivatives %>%
  dplyr::select(!c(temperature_of_isolation)) %>%
  rename("Site_initials" = "Reserve")%>%
  rename("temperature_of_isolation" = "Isolation.temperature")

write.csv(y_gcplyr_with_derivatives,"y_gcplyr_with_derivatives.csv",row.names = FALSE)


y_gcplyr_with_derivatives <- fread("y_gcplyr_with_derivatives.csv")



pdf("~/growthcurves_gcplyr_250324.pdf")
for(i in 1:120){
  print(ggplot(y_gcplyr_with_derivatives, aes(x = t, y = movavg_21,col=as.factor(replicate),group= interaction(seqID,well,tmp,replicate) )) +
          geom_line(lwd = 0.75)+
          facet_wrap_paginate(seqID~tmp, scales = 'free_y',nrow = 4,ncol=6, page = i)+
          theme(
            legend.position = 'none',
            panel.background = element_rect(fill='transparent'), #transparent panel bg
            plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank(), #remove minor gridlines
            legend.background = element_rect(fill='transparent'), #transparent legend bg
            legend.box.background = element_rect(fill='transparent'))+
          coord_cartesian(ylim = c(0, 0.75))
  )
}

dev.off()



#### Summarize parameters after filtering out very small absorbances ####
#
# 
# Here we calculate the growth parameters per temperature, ID, and replicate, 
# with different smoothing levels. We provide only a few examples for script conciseness
#
#

y_params_gcplyr_29_21 <- y_gcplyr_with_derivatives %>%
  group_by(seqID,tmp,replicate,Domain,Phylum,Class,Order,Family,Genus,Species,temperature_of_isolation,Area,Site_initials) %>%
  summarize(
    lag_time = lag_time(x = t, y = movavg_21, deriv = derivpercap_29_21),
    r = max(derivpercap_29_21[movavg_21 > 0.01 & t > 3], na.rm = TRUE),
    odmax = max(movavg_21[t > 3], na.rm = TRUE),
    auc = auc(y = movavg_21[t > 3], x = t[t > 3]),
    max_percap_time = extr_val(t, which_max_gc(derivpercap_29_21)),
    time_to_stationary = first_above(x = t, y = movavg_21, threshold = max(movavg_21[t > max_percap_time & derivpercap_29_21 == r*0.5 ], na.rm = TRUE),return = "x"),
    exponential_interval = (max_percap_time+time_to_stationary-lag_time)
  )

y_params_gcplyr_29_21$r[is.infinite(y_params_gcplyr_29_21$r)] <- 0
y_params_gcplyr_29_21$r[y_params_gcplyr_29_21$r<0] <- 0


y_params_gcplyr_9_21 <- y_gcplyr_with_derivatives %>%
  group_by(seqID,tmp,replicate,Domain,Phylum,Class,Order,Family,Genus,Species,temperature_of_isolation,Area,Site_initials) %>%
  summarize(
    lag_time = lag_time(x = t, y = movavg_21, deriv = derivpercap_9_21),
    r = max(derivpercap_9_21[movavg_21 > 0.01 & t > 3], na.rm = TRUE),
    odmax = max(movavg_21[t > 3], na.rm = TRUE),
    auc = auc(y = movavg_21[t > 3], x = t[t > 3]),
    max_percap_time = extr_val(t, which_max_gc(derivpercap_29_21)),
    time_to_stationary = first_above(x = t, y = movavg_21, threshold = max(movavg_21[t > max_percap_time & derivpercap_29_21 == r*0.5 ], na.rm = TRUE),return = "x"),
    exponential_interval = (max_percap_time+time_to_stationary-lag_time)
  )

y_params_gcplyr_9_21$r[is.infinite(y_params_gcplyr_9_21$r)] <- 0
y_params_gcplyr_9_21$r[y_params_gcplyr_9_21$r<0] <- 0



#### We merge the datasets again ####

y_params_gcplyr_29_21$derivative <- "29"
y_params_gcplyr_9_21$derivative <- "9"

y_params_gcplyr_29_21$moving_avg <- "21"
y_params_gcplyr_9_21$moving_avg <- "21"


y_params_gcplyr <- bind_rows(y_params_gcplyr_9_21,y_params_gcplyr_29_21)

y_params_gcplyr$derivative <- factor(y_params_gcplyr$derivative, levels = c("9","29"))


#### Save data for input for TPCs ####


write.csv(y_params_gcplyr, "240820_params_gcplyr.csv",row.names = FALSE)
write.csv(y_gcplyr_with_derivatives, "240820_fits_gcplyr.csv",row.names = FALSE)



