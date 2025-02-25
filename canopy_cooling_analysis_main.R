# Introduction ----------------------------

#Purpose: Analyze temperature anomaly data for air or LST using GAMMs

#Created: MGA 9/10/2020
#Updated: MGA 11/2024
#
#Inputs: Point shapefiles for each city containing temperature data for air
# (anomaly) and LST (LST_avg). Also, a suite of predictor variables summarized
# within 90 m buffers

#Output: numerical results, plots, and tables
#
#Notes:  It is helpful to use RStudio "document outline" functionality.
#TCF means "tree canopy fraction" which is synonymous with "%canopy"
#The primary workflow for replicability is as follows:
# 1. Run from 00 through 03 to load data
# 2. in 04 make choices about subsetting data by time of day
# 3. Go into 6.1, 6.2, 6.3 to see which GAMM you want to run
# 4. Run the loop once or more (depends on time of day goals)
# 5. Save data for plotting or other downstream uses


# 00 libraries ---------------------------- 

#possible not all of these are necessary....

library(mgcv)
library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
library(forcats)
library(sf)
library(gamm4)
library(doSNOW)
library(rgdal)
library(rCAT)
library(ape)
library(matrixStats)
library(spdep)

# 01 USER: Set initial parameters  -------------------------------------------------------

#This is the first place to tinker with the main parameters for running this code
#But for subsetting data and model choice there are more choices to be made below

#set working directory (where all input data such as shapefiles should live)
wd <- "F:/OneDriveAU/Projects2/Urban_Heat/share/" #full filepath string

#directory for saving analysis plots and other variables (as needed)
#will be created automatically
plot_dir <- "./analysis/r-plots/"


# 02 Load city data -------------------------------------------------------

#set working directory
setwd(wd)

#Baltimore
dat_orig <- st_read("Baltimore_results_v4.shp")
dat_balt <- dat_orig
dat_balt <- dat_balt %>% mutate("city"="balt")

#Boston
dat_orig <- st_read("Boston_results_v4.shp")
dat_bos <- dat_orig
dat_bos <- dat_bos %>% mutate("city"="bos")

#Burlington
dat_orig <- st_read("Burlington_results_v4.shp")
dat_burl <- dat_orig
dat_burl <- dat_burl %>% mutate("city"="burl")

#Cincinnati
dat_orig <- st_read("Cincinnati_result_v4.shp")
dat_cinc <- dat_orig
dat_cinc <- dat_cinc %>% mutate("city"="cinc")

#Durham
dat_orig <- st_read("Durham_result_v4.shp")
dat_durh <- dat_orig
dat_durh <- dat_durh %>% mutate("city"="durh")
#set watrdis_90 to 13000 for the whole city. That is mostly just a placeholder
#since water wont be relevant but it is the distance to the nearest, large 
#lake in the NE
dat_durh$watrdis_90 <- rep(13000, nrow(dat_durh))

#Detroit
dat_orig <- st_read("Detroit_result_v4.shp")
#first correct detroit's elev column
dat_orig <- dat_orig %>% rename(elev_90=elev1m_90)
dat_det <- dat_orig
dat_det <- dat_det %>% mutate("city"="det")

#Newark
dat_orig <- st_read("Newark_results_v4.shp")
dat_new <- dat_orig
dat_new <- dat_new %>% mutate("city"="new")

#NYC
dat_orig <- st_read("NewYork_result_v4.shp")
dat_nyc <- dat_orig
dat_nyc <- dat_nyc %>% mutate("city"="nyc")

#Richmond
dat_orig <- st_read("Richmond_results_v4.shp")
dat_rich <- dat_orig
dat_rich <- dat_rich %>% mutate("city"="rich")

#Virgina Beach
dat_orig <- st_read("virginiabeach_results_v4.shp")
dat_vb <- dat_orig
dat_vb <- dat_vb %>% mutate("city"="vb")

#Wash_DC
dat_orig <- st_read("dc_full_dataset_match_col_sel_v4.shp")
dat_dc <- dat_orig
#ditch the extra station columns (can swap out if needed)
dat_dc <- dat_dc %>% select(-c("west_temp", "wapo_temp"))


# 03 Combine city data --------------------------

#list of all the cities
city_list <- list(dat_balt, dat_bos, dat_burl, dat_cinc, dat_det, dat_durh, dat_new, dat_nyc, dat_rich, dat_vb, dat_dc)

#get all cities on same crs and rbind
for (i in 1:length(city_list)) {
  
  this_city <- city_list[[i]]
  this_city_wgs84 <- st_transform(this_city, 4326, CRS(this_city))
  
  if (i==1) {
    dat_orig <- this_city_wgs84
  } else {
    dat_orig <- bind_rows(dat_orig, this_city_wgs84)
  }
  print(i)
}

#set city name to factor 
dat_orig$city <- as.factor(dat_orig$city)

#make a locally-normalized elevation variable. Using mean centering only. 
#We want to take the regional trend away but allow elevation ranges to be larger
#or smaller by city
dat_orig <- dat_orig %>% group_by(city) %>% mutate(elev_norm=(elev_90 - mean(elev_90, na.rm=T))) %>% ungroup()

#change the fractions to percentages
dat_orig$tcf_90 <- dat_orig$tcf_90*100
dat_orig$imp_90 <- dat_orig$imp_90*100


# 04 USER: Choose subset of data to analyze  ----------------------------------------

#User slices and dices the data and models here

#AM AF PM data entry or subset criteria
num_fixed_effects <- 5 #usually 5 for full model or 1 for canopy only

#are we subsetting the data just to reduce data loads?
#for the main analysis should put "1" which means "use all data"
data_reduction_factor <- 1 #"5" means "take 1 in every 5 points systematically"

#which time period are you modeling?
#am = morning, af = afternoon, pm = evening
user_time_period <- "af" #"am", "af", "pm"

#other data exclusions 
speed_low <- 3 #take out very low speeds (default = 3)
speed_high <- 40 #take out very high speeds (default = 40). Speeds should all be mph.


# 05 Prep data for analysis loop  -------------------------------------------------------

#user rarely alters this code

#re-assemble a true date_time column from the as.character version
dat_orig$date_tm <- as.POSIXct(dat_orig$date_tm, "%Y-%m-%d %H:%M:%OS", tz="EST5EDT")
#make sure that wind_dr is factor
dat_orig$wind_dr <- as.factor(dat_orig$wind_dr)
#make the geometry into lon and lat variables
coords_df <- st_coordinates(dat_orig)
dat_orig <- dat_orig %>% mutate(lon=coords_df[,1], lat=coords_df[,2])
#Examine data types for the dataset
str(dat_orig)


# * 5.1 Further data subsetting  -------------------------

#option to subset data systematically
#but if data_reduction_factor is 1, this doesn't do anything
#user rarely alters this code

#reduce this dataset SYSTEMATICALLY by factor of data_reduction_factor
seq_idx <- seq(1,nrow(dat_orig), by=data_reduction_factor)
dat <- dat_orig[seq_idx,]
coords_df <- coords_df[seq_idx,]
#generic subsetting for various reasons (e.g., time period to model, car speed, all data inside study area)
sub_idx <- which(dat$trav_id==user_time_period & dat$speed>speed_low & dat$speed<speed_high & 
                  dat$tcf_90 < 99 & dat$anomaly > -100 & dat$lcbin_90>0 & dat$elev_90>-5 &
                   dat$watrdis_90>=0)
dat <- dat[sub_idx,]
coords_df <- coords_df[sub_idx,]
#get rid of points at exact same location
not_dupes_idx <- which(duplicated(dat[,c("lat","lon")])==FALSE)
dat <- dat[not_dupes_idx,]
coords_df <- coords_df[not_dupes_idx,]

#is your current data frame the expected size?
nrow(dat)

# 06 Analysis loop  -----------------------------------------

#the user can make some choices here regarding subsampling by city
#or, mostly commonly, which model to run in the loop (e.g., full model,
#canopy-only model, which heat metric as response)

#single city??
#dat <- dat %>% filter(city=="vb")

#in advance of within-city random sampling, make city sample
#sizes the same
city_n <- dat %>% group_by(city) %>% summarise(n=n())
sample_city <- min(city_n$n)
sample_city
#now subset the whole dataset
dat <- dat %>%
  group_by(city) %>%
  sample_n(sample_city) %>%
  filter(wdistGSWD>=0000)

#set up for parallel processing
cl <- makeCluster(16,outfile="") # Make clusters
registerDoSNOW(cl)
#set up for a doSNOW progress bar
iterations <- 50
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
start_time <- Sys.time()
#Initialize the holder for generic model outputs
fits_df <- data.frame("var1_x"=rep(NA, 100), "var1_fit"=rep(NA,100), "var1_fit_se"=rep(NA,100),
                      "var2_x"=rep(NA, 100), "var2_fit"=rep(NA,100), "var2_fit_se"=rep(NA,100),
                      "var3_x"=rep(NA, 100), "var3_fit"=rep(NA,100), "var3_fit_se"=rep(NA,100),
                      "var4_x"=rep(NA, 100), "var4_fit"=rep(NA,100), "var4_fit_se"=rep(NA,100))
                      

#the parallel loop
results_list <- foreach (i = 1:iterations, .options.snow = opts) %dopar% { 
  
  library(dplyr)
  library(tidyr)
  library(mgcv)
  library(spdep)
  library(gstat)
  
  #get random subsample
  #dat_sub <- dat %>% sample_frac(0.1)
  
  #systematic sample to maintain constant distance between points
  # seq_idx <- seq(i,(nrow(dat)-(iterations-i)), by=iterations)
  # dat_sub <- dat[seq_idx,]
  
  #for all cities, strat random by city
  # dat_sub <- dat %>%
  #   group_by(city) %>%
  #   filter(watrdis_90>0) %>%
  #   sample_n(500)
  
  #non stratified sampling to allow more inclusion/exclusion of a given city
  #since the main variability is cross-city
  #THIS IS THE MAIN WAY THAT DATA WERE SUBSET IN EACH ITERATION
  dat_sub <- dat %>% slice_sample(prop = 0.02)
  
  
  # *6.1 multivariate model (with city as random effect)------------------------------
  
  #this is the main, multivariate model
  #user can switch between response variables of "anomaly" for air
  #and "LST_avg" for LST
  gam_mod_90 <-  gam(data = dat_sub, LST_avg ~ s(tcf_90, k=3) +
                       s(imp_90, k=3) + s(wdistGSWD, k=3) + s(elev_norm,k=3) +
                       s(lon,lat,k=5) + s(city, bs="re"),
                              method = "REML",
                              select = TRUE)
  
  # *6.2 canopy-only model (with city as random effect)------------------------------
  
  #this is the main, univariate model
  #user can switch between response variables of "anomaly" for air
  #and "LST_avg" for LST
  # gam_mod_90 <-  gam(data = dat_sub, LST_avg ~ s(tcf_90, k=3) + s(city, bs="re"),
  #                    method = "REML",
  #                    select = TRUE)
  
  # *6.3 single city models -----------------------------------------------------
  
  #single city model
  # gam_mod_90 <-  gam(data = dat_sub, anomaly ~ s(patch_90, k=3) + s(dstrb_90, k=3) +
  #                      s(imp_90, k=3) +
  #                      s(wdistGSWD, k=3) + s(elev_norm,k=3) +
  #                      ti(tcf_90, imp_90, k=5),
  #                    method = "REML",
  #                    select = TRUE)
  
  #single city model replacing elev and water with s(lon,lat)
  # gam_mod_90 <-  gam(data = dat_sub, anomaly ~ s(tcf_90, k=3) + s(imp_90, k=3) +
  #                     ti(tcf_90, imp_90, k=5) + s(lon,lat,k=27),
  #                    method = "REML",
  #                    select = TRUE)

  
  #assemble model outputs
  gam_mod <- gam_mod_90
  gam_info <- plot.gam(gam_mod, seWithMean = TRUE, scheme = 3, n=100, n2 = 100)
  
  # # #get mean fits for each term (TCF model)
  #right now this is just set up to deal with FULL model or TCF-only model
  if (num_fixed_effects>=5) {
    fits_df$var1_x<- gam_info[[1]]$x
    fits_df$var1_fit <- gam_info[[1]]$fit
    fits_df$var1_fit_se <- gam_info[[1]]$se
    fits_df$var2_x<- gam_info[[2]]$x
    fits_df$var2_fit <- gam_info[[2]]$fit
    fits_df$var2_fit_se <- gam_info[[2]]$se
    fits_df$var3_x<- gam_info[[3]]$x
    fits_df$var3_fit <- gam_info[[3]]$fit
    fits_df$var3_fit_se <- gam_info[[3]]$se
    fits_df$var4_x<- gam_info[[4]]$x
    fits_df$var4_fit <- gam_info[[4]]$fit
    fits_df$var4_fit_se <- gam_info[[4]]$se
  } else if (num_fixed_effects==3) { # tcf, imp, interaction
    fits_df$var1_x<- gam_info[[1]]$x
    fits_df$var1_fit <- gam_info[[1]]$fit
    fits_df$var1_fit_se <- gam_info[[1]]$se
    fits_df$var2_x<- gam_info[[2]]$x
    fits_df$var2_fit <- gam_info[[2]]$fit
    fits_df$var2_fit_se <- gam_info[[2]]$se
  } else { 
    fits_df$var1_x<- gam_info[[1]]$x
    fits_df$var1_fit <- gam_info[[1]]$fit
    fits_df$var1_fit_se <- gam_info[[1]]$se
  }
  
  #retrieve items from summary object
  gam_mod_summary <- summary(gam_mod)
  r2 <- gam_mod_summary$r.sq
  gam_mod_summary_table <- gam_mod_summary$s.table
  
 #assemble results from this iteration
  results_this_iter <- list(fits_df, gam_mod_summary_table, r2)
  
  return(results_this_iter)
  
}

#clean up doSNOW stuff, shutting down cluster
close(pb)
stopCluster(cl)
registerDoSEQ()

#makes a list that is 7*iterations long to separate out each sub list
results_list_flat <- unlist(results_list, recursive = FALSE)
fit_results_list <- results_list_flat[seq(1,length(results_list_flat)-3,by=3)]
smooth_terms_table_results_list <- results_list_flat[seq(2,length(results_list_flat)-1,by=3)]
r2_results_list <- results_list_flat[seq(3,length(results_list_flat),by=3)]

#final results
smooth_terms_table <- Reduce(`+`, smooth_terms_table_results_list) / length(smooth_terms_table_results_list)
smooth_terms_table
fit_df_mean <- Reduce(`+`, fit_results_list) / length(fit_results_list)

#print out some basic stats from this run
r2_vals <- unlist(r2_results_list)
hist(r2_vals); r2_mean <- mean(r2_vals); r2_mean


# 07 Save loop data by time of day ---------------------------------------------

#User can save results by time of day here. This requires manual uncommenting
#of the time of day in question.

#After saving data from each time of day, you can then run the main GAMM output plot

# *7.1 Morning ("am") ------------------------------------------------
# #
  # var1_fit_df_mean_90_pd <- fit_df_mean
  # tcf_r2_vals_90_pd <- r2_vals
  # tcf_smooth_terms_table_90_pd <- smooth_terms_table
 

# *7.2 Afternoon ("af") -----------------------------------------------------------
# # 
 # var1_fit_df_mean_90 <- fit_df_mean
 # tcf_r2_vals_90 <- r2_vals
 # tcf_smooth_terms_table_90 <- smooth_terms_table

# *7.3 Evening ("pm") --------------------------------------------------------------
# 
# tcf_smooth_terms_table_90_eve <- smooth_terms_table
# var1_fit_df_mean_90_eve <- fit_df_mean
# tcf_r2_vals_90_eve <- r2_vals


#report r2 values and smooth curve params
#am
mean(tcf_r2_vals_90_pd)
tcf_smooth_terms_table_90_pd
#mean(tcf_moran_90_pd)
#af
mean(tcf_r2_vals_90)
tcf_smooth_terms_table_90
#mean(tcf_moran_90)
#pm
mean(tcf_r2_vals_90_eve)
tcf_smooth_terms_table_90_eve
#mean(tcf_moran_90_eve)


#Can save off an R object with all times of day data
am_af_pm_AIR_model_20240917 <- list(var1_fit_df_mean_90_pd, tcf_r2_vals_90_pd,tcf_smooth_terms_table_90_pd,
                                        var1_fit_df_mean_90, tcf_r2_vals_90, tcf_smooth_terms_table_90,
                                    var1_fit_df_mean_90_eve, tcf_smooth_terms_table_90_eve, tcf_r2_vals_90_eve
                                        )


# 08 Make one nice plot ----------------------------


#TCF plots but only 3 panesl total by time of day
#change these depending on LST_avg or anomaly
fit_max <- 3
fit_min <- -10

#plotting tcf with averaged mean and quantiles of the n=50 or whatever fits
tcf_plot <- ggplot() +
  geom_ribbon(data=var1_fit_df_mean_90_pd, aes(x=var1_x, y=var1_fit, ymin = var1_fit - var1_fit_se, ymax = var1_fit + var1_fit_se), fill = "blue", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90_pd, aes(x=var1_x, y = var1_fit), color="blue") +
  geom_ribbon(data=var1_fit_df_mean_90_eve, aes(x=var1_x, y=var1_fit,  ymin = var1_fit - var1_fit_se, ymax = var1_fit + var1_fit_se), fill = "green", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90_eve, aes(x=var1_x, y = var1_fit), color="forestgreen") +
  geom_ribbon(data=var1_fit_df_mean_90, aes(x=var1_x, y=var1_fit,  ymin = var1_fit - var1_fit_se, ymax = var1_fit + var1_fit_se), fill = "black", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90, aes(x=var1_x, y = var1_fit), color="black") + 
  ylim(c(fit_min, fit_max)) + xlab("% cover") + ylab("anomaly (C)") +
  theme_light() + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =1))

imp_plot <- ggplot() +
  geom_ribbon(data=var1_fit_df_mean_90_pd, aes(x=var2_x, y=var2_fit, ymin = var2_fit - var2_fit_se, ymax = var2_fit + var2_fit_se), fill = "blue", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90_pd, aes(x=var2_x, y = var2_fit), color="blue") +
  geom_ribbon(data=var1_fit_df_mean_90_eve, aes(x=var2_x, y=var2_fit, ymin = var2_fit - var2_fit_se, ymax = var2_fit + var2_fit_se), fill = "green", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90_eve, aes(x=var2_x, y = var2_fit), color="forestgreen") +
  geom_ribbon(data=var1_fit_df_mean_90, aes(x=var2_x, y=var2_fit, ymin = var2_fit - var2_fit_se, ymax = var2_fit + var2_fit_se), fill = "black", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90, aes(x=var2_x, y = var2_fit), color="black") + 
  ylim(c(fit_min, fit_max)) + xlab("% cover") + ylab("") +
  theme_light() + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =1))

wtr_plot <- ggplot() +
  geom_ribbon(data=var1_fit_df_mean_90_pd, aes(x=var3_x, y=var3_fit, ymin = var3_fit - var3_fit_se, ymax = var3_fit + var3_fit_se), fill = "blue", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90_pd, aes(x=var3_x, y = var3_fit), color="blue") +
  geom_ribbon(data=var1_fit_df_mean_90_eve, aes(x=var3_x, y=var3_fit,ymin = var3_fit - var3_fit_se, ymax = var3_fit + var3_fit_se), fill = "green", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90_eve, aes(x=var3_x, y = var3_fit), color="forestgreen") +
  geom_ribbon(data=var1_fit_df_mean_90, aes(x=var3_x, y=var3_fit, ymin = var3_fit - var3_fit_se, ymax = var3_fit + var3_fit_se), fill = "black", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90, aes(x=var3_x, y = var3_fit), color="black") + 
  ylim(c(fit_min, fit_max)) + xlab("meters") + ylab("") +
  theme_light() + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =1))

elev_plot <- ggplot() +
  geom_ribbon(data=var1_fit_df_mean_90_pd, aes(x=var4_x, y=var4_fit, ymin = var4_fit - var4_fit_se, ymax = var4_fit + var4_fit_se), fill = "blue", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90_pd, aes(x=var4_x, y = var4_fit), color="blue") +
  geom_ribbon(data=var1_fit_df_mean_90_eve, aes(x=var4_x, y=var4_fit,ymin = var4_fit - var4_fit_se, ymax = var4_fit + var4_fit_se), fill = "green", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90_eve, aes(x=var4_x, y = var4_fit), color="forestgreen") +
  geom_ribbon(data=var1_fit_df_mean_90, aes(x=var4_x, y=var4_fit, ymin = var4_fit - var4_fit_se, ymax = var4_fit + var4_fit_se), fill = "black", alpha=0.3) +
  geom_line(data=var1_fit_df_mean_90, aes(x=var4_x, y = var4_fit), color="black") + 
  ylim(c(fit_min, fit_max)) + xlab("meters") + ylab("") +
  theme_light() + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =1))


#multivariable plot
p <- plot_grid(
  tcf_plot + ggtitle("a) Tree canopy fraction") + theme(plot.title = element_text(hjust = -0.3, vjust = 2.0)),
  imp_plot + ggtitle("b) Impervious") + theme(plot.title = element_text(hjust = -0.3, vjust = 2.0)),
  wtr_plot + ggtitle("c) Water Distance") + theme(plot.title = element_text(hjust = -0.3, vjust = 2.0)),
  elev_plot + ggtitle("d) Local elevation") + theme(plot.title = element_text(hjust = -0.3, vjust = 2.0)),
  label_fontface = "bold",
  label_colour = "black",
  nrow = 2
)
dev.new(width =6, height = 6 , noRStudioGD = TRUE); p



