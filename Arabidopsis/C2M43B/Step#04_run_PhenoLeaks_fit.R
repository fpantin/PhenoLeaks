################################################################################
#                                                                              #
#              PhenoLeaks - Step #04 - Running the PhenoLeaks fit              #
#                                                                              #
#             Script to fit a sine model to the transpiration data             #
#          and extract parameters on the observed and fitted kinetics          #
#                                                                              #
#                             Florent Pantin, 2022                             #
#                                                                              #
################################################################################

#>                                                                            <#
#     >                                                                  <     #
#           >                                                      <           #
#                 >                                          <                 #
#                       >                              <                       #
#                             >                  <                             #
#                                    C2M43B                                    #
#                             >                  <                             #
#                       >                              <                       #
#                 >                                          <                 #
#           >                                                      <           #
#     >                                                                  <     #
#>                                                                            <#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                              (1)  Load libraries                             #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#options(repos = "https://cran.rstudio.com/") # RStudio
for (pkg in c("here"))
  {
  if (!pkg %in% installed.packages()[, "Package"]) { install.packages(pkg) }
  #update.packages(pkg)
  library(pkg, character.only = T)
  }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                             (2)  Source functions                            #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Note that paths are built relative to the root directory of the PhenoLeaks project.
dir_PhenoLeaks <- here::here("_core")
source(file.path(dir_PhenoLeaks, "PhenoLeaks_generic.R"))
source(file.path(dir_PhenoLeaks, "PhenoLeaks_graphics.R"))
source(file.path(dir_PhenoLeaks, "PhenoLeaks_fit.R"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                 (3)  Retrieve the features of the experiment                 #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Here the file "Step#00_define_experiment.R" should be sourced.

  ## OPTION 1: open the file and source it manually

  ## OPTION 2: set the correct file path and source it from this script, e.g.:
  source(file.path(here::here(), "Arabidopsis", "C2M43B", "Step#00_define_experiment.R"))

# Now the species and ID of the experiment can be found by entering:
#c(spcs, idExp)

# so that the directory of the experiment is now explicitly defined as:
dir_Exp <- file.path(here::here(), spcs, idExp)

# To check all constants, enter:
#set_constants_C2M43B()

# To check the colors, enter:
#set_colors_C2M43B()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                             (4)  Run the script                              #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



#------------------------------------------------------------------------------#
#                       Import and manage processed data                       #
#------------------------------------------------------------------------------#

# Import corrected transpiration data
df <- read.csv(file.path(dir_Exp, "Corrected_data", "E_corr.csv"))

# Make sure the fit line starts and ends at the same time as the data (useful here for graphical display)
ColorsPeriod$Time1[1] <- Time_start_exp
ColorsPeriod$Time2[nrow(ColorsPeriod)] <- Time_end_exp

# Assign the severity of SWC drop as in Step#03
minSWC <- aggregate(list(SWC = df$SWC), by = list(idPot = df$idPot), FUN = min, na.rm = T)
minSWC[order(minSWC$SWC), ]
df$SWC_drop <- "moderate"
fast <- minSWC$idPot[minSWC$SWC <= 0.8]
df$SWC_drop[df$idPot %in% fast] <- "fast"
slow <- minSWC$idPot[minSWC$SWC > 1.2]
df$SWC_drop[df$idPot %in% slow] <- "slow"


#------------------------------------------------------------------------------#
#                     Run the fit and plot it for each pot                     #
#------------------------------------------------------------------------------#

# Run the fit and extract parameters (takes a few seconds)
results_fit <- run_fit(df,
                       use_VPD = F,
                       Time_var = "Time",
                       E_var = "E_corr",
                       Trt_var = c("idGenotype", "idWatering", "SWC_drop"),
                       col.per = ColorsPeriod,
                       linear_detrend = T,
                       period_for_detrend = "Control") 

# Plot the fits
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#04a_fit_all_pots.pdf", sep = "_")), width = 8, height = 5)
for (geno in sort(unique(df$idGenotype)))
  {
  for (pot in unique(df$idPot[df$idGenotype == geno]))
    {
    # Prepare the plot
    dat <- df[df$idPot == pot, ]
    prepare_kin(dat,
                use_VPD = F,
                Time_var = "Time",
                E_var = "E_corr",
                main = paste(geno, sub(paste(idExp, "-", sep = ""), "Pot ", pot), sep = " - "),
                inside = F)
      
    # Show the data
    points(E_corr ~ Time, data = dat, pch = 16, cex = 0.5, col = "grey") # includes the transition estimates, if any (also used for the fit)

    # Show the diel trend used for the fit
    diel_trend <- results_fit$diel_trend[results_fit$idPot == pot][1]
    diel_trend_baseline <- results_fit$diel_trend_baseline[results_fit$idPot == pot][1]
    abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
    
    # Add the fit and the data (here, not detrended) that have been used for the fit
    curve_fit(dat, use_VPD = F, Time_var = "Time", E_var = "E_corr", results = results_fit, col.per = ColorsPeriod)
    }
  }  
dev.off()

# Save the results
if (!dir.exists(file.path(dir_Exp, "Fitted_data"))) { dir.create(file.path(dir_Exp, "Fitted_data")) }
write.csv(results_fit, file.path(dir_Exp, "Fitted_data", "results_fit.csv"), row.names = F)


#------------------------------------------------------------------------------#
#              Same procedure with transpiration normalized by VPD             #
#------------------------------------------------------------------------------#

# Run the fit and extract parameters (takes a few seconds)
results_fit_use_VPD <- run_fit(df,
                               use_VPD = T,
                               Time_var = "Time",
                               E_var = "E_corr_per_kPa",
                               Trt_var = c("idGenotype", "idWatering", "SWC_drop"),
                               col.per = ColorsPeriod,
                               linear_detrend = T,
                               period_for_detrend = "Control") 

# Plot the fits
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#04b_fit_all_pots_use_VPD.pdf", sep = "_")), width = 8, height = 5)
for (geno in sort(unique(df$idGenotype)))
  {
  for (pot in unique(df$idPot[df$idGenotype == geno]))
    {
    # Prepare the plot
    dat <- df[df$idPot == pot, ]
    prepare_kin(dat,
                use_VPD = T,
                Time_var = "Time",
                E_var = "E_corr_per_kPa",
                main = paste(geno, sub(paste(idExp, "-", sep = ""), "Pot ", pot), sep = " - "),
                inside = F)
    
    # Show the data
    points(E_corr_per_kPa ~ Time, data = dat, pch = 16, cex = 0.5, col = "grey") # includes the transition estimates, if any (also used for the fit)
    
    # Show the diel trend used for the fit
    diel_trend <- results_fit_use_VPD$diel_trend[results_fit_use_VPD$idPot == pot][1]
    diel_trend_baseline <- results_fit_use_VPD$diel_trend_baseline[results_fit_use_VPD$idPot == pot][1]
    abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
                                                                                          
    # Add the fit
    curve_fit(dat, use_VPD = T, Time_var = "Time", E_var = "E_corr_per_kPa", results = results_fit_use_VPD, col.per = ColorsPeriod)
    }
  }  
dev.off()

# Save the results
write.csv(results_fit_use_VPD, file.path(dir_Exp, "Fitted_data", "results_fit_use_VPD.csv"), row.names = F)


#------------------------------------------------------------------------------#
# Same procedure with custom time intervals for comparison between experiments #
#------------------------------------------------------------------------------#

# Create the customized 'ColorsPeriod_custom' dataframe
# The "Control" range is restricted to a short range because:
#    - "Control_exclude 1": the first day period is incomplete and introduces bias in the early morning
#    - "Control_exclude 2": the data after the second night are subjected to excessive soil drying for some pots
ColorsPeriod_custom <- data.frame(idPeriod = c("Control_exclude 1", "Control", "Control_exclude 2", "Low light"),
                                  Time1 = c(Time_start_exp, -0.5, 1, 3),
                                  Time2 = c(-0.5, 1, 3, 4),
                                  col = c(NA, "black", NA, "burlywood4"))

# The control range will be extended by one diel cycle for two pots, 116 and 120, which are not fitted otherwise due to missing data.
# This is acceptable because both pots have a slow soil drying.
ColorsPeriod_custom_116_120 <- data.frame(idPeriod = c("Control_exclude 1", "Control", "Control_exclude 2", "Low light"),
                                          Time1 = c(Time_start_exp, -0.5, 2, 3),
                                          Time2 = c(-0.5, 2, 3, 4),
                                          col = c(NA, "black", NA, "burlywood4"))


## 1 - Non-normalized data

# Run the fit and extract parameters (takes a few seconds)
results_fit_custom <- run_fit(df,
                              use_VPD = F,
                              Time_var = "Time",
                              E_var = "E_corr",
                              Trt_var = c("idGenotype", "idWatering", "SWC_drop"),
                              col.per = ColorsPeriod_custom,
                              linear_detrend = T,
                              period_for_detrend = c("Control_exclude 1", "Control", "Control_exclude 2"),
                              remove_mode = "day_or_night") # allows keeping several control observations that would be lost with the more stringent default ("day_and_night")

# Again for the two pots 116 and 120
results_fit_custom_116_120 <- run_fit(df[df$idPot %in% c("116", "120"), ],
                                      use_VPD = F,
                                      Time_var = "Time",
                                      E_var = "E_corr",
                                      Trt_var = c("idGenotype", "idWatering", "SWC_drop"),
                                      col.per = ColorsPeriod_custom_116_120,
                                      linear_detrend = T,
                                      period_for_detrend = c("Control_exclude 1", "Control", "Control_exclude 2"),
                                      remove_mode = "day_or_night") # allows keeping several control observations that would be lost with the more stringent default ("day_and_night")

# Merge the results in 'results_fit_custom'
results_fit_custom[results_fit_custom$idPot == "116", ] <- results_fit_custom_116_120[results_fit_custom_116_120$idPot == "116", ]
results_fit_custom[results_fit_custom$idPot == "120", ] <- results_fit_custom_116_120[results_fit_custom_116_120$idPot == "120", ]

# Prevent inclusion of EoN fitted data of the previous period which was excluded
results_fit_custom$E_EON_before_mod[results_fit_custom$idPeriod == "Low light"] <- NA

# Assign the correct 'A_rapid_op' to the control (which is here computed based on the previous excluded period, see 'E_val_obs()' in 'Phenoleaks_fit.R')
results_fit_custom$A_rapid_op_obs[results_fit_custom$idPeriod == "Control"] <- results_fit_custom$A_rapid_op_stable_obs[results_fit_custom$idPeriod == "Control"]

# Plot the fits
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#04c_fit_all_pots_custom.pdf", sep = "_")), width = 8, height = 5)
for (geno in sort(unique(df$idGenotype)))
  {
  for (pot in unique(df$idPot[df$idGenotype == geno]))
    {
    # Prepare the plot
    dat <- df[df$idPot == pot, ]
    prepare_kin(dat,
                use_VPD = F,
                Time_var = "Time",
                E_var = "E_corr",
                main = paste(geno, sub(paste(idExp, "-", sep = ""), "Pot ", pot), sep = " - "),
                inside = F)
    
    # Show all data
    points(E_corr ~ Time, data = dat, pch = 16, cex = 0.5, col = "grey") # includes the transition estimates, if any (also used for the fit)
    
    # Show the diel trend used for the fit
    diel_trend <- results_fit_custom$diel_trend[results_fit_custom$idPot == pot][1]
    diel_trend_baseline <- results_fit_custom$diel_trend_baseline[results_fit_custom$idPot == pot][1]
    abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
    
    # Overlay the points used for each period and add the fit
    curve_fit(dat, use_VPD = F, Time_var = "Time", E_var = "E_corr", results = results_fit_custom, col.per = if (!pot %in% c("116", "120")) ColorsPeriod_custom else ColorsPeriod_custom_116_120, remove_mode = "day_or_night")
    }
  }  
dev.off()

# Save the results
write.csv(results_fit_custom[results_fit_custom$idPeriod %in% c("Control", "Low light"), ],
          file.path(dir_Exp, "Fitted_data", "results_fit_custom.csv"), row.names = F)


## 2 - Data normalized by the VPD

# Run the fit and extract parameters (takes a few seconds)
results_fit_use_VPD_custom <- run_fit(df,
                                      use_VPD = T,
                                      Time_var = "Time",
                                      E_var = "E_corr_per_kPa",
                                      Trt_var = c("idGenotype", "idWatering", "SWC_drop"),
                                      col.per = ColorsPeriod_custom,
                                      linear_detrend = T,
                                      period_for_detrend = c("Control_exclude 1", "Control", "Control_exclude 2"),
                                      remove_mode = "day_or_night") # allows keeping several control observations that would be lost with the more stringent default ("day_and_night")

# Again for the two pots 116 and 120
results_fit_use_VPD_custom_116_120 <- run_fit(df[df$idPot %in% c("116", "120"), ],
                                              use_VPD = T,
                                              Time_var = "Time",
                                              E_var = "E_corr_per_kPa",
                                              Trt_var = c("idGenotype", "idWatering", "SWC_drop"),
                                              col.per = ColorsPeriod_custom_116_120,
                                              linear_detrend = T,
                                              period_for_detrend = c("Control_exclude 1", "Control", "Control_exclude 2"),
                                              remove_mode = "day_or_night") # allows keeping several control observations that would be lost with the more stringent default ("day_and_night")

# Merge the results in 'results_fit_use_VPD_custom'
results_fit_use_VPD_custom[results_fit_use_VPD_custom$idPot == "116", ] <- results_fit_use_VPD_custom_116_120[results_fit_use_VPD_custom_116_120$idPot == "116", ]
results_fit_use_VPD_custom[results_fit_use_VPD_custom$idPot == "120", ] <- results_fit_use_VPD_custom_116_120[results_fit_use_VPD_custom_116_120$idPot == "120", ]

# Prevent inclusion of EoN fitted data of the previous period which was excluded
results_fit_use_VPD_custom$E_EON_before_mod[results_fit_use_VPD_custom$idPeriod == "Low light"] <- NA

# Assign the correct 'A_rapid_op' to the control (which is here computed based on the previous excluded period, see 'E_val_obs()' in 'Phenoleaks_fit.R')
results_fit_use_VPD_custom$A_rapid_op_obs[results_fit_use_VPD_custom$idPeriod == "Control"] <- results_fit_use_VPD_custom$A_rapid_op_stable_obs[results_fit_use_VPD_custom$idPeriod == "Control"]

# Plot the fits
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#04d_fit_all_pots_use_VPD_custom.pdf", sep = "_")), width = 8, height = 5)
for (geno in sort(unique(df$idGenotype)))
  {
  for (pot in unique(df$idPot[df$idGenotype == geno]))
    {
    # Prepare the plot
    dat <- df[df$idPot == pot, ]
    prepare_kin(dat,
                use_VPD = T,
                Time_var = "Time",
                E_var = "E_corr_per_kPa",
                main = paste(geno, sub(paste(idExp, "-", sep = ""), "Pot ", pot), sep = " - "),
                inside = F)
    
    # Show all data
    points(E_corr_per_kPa ~ Time, data = dat, pch = 16, cex = 0.5, col = "grey") # includes the transition estimates, if any (also used for the fit)

    # Show the diel trend used for the fit
    diel_trend <- results_fit_use_VPD_custom$diel_trend[results_fit_use_VPD_custom$idPot == pot][1]
    diel_trend_baseline <- results_fit_use_VPD_custom$diel_trend_baseline[results_fit_use_VPD_custom$idPot == pot][1]
    abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
    
    # Overlay the points used for each period and add the fit
    curve_fit(dat, use_VPD = T, Time_var = "Time", E_var = "E_corr_per_kPa", results = results_fit_use_VPD_custom, col.per = if (!pot %in% c("116", "120")) ColorsPeriod_custom else ColorsPeriod_custom_116_120, remove_mode = "day_or_night")
    }
  }  
dev.off()

# Save the results
write.csv(results_fit_use_VPD_custom[results_fit_use_VPD_custom$idPeriod %in% c("Control", "Low light"), ],
          file.path(dir_Exp, "Fitted_data", "results_fit_use_VPD_custom.csv"), row.names = F)


#------------------------------------------------------------------------------#
#          Add a 24-h 'acclimation pulse' when the environment changes         #
#------------------------------------------------------------------------------#

# Set a new 'ColorsPeriod' table with acclimation ranges
ColorsPeriod_acclim_custom <- data.frame(idPeriod = c("Control_exclude 1", "Control", "Control_exclude 2", "Low light"),
                                  Time1 = c(Time_start_exp, -0.5, 1, 3),
                                  Time2 = c(-0.5, 1, 3, Time_end_exp),
                                  Time1_acclim = c(NA, NA, NA, 3),
                                  Time2_acclim = c(NA, NA, NA, 4),
                                  col = c(NA, "black", NA, "burlywood4"))


## 1 - Non-normalized data

# Run the fit and extract parameters (takes a few seconds) - will be useful only for low light
results_fit_acclim_custom_low_light <- run_fit(df,
                                               use_VPD = F,
                                               Time_var = "Time",
                                               E_var = "E_corr",
                                               Trt_var = c("idGenotype", "idWatering", "SWC_drop"),
                                               use_acclim = T,
                                               col.per = ColorsPeriod_acclim_custom,
                                               linear_detrend = T,
                                               period_for_detrend = c("Control_exclude 1", "Control", "Control_exclude 2"),
                                               remove_mode = "global", max_skew = 2*60) # necessary so that the start of the second day is not discarded

# Prevent inclusion of EoN fitted data of the previous period which was excluded
results_fit_acclim_custom_low_light$E_EON_before_mod[results_fit_acclim_custom_low_light$idPeriod == "Low light acclimating"] <- NA

# Keep relevant data and merge with the control
results_fit_acclim_custom <- rbind(results_fit_custom[results_fit_custom$idPeriod == "Control", ],
                                   results_fit_acclim_custom_low_light[results_fit_acclim_custom_low_light$idPeriod %in% c("Low light acclimating", "Low light acclimated"), ])

# Plot the fits
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#04f_fit_all_pots_acclim_custom.pdf", sep = "_")), width = 8, height = 5)
for (irr in c("WW", "WS"))
  {
  for (geno in sort(unique(df$idGenotype[df$idWatering == irr])))
    {
    for (pot in unique(df$idPot[df$idWatering == irr & df$idGenotype == geno]))
      {
      # Prepare the plot
      dat <- df[df$idPot == pot, ]
      prepare_kin(dat,
                  use_VPD = F,
                  Time_var = "Time",
                  E_var = "E_corr",
                  main = paste(geno, irr, sub(paste(idExp, "-", sep = ""), "Pot ", pot), sep = " - "),
                  inside = F)
      
      # Show the data
      points(E_corr ~ Time, data = dat, pch = 16, cex = 0.5, col = "grey") # includes the transition estimates, if any (also used for the fit)
      
      # Show the diel trend used for the fit
      diel_trend <- results_fit_acclim_custom$diel_trend[results_fit_acclim_custom$idPot == pot][1]
      diel_trend_baseline <- results_fit_acclim_custom$diel_trend_baseline[results_fit_acclim_custom$idPot == pot][1]
      abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
      
      # Overlay the points used for each period and add the fit
      curve_fit(dat, use_VPD = F, Time_var = "Time", E_var = "E_corr", results = results_fit_custom, col.per = if (!pot %in% c("116", "120")) ColorsPeriod_custom[1:2, ] else ColorsPeriod_custom_116_120[1:2, ], remove_mode = "day_or_night") # note that idPeriod 1 ("Control_exclude 1", color = NA) is kept here just to make sure the display of idPeriod 2 ("Control") starts at the right time, i.e. not including t1
      curve_fit(dat, use_VPD = F, Time_var = "Time", E_var = "E_corr", results = results_fit_acclim_custom_low_light, use_acclim = T, col.per = ColorsPeriod_acclim_custom[3:4, ], remove_mode = "global", max_skew = 2*60) # note that idPeriod 3 ("Control_exclude 2", color = NA) is kept here just to make sure the display of idPeriod 4 ("Low light") starts at the right time, i.e. not including t1
      }
    }  
  }
dev.off()

# Save the results
write.csv(results_fit_acclim_custom, file.path(dir_Exp, "Fitted_data", "results_fit_acclim_custom.csv"), row.names = F)


## 2 - Data normalized by the VPD

# Run the fit and extract parameters (takes a few seconds) - will be useful only for low light
results_fit_acclim_use_VPD_custom_low_light <- run_fit(df,
                                                       use_VPD = T,
                                                       Time_var = "Time",
                                                       E_var = "E_corr_per_kPa",
                                                       Trt_var = c("idGenotype", "idWatering", "SWC_drop"),
                                                       use_acclim = T,
                                                       col.per = ColorsPeriod_acclim_custom,
                                                       linear_detrend = T,
                                                       period_for_detrend = c("Control_exclude 1", "Control", "Control_exclude 2"),
                                                       remove_mode = "global", max_skew = 2*60) # necessary so that the start of the second day is not discarded

# Prevent inclusion of EoN fitted data of the previous period which was excluded
results_fit_acclim_use_VPD_custom_low_light$E_EON_before_mod[results_fit_acclim_use_VPD_custom_low_light$idPeriod == "Low light acclimating"] <- NA

# Keep relevant data and merge with the control
results_fit_acclim_use_VPD_custom <- rbind(results_fit_use_VPD_custom[results_fit_use_VPD_custom$idPeriod == "Control", ],
                                           results_fit_acclim_use_VPD_custom_low_light[results_fit_acclim_use_VPD_custom_low_light$idPeriod %in% c("Low light acclimating", "Low light acclimated"), ])

# Plot the fits
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#04f_fit_all_pots_acclim_use_VPD_custom.pdf", sep = "_")), width = 8, height = 5)
for (irr in c("WW", "WS"))
  {
  for (geno in sort(unique(df$idGenotype[df$idWatering == irr])))
    {
    for (pot in unique(df$idPot[df$idWatering == irr & df$idGenotype == geno]))
      {
      # Prepare the plot
      dat <- df[df$idPot == pot, ]
      prepare_kin(dat,
                  use_VPD = T,
                  Time_var = "Time",
                  E_var = "E_corr_per_kPa",
                  main = paste(geno, irr, sub(paste(idExp, "-", sep = ""), "Pot ", pot), sep = " - "),
                  inside = F)
      
      # Show the data
      points(E_corr_per_kPa ~ Time, data = dat, pch = 16, cex = 0.5, col = "grey") # includes the transition estimates, if any (also used for the fit)
      
      # Show the diel trend used for the fit
      diel_trend <- results_fit_acclim_use_VPD_custom$diel_trend[results_fit_acclim_use_VPD_custom$idPot == pot][1]
      diel_trend_baseline <- results_fit_acclim_use_VPD_custom$diel_trend_baseline[results_fit_acclim_use_VPD_custom$idPot == pot][1]
      abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
      
      # Overlay the points used for each period and add the fit
      curve_fit(dat, use_VPD = T, Time_var = "Time", E_var = "E_corr_per_kPa", results = results_fit_use_VPD_custom, col.per = if (!pot %in% c("116", "120")) ColorsPeriod_custom[1:2, ] else ColorsPeriod_custom_116_120[1:2, ], remove_mode = "day_or_night") # note that idPeriod 1 ("Control_exclude 1", color = NA) is kept here just to make sure the display of idPeriod 2 ("Control") starts at the right time, i.e. not including t1
      curve_fit(dat, use_VPD = T, Time_var = "Time", E_var = "E_corr_per_kPa", results = results_fit_acclim_use_VPD_custom_low_light, use_acclim = T, col.per = ColorsPeriod_acclim_custom[3:4, ], remove_mode = "global", max_skew = 2*60) # note that idPeriod 3 ("Control_exclude 2", color = NA) is kept here just to make sure the display of idPeriod 4 ("Low light") starts at the right time, i.e. not including t1
      }
    }  
  }
dev.off()

# Save the results
write.csv(results_fit_acclim_use_VPD_custom, file.path(dir_Exp, "Fitted_data", "results_fit_acclim_use_VPD_custom.csv"), row.names = F)

