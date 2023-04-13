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
#                                     C3M31                                    #
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
  source(file.path(here::here(), "Arabidopsis", "C3M31", "Step#00_define_experiment.R"))

# Now the species and ID of the experiment can be found by entering:
#c(spcs, idExp)

# so that the directory of the experiment is now explicitly defined as:
dir_Exp <- file.path(here::here(), spcs, idExp)

# To check all constants, enter:
#set_constants_C3M31()

# To check the colors, enter:
#set_colors_C3M31()



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

# Keep only the time interval that will be used for fitting
df <- df[df$Time >= Time_start_exp, ]

# Make sure the fit line starts and ends at the same time as the data (useful here for graphical display)
ColorsPeriod$Time1[1] <- Time_start_exp
ColorsPeriod$Time2[nrow(ColorsPeriod)] <- Time_end_exp


#------------------------------------------------------------------------------#
#                     Run the fit and plot it for each pot                     #
#------------------------------------------------------------------------------#

# Run the fit and extract parameters (takes a few seconds)
results_fit <- run_fit(df,
                       use_VPD = F,
                       Time_var = "Time",
                       E_var = "E_corr",
                       Trt_var = c("idGenotype", "idWatering"),
                       col.per = ColorsPeriod,
                       linear_detrend = T,
                       period_for_detrend = c("Control 1", "Control 2"),
                       remove_mode = "global", max_skew = 2*60) # otherwise no fit is computed for "Control 2" that overlaps two consecutive cycles

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
                inside = F, irrig_show_mode = "pot", pot = pot)
      
    # Show the data
    points(E_corr ~ Time, data = dat, pch = 16, cex = 0.5, col = "grey") # includes the transition estimates, if any (also used for the fit)
    
    # Show the diel trend used for the fit
    diel_trend <- results_fit$diel_trend[results_fit$idPot == pot][1]
    diel_trend_baseline <- results_fit$diel_trend_baseline[results_fit$idPot == pot][1]
    abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
    
    # Add the fit and the data (here, not detrended) that have been used for the fit
    curve_fit(dat, use_VPD = F, Time_var = "Time", E_var = "E_corr", results = results_fit, col.per = ColorsPeriod, remove_mode = "global", max_skew = 2*60)
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
                               Trt_var = c("idGenotype", "idWatering"),
                               col.per = ColorsPeriod,
                               linear_detrend = T,
                               period_for_detrend = c("Control 1", "Control 2"),
                               remove_mode = "global", max_skew = 2*60) # otherwise no fit is computed for "Control 2" that overlaps two consecutive cycles

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
                inside = F, irrig_show_mode = "pot", pot = pot)
    
    # Show the data
    points(E_corr_per_kPa ~ Time, data = dat, pch = 16, cex = 0.5, col = "grey") # includes the transition estimates, if any (also used for the fit)
    
    # Show the diel trend used for the fit
    diel_trend <- results_fit_use_VPD$diel_trend[results_fit_use_VPD$idPot == pot][1]
    diel_trend_baseline <- results_fit_use_VPD$diel_trend_baseline[results_fit_use_VPD$idPot == pot][1]
    abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
    
    # Add the fit
    curve_fit(dat, use_VPD = T, Time_var = "Time", E_var = "E_corr_per_kPa", results = results_fit_use_VPD, col.per = ColorsPeriod, remove_mode = "global", max_skew = 2*60)
    }
  }  
dev.off()

# Save the results
write.csv(results_fit_use_VPD, file.path(dir_Exp, "Fitted_data", "results_fit_use_VPD.csv"), row.names = F)


#------------------------------------------------------------------------------#
# Same procedure with custom time intervals for comparison between experiments #
#------------------------------------------------------------------------------#

# Create the customized 'ColorsPeriod_custom' dataframe
ColorsPeriod_custom <- data.frame(idPeriod = c("Control", "Control_exclude 1", "High VPD", "Control_exclude 2"),
                                  Time1 = c(Time_start_exp, 1+Pho_Per/24/2, 2+3/24, 2+5/24),
                                  Time2 = c(1+Pho_Per/24/2, 2+3/24, 2+5/24, 3.5),
                                  col = c("black", rep(NA, 3)))

# Run the fit and extract parameters (takes a few seconds)
results_fit_use_VPD_custom <- run_fit(df,
                                      use_VPD = T,
                                      Time_var = "Time",
                                      E_var = "E_corr_per_kPa",
                                      Trt_var = c("idGenotype", "idWatering"),
                                      col.per = ColorsPeriod_custom,
                                      linear_detrend = T,
                                      period_for_detrend = c("Control", "Control_exclude 1", "Control_exclude 2"),
                                      remove_mode = "day_or_night") # better use the two night periods but only one day period (the first is incomplete, the third may be affected by the double irrigation)

# Keep relevant data
results_fit_use_VPD_custom <- results_fit_use_VPD_custom[results_fit_use_VPD_custom$idPeriod == "Control", ]
#ColorsPeriod_custom <- ColorsPeriod_custom[ColorsPeriod_custom$idPeriod == "Control", ]

# Plot the fits
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#04c_fit_all_pots_use_VPD_custom.pdf", sep = "_")), width = 8, height = 5)
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
    curve_fit(dat, use_VPD = T, Time_var = "Time", E_var = "E_corr_per_kPa", results = results_fit_use_VPD_custom, col.per = ColorsPeriod_custom[1, ], remove_mode = "day_or_night")
    }
  }  
dev.off()

# Save the results
write.csv(results_fit_use_VPD_custom, file.path(dir_Exp, "Fitted_data", "results_fit_use_VPD_custom.csv"), row.names = F)


