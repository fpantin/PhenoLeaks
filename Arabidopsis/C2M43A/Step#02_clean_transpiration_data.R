################################################################################
#                                                                              #
#            PhenoLeaks - Step #02 - Cleaning the transpiration data           #
#                                                                              #
#               Script to clean the processed transpiration data               #
#                        through outlier identification                        #
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
#                                    C2M43A                                    #
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
for (pkg in c("scales", "here"))
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
source(file.path(dir_PhenoLeaks, "PhenoLeaks_outliers.R"))
source(file.path(dir_PhenoLeaks, "PhenoLeaks_graphics.R"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                 (3)  Retrieve the features of the experiment                 #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# Here the file "Step#00_define_experiment.R" should be sourced.

  ## OPTION 1: open the file and source it manually

  ## OPTION 2: set the correct file path and source it from this script, e.g.:
  source(file.path(here::here(), "Arabidopsis", "C2M43A", "Step#00_define_experiment.R"))

# Now the species and ID of the experiment can be found by entering:
#c(spcs, idExp)

# so that the directory of the experiment is now explicitly defined as:
dir_Exp <- file.path(here::here(), spcs, idExp)

# To check all constants, enter:
#set_constants_C2M43A()

# To check the colors, enter:
#set_colors_C2M43A()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                             (4)  Run the script                              #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#------------------------------------------------------------------------------#
#                       Import and manage processed data                       #
#------------------------------------------------------------------------------#


# Transpiration data calculated from Step#01
df <- read.csv(file.path(dir_Exp, "Processed_data", "C2M43A_pot_transpiration.csv"))
df <- df[, c("idObs", "idPot", "idGenotype", "idWatering", "Time", "Exact_Time_min", "E_mmol_per_m2_s", "meanVPD", "E_mmol_per_m2_s_kPa", "SWC")]



#------------------------------------------------------------------------------#
#                     Identify series of outliers visually                     #
#------------------------------------------------------------------------------#

# Inspection of the PDF files generated at Step#01 (transpiration) allows identifying outliers visually,
# based on  series of inconsistent data (that will be identified as grey points in the next plots).

  ## Create the column 'outlier' for identification (TRUE or FALSE)
  df$outlier <- F

  ## All pots
  #none
  
  ## Col-0
  df <- df[!df$idPot %in% c("1"), ]


# Check the percent of points that are visually considered as outliers (beyond 'Time_start_exp' and not considering fully removed pots)
print(paste(round(length(which(df$outlier[df$Time >= Time_start_exp])) / nrow(df[df$Time >= Time_start_exp, ]) * 100, 1), "%", sep = ""))
  
# Rename lines and keep order in the 'idObs' column.
# This is required because lines have been deleted for plants on which all data were discarded (pot 22).
row.names(df) <- 1:nrow(df)
df$idObs <- 1:nrow(df)
  
  
#------------------------------------------------------------------------------#
#                        Identify outliers automatically                       #
#------------------------------------------------------------------------------#

# Create the variable 'E_clean' and 'E_clean_kPa' with NA instead of outlier values identified manually,
# to make subsequent outlier detection more sensitive
df$E_clean <- df$E_mmol_per_m2_s
df$E_clean[df$outlier] <- NA
df$E_clean_per_kPa <- df$E_mmol_per_m2_s_kPa
df$E_clean_per_kPa[df$outlier] <- NA

# Assign "penalties" to each point depending on the local transpiration variation
# compared to that of the other plants of the same treatment(s) at the same time point
# in the new column 'sum_out' (takes less than 1 min)
df <- detect_outliers(df,
                      E_var = "E_clean_per_kPa",
                      Trt_var = "idGenotype",
                      lags = 1:4)

# Check the range of 'sum_out' (higher value = more chance for the point to be an outlier)
hist(df$sum_out)

# Visualize individual curves with outliers
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#02a_visualize_potential_outliers_all_pots.pdf", sep = "_")), width = 10, height = 5)
for (geno in sort(unique(df$idGenotype)))
  {
  for (pot in sort(unique(df$idPot[df$idGenotype == geno])))
    {
    # Prepare the plot
    dat <- df[df$idPot == pot, ]
    prepare_kin(dat, use_VPD = T, Time_var = "Time", E_var = "E_mmol_per_m2_s_kPa",
                main = paste(geno, sub(paste(idExp, "-", sep = ""), "Pot ", pot), sep = " - "),
                inside = F, irrig_show_mode = "pot", pot = pot)
    
    # Print data
    points(E_clean_per_kPa ~ Time, data = dat, pch = 16, cex = 0.5)
    # which - at this point - is the same as:
    #points(E_mmol_per_m2_s_kPa ~ Time, data = dat[!dat$outlier, ] pch = 16, cex = 0.5)
    
    # Print outliers identified manually (whole periods)
    points(E_mmol_per_m2_s_kPa ~ Time, data = dat[dat$outlier, ], pch = 1, cex = 0.5, col = "grey")
    
    # Visualize different thresholds for the sum of penalties
    points(E_mmol_per_m2_s_kPa ~ Time, data = dat[!dat$outlier & dat$sum_out >= 25, ], pch = 16, col = 4, cex = 0.5)
    points(E_mmol_per_m2_s_kPa ~ Time, data = dat[!dat$outlier & dat$sum_out >= 35, ], pch = 16, col = 3, cex = 0.5)
    points(E_mmol_per_m2_s_kPa ~ Time, data = dat[!dat$outlier & dat$sum_out >= 45, ], pch = 16, col = 2, cex = 0.5)
    
    # Add the legend
    legend("topright", inset = c(0.1, 0), pch = 1, col = "grey", cex = 0.5, bty = "n",
           legend = "Outlier", title = "VISUALLY")
    legend("topright", inset = c(0, 0), pch = 16, col = 4:2, cex = 0.5, bty = "n",
           legend = c("25+", "35+", "45+"), title = "STD.RESID")

    rm(dat)
    }
  }

dev.off()


# Look at the kinetics and define the rule for setting points as outliers
df$rule_out <- df$sum_out >= 35

# Check the percent of points that are considered as outliers using this threshold
print(paste(round(length(which(df$rule_out[!df$outlier])) / length(which(!df$outlier)) * 100, 1), "%", sep = ""))

# Assign NA for identified outliers
df$E_clean[df$rule_out] <- NA
df$E_clean_per_kPa[df$rule_out] <- NA



#------------------------------------------------------------------------------#
#              Re-estimate transpiration values at the transitions             #
#------------------------------------------------------------------------------#

# Transpiration values at the transitions and just after are required
# for parameter estimation by the function 'E_val_obs()' of the script 'PhenoLeaks_fit.R'.
# Re-estimation is done for each pot and each period, since non-missing transition values will be used
# to check that the prediction is correct.

# Run the function (takes a few seconds)
df <- estim_transitions(df,
                        Time_var = "Time",
                        E_var = "E_clean_per_kPa",
                        time_start = Time_start_exp,
                        time_end = Time_end_exp,
                        interval_est = 0.15)

df$E_estimate_per_kPa <- df$E_estimate

# Recompute the data non-normalized by VPD
#df$E_estimate <- df$E_estimate * df$meanVPD # does not work where missing VPD values so we have to run 'estim_transitions()' using "E_clean"
df <- df[names(df) != "E_estimate"]
df <- estim_transitions(df,
                        Time_var = "Time",
                        E_var = "E_clean",
                        time_start = Time_start_exp,
                        time_end = Time_end_exp,
                        interval_est = 0.15)


# Check that the estimation is correct based on non-missing data
#require(ggplot2)
#require(ggpmisc)
  
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#02b_check_transition_estimates.pdf", sep = "_")), width = 8, height = 5)

maxi <- max(c(df$E_clean_per_kPa[!is.na(df$E_estimate_per_kPa)], df$E_estimate_per_kPa), na.rm = T)

ggplot(df, aes(E_clean_per_kPa, E_estimate_per_kPa, color = idGenotype)) +
  ggtitle("All data") + xlab("Observed E") + ylab("Estimated E") +
  geom_point() + xlim(-1, maxi) + ylim(-1, maxi) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
ggplot(df[df$after_dark2light,], aes(E_clean_per_kPa, E_estimate_per_kPa, color = idGenotype)) +
  ggtitle("Start of day") + xlab("Observed E") + ylab("Estimated E") +
  geom_point() + xlim(-1, maxi) + ylim(-1, maxi) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
ggplot(df[df$light2dark,], aes(E_clean_per_kPa, E_estimate_per_kPa, color = idGenotype)) +
  ggtitle("End of day") + xlab("Observed E") + ylab("Estimated E") +
  geom_point() + xlim(-1, maxi) + ylim(-1, maxi) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
ggplot(df[df$after_light2dark,], aes(E_clean_per_kPa, E_estimate_per_kPa, color = idGenotype)) +
  ggtitle("Start of night") + xlab("Observed E") + ylab("Estimated E") +
  geom_point() + xlim(-1, maxi) + ylim(-1, maxi) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
ggplot(df[df$dark2light, ], aes(E_clean_per_kPa, E_estimate_per_kPa, color = idGenotype)) +
  ggtitle("End of night") + xlab("Observed E") + ylab("Estimated E") +
  geom_point() + xlim(-1, maxi) + ylim(-1, maxi) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
  
dev.off()


#------------------------------------------------------------------------------#
#                       Plot the corrected data and save                       #
#------------------------------------------------------------------------------#
  
# Visualize individual curves with outliers and transition estimates
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#02c_visualize_corrected_data_all_pots.pdf", sep = "_")), width = 10, height = 5)
for (geno in sort(unique(df$idGenotype)))
  {
  for (pot in sort(unique(df$idPot[df$idGenotype == geno])))
    {
    # Prepare the plot
    dat <- df[df$idPot == pot, ]
    prepare_kin(dat,
                use_VPD = T,
                Time_var = "Time",
                E_var = "E_mmol_per_m2_s_kPa",
                ylim = c(0, 4),
                main = paste(geno, sub(paste(idExp, "-", sep = ""), "Pot ", pot), sep = " - "),
                inside = F, irrig_show_mode = "pot", pot = pot)
      
    # Print robust data
    points(E_clean_per_kPa ~ Time, data = dat, pch = 16, cex = 0.5)
      
    # Print outliers identified manually (whole periods)
    points(E_mmol_per_m2_s_kPa ~ Time, data = dat[dat$outlier, ], pch = 1, cex = 0.5, col = "grey")
        
    # Print outliers identified automatically
    points(E_mmol_per_m2_s_kPa ~ Time, data = dat[!dat$outlier & dat$rule_out, ], pch = 16, col = "red", cex = 0.5)
      
    # Print transition estimates
    points(E_estimate_per_kPa ~ Time, data = dat[is.na(dat$E_clean_per_kPa), ], pch = "*", cex = 1, col = "goldenrod2")
      
    # Add the legend
    legend("topright", ncol = 3, pch = c(1, 16, charToRaw("*")), col = c("grey", "red", "goldenrod2"), cex = 0.5, pt.cex = c(0.5, 0.5, 1), bty = "n",
           legend = c("Outlier (visual)", "Outlier (auto)", "Transition estimate"))
      
    rm(dat)
    }
  }

dev.off()
    

# Create a column 'E_corr' that corresponds to 'E_clean' (where outliers are NA)
# but updated for the missing transition estimates ('E_estimate')
df$E_corr <- df$E_clean
df$E_corr[is.na(df$E_corr)] <- df$E_estimate[is.na(df$E_corr)]
df$E_corr_per_kPa <- df$E_clean_per_kPa
df$E_corr_per_kPa[is.na(df$E_corr_per_kPa)] <- df$E_estimate_per_kPa[is.na(df$E_corr_per_kPa)]

# Save the results in 'E_corr.csv'
df$Trt <- paste(df$idGenotype, df$idWatering, sep = " - ") # just for easier use of the next script
df.save <- df[, c("Trt", "idGenotype", "idWatering", "idPot", "Time", "Exact_Time_min", "SWC", "E_mmol_per_m2_s", "E_clean", "E_estimate", "E_corr", "meanVPD", "E_mmol_per_m2_s_kPa", "E_clean_per_kPa", "E_estimate_per_kPa", "E_corr_per_kPa")]
if (!dir.exists(file.path(dir_Exp, "Corrected_data"))) { dir.create(file.path(dir_Exp, "Corrected_data")) }
write.csv(df.save, file.path(dir_Exp, "Corrected_data", "E_corr.csv"), row.names = F)

