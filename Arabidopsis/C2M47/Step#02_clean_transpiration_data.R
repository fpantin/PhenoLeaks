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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                              (1)  Load libraries                             #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#options(repos = "https://cran.rstudio.com/") # RStudio
for (pkg in c("scales"))
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


dir_PhenoLeaks <- file.path(dirname(getwd()), "PhenoLeaks")
source(file.path(dir_PhenoLeaks, "PhenoLeaks_generic.R"))
source(file.path(dir_PhenoLeaks, "PhenoLeaks_outliers.R"))
source(file.path(dir_PhenoLeaks, "PhenoLeaks_graphics.R"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                 (3)  Retrieve the features of the experiment                 #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


source(file.path(getwd(), "Step#00_define_experiment.R"))

# To check the constants, run:
#set_constants_C2M47()

# To check the colors, enter:
#set_colors_C2M47()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                             (4)  Run the script                              #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#------------------------------------------------------------------------------#
#                       Import and manage processed data                       #
#------------------------------------------------------------------------------#


# Transpiration data
# calculated from "cleaned_outliers.csv" and "Transpiration_calculations_v2.Rmd"
df <- read.csv(file.path(getwd(), "Processed_data", "C2M47_starch_pot_transpiration_SWC.csv"))
df <- df[, c("idPot", "idGenotype", "idWatering", "decimalDay", "E_mmol_per_m2_s", "meanVPD", "E_mmol_per_m2_s_kPa", "SWC")]

# Format the dataframe with a uniform time sequence for all pots (if not done at Step#01, add lines and print warning messages)
df <- format_time(df, Time_var = "decimalDay", Trt_var = c("idGenotype", "idWatering"), time_step = 30)

# Create the column 'Trt' (concatenation of the statistical treatments)
df$Trt <- paste(df$idGenotype, df$idWatering, sep = " - ")



#------------------------------------------------------------------------------#
#                       Plot kinetics of individual pots                       #
#------------------------------------------------------------------------------#

pdf(file.path(getwd(), "Figures", paste(idExp, "Step#02a_plot_all_pots.pdf", sep = "_")), width = 10, height = 5)
for (irr in c("WW", "WS"))
  {
  for (geno in sort(unique(df$idGenotype[df$idWatering == irr])))
    {
    for (pot in sort(unique(df$idPot[df$idWatering == irr & df$idGenotype == geno])))
      {
      dat <- df[df$idPot == pot, ]
      prepare_kin(dat, use_VPD = T, Time_var = "Time", E_var = "E_mmol_per_m2_s_kPa",
                  main = paste(geno, irr, paste("Pot", pot, sep = " "), sep = " - "),
                  inside = F, irrig_show_mode = "pot", pot = pot,
                  add_SWC = T, mar = c(2.5, 3.5, 2.5, 7.5))
      points(E_mmol_per_m2_s_kPa ~ Time, data = dat, type = "o", cex = 0.5)
      rm(dat)
      }
    }
  }

dev.off()


pdf(file.path(getwd(), "Figures", paste(idExp, "Step#02b_plot_all_pots_per_condition.pdf", sep = "_")), width = 10, height = 5)
n.max <- max(aggregate(df$idPot[!duplicated(df$idPot)], by = list(df$Trt[!duplicated(df$idPot)]), FUN = length)$x) # maximum number of replicates per treatment
for (irr in c("WW", "WS"))
  {
  for (geno in sort(unique(df$idGenotype[df$idWatering == irr])))
    {
    df.geno <- df[df$idWatering == irr & df$idGenotype == geno, ]
    prepare_kin(df.geno, Time_var = "Time", E_var = "E_mmol_per_m2_s",
                main = paste(geno, irr, sep = " - "),
                inside = F)
    color = 0
    for (pot in sort(unique(df.geno$idPot[df.geno$idGenotype == geno])))
      {
      #require(scales)
      color <- color + 1
      dat <- df[df$idPot == pot, ]
      points(E_mmol_per_m2_s ~ Time, data = dat, type = "o", col = hue_pal()(n.max)[color], cex = 0.5)
      rm(dat)
      }
    legend("top", ncol = 4, bty = "n",
           legend = paste("Pot", sort(unique(df.geno$idPot[df.geno$idGenotype == geno])), sep = " "),
           col = hue_pal()(n.max)[1:color], lty = 1, pch = 21, pt.cex = 0.5)
    rm(df.geno)
    }
  }

dev.off()



#------------------------------------------------------------------------------#
#                     Identify series of outliers visually                     #
#------------------------------------------------------------------------------#

# Inspection of the PDF files generated above allows identifying outliers visually,
# based on  series of inconsistent data (that will be identified as grey points in the next plots).

  ## Create the column 'outlier' for identification (TRUE or FALSE)
  df$outlier <- F

  
  ## Initial set-up - All pots
  
  # Data of the first 1.25 days (before 'Time_start_exp') are systematically discarded
  # because there were noisy (due to mechanical stability of the mobile scale that was then improved)
  df$outlier[df$Time < Time_start_exp] <- T
  
  
  ## Well-watered ("WW")
  
  # Col-0
  df$outlier[df$idPot == "295" & df$Time > 3 & df$Exact_Time_min <= 3*60*24 + 4*60] <- T
  df$outlier[df$idPot == "295" & df$Time > 5] <- T
  df$outlier[df$idPot == "350" & df$Time > 5] <- T
  
  # Col-0 ProStarv::LUC
  df$outlier[df$idPot == "222" & df$Time > Pho_Per/24 & df$Time <= 1] <- T
  
  # pgm
  df$outlier[df$idPot == "250" & df$Exact_Time_min > 4*60*24 + Pho_Per*60 & df$Exact_Time_min <= 4*60*24 + Pho_Per*60 + 3.5*60] <- T
  
  # isa1
  df$outlier[df$idPot == "320" & df$Exact_Time_min >= 3*60*24 - 1.5*60 & df$Exact_Time_min <= 3*60*24 + 3.5*60] <- T
  ##############df$outlier[df$idPot == "320" & df$Time > 3 & df$Time <= 4] <- T
  df$outlier[df$idPot == "320" & df$Time > 5] <- T
  df$outlier[df$idPot == "358" & df$Time > 5] <- T
  
  # ss4
  df$outlier[df$idPot == "221" & df$Time > 2 & df$Time <= 3] <- T
  df$outlier[df$idPot == "221" & df$Time > 5 + Pho_Per/24] <- T
  df$outlier[df$idPot == "376" & df$Time > 5 + Pho_Per/24] <- T
  
  # sex1
  df$outlier[df$idPot == "290" & df$Time > 5 + Pho_Per/24] <- T
  
  # mex1
  df$outlier[df$idPot == "219" & df$Time > 5 + Pho_Per/24] <- T
  df$outlier[df$idPot == "269" & df$Time > 4 & df$Time <= 5] <- T
  df$outlier[df$idPot == "279" & df$Time > 3 & df$Time <= 4] <- T
  
  # dpe1
  df$outlier[df$idPot == "244" & df$Time > 5] <- T
  df$outlier[df$idPot == "262" & df$Time > 2 & df$Time <= 3] <- T
  df$outlier[df$idPot == "352" & df$Time > 5] <- T
  
  # amy3
  df$outlier[df$idPot == "302" & df$Time > 5] <- T
  df$outlier[df$idPot == "334" & df$Time > 2 & df$Time <= 3 + 1.25/24] <- T
  ##################################df$outlier[df$idPot == "334" & df$Time > 2 & df$Time <= 3] <- T
  
  # amy3 bam1
  df$outlier[df$idPot == "232" & df$Time < - 3/24] <- T
  df$outlier[df$idPot == "248" & df$Time > 4 & df$Time <= 5] <- T
  df$outlier[df$idPot == "325" & df$Time > 5] <- T
  
  # bam3
  df$outlier[df$idPot == "370" & df$Time > 3 & df$Time <= 4] <- T
  
  # abcb14-2
  df$outlier[df$idPot == "294" & df$Time > 1.5 & df$Time <= 3] <- T
  
  
  ## Water stress ("WS")
  
  # Col-0 WS
  df$outlier[df$idPot == "224" & df$Time < - Sko_Per/24] <- T
  df$outlier[df$idPot == "365" & df$Time > - Sko_Per/24 & df$Time <= 0.25] <- T
  df$outlier[df$idPot == "374" & df$Time > Pho_Per/2/24 & df$Exact_Time_min <= Pho_Per*60 + 30] <- T
  df$outlier[df$idPot == "374" & df$Time >= 1 & df$Time <= 1 + Pho_Per/2/24] <- T

  # pgm WS
  df$outlier[df$idPot == "263" & df$Time > 4 & df$Time <= 5] <- T
  df$outlier[df$idPot == "326" & df$Time > 5 + Pho_Per/24] <- T


# Check the percent of points that are visually considered as outliers (beyond 'Time_start_exp')
print(paste(round(length(which(df$outlier[df$Time >= Time_start_exp])) / nrow(df[df$Time >= Time_start_exp, ]) * 100, 1), "%", sep = ""))
  
# It is also possible to keep every data and just run the automatic outlier detection below.
# It can then be easily verified that the most of the data of the periods visually identified as outliers above
# are also automatically detected as suspicious.
# However, once verified, it is strongly advised to come back to this manual removal of strong outliers
# because removing these strong outliers will make the function 'detect_outliers()' more sensitive for the detection
# of more subtle outliers.

  

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
                      E_var = "E_clean",
                      Trt_var = c("idGenotype", "idWatering"),
                      lags = 1:4)

# Check the range of 'sum_out' (higher value = more chance for the point to be an outlier)
hist(df$sum_out)

# Visualize individual curves with outliers
pdf(file.path(getwd(), "Figures", paste(idExp, "Step#02c_visualize_potential_outliers_all_pots.pdf", sep = "_")), width = 10, height = 5)
for (irr in c("WW", "WS"))
  {
  for (geno in sort(unique(df$idGenotype[df$idWatering == irr])))
    {
    for (pot in sort(unique(df$idPot[df$idWatering == irr & df$idGenotype == geno])))
      {
      # Prepare the plot
      dat <- df[df$idPot == pot, ]
      prepare_kin(dat, use_VPD = F, Time_var = "Time", E_var = "E_mmol_per_m2_s",
                  main = paste(geno, irr, paste("Pot", pot, sep = " "), sep = " - "),
                  inside = F, irrig_show_mode = "pot", pot = pot)
    
      # Print data
      points(E_clean ~ Time, data = dat, pch = 16, cex = 0.5)
      # which - at this point - is the same as:
      #points(E_mmol_per_m2_s ~ Time, data = dat[!dat$outlier, ] pch = 16, cex = 0.5)
    
      # Print outliers identified manually (whole periods)
      points(E_mmol_per_m2_s ~ Time, data = dat[dat$outlier, ], pch = 1, cex = 0.5, col = "grey")
    
      # Visualize different thresholds for the sum of penalties
      points(E_mmol_per_m2_s ~ Time, data = dat[!dat$outlier & dat$sum_out >= 25, ], pch = 16, col = 4, cex = 0.5)
      points(E_mmol_per_m2_s ~ Time, data = dat[!dat$outlier & dat$sum_out >= 35, ], pch = 16, col = 3, cex = 0.5)
      points(E_mmol_per_m2_s ~ Time, data = dat[!dat$outlier & dat$sum_out >= 45, ], pch = 16, col = 2, cex = 0.5)
    
      # Add the legend
      legend("topright", inset = c(0.1, 0), pch = 1, col = "grey", cex = 0.5, bty = "n",
             legend = "Outlier", title = "VISUALLY")
      legend("topright", inset = c(0, 0), pch = 16, col = 4:2, cex = 0.5, bty = "n",
             legend = c("25+", "35+", "45+"), title = "STD.RESID")
    
      rm(dat)
      }
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
                        E_var = "E_clean",
                        time_start = Time_start_exp,
                        time_end = Time_end_exp,
                        interval_est = 0.15)

# Recompute the data normalized by VPD
df$E_estimate_per_kPa <- df$E_estimate / df$meanVPD


# Check that the estimation is correct based on non-missing data
#require(ggplot2)
#require(ggpmisc)
  
pdf(file.path(getwd(), "Figures", paste(idExp, "Step#02d_check_transition_estimates.pdf", sep = "_")), width = 8, height = 5)

maxi <- max(c(df$E_clean[!is.na(df$E_estimate)], df$E_estimate), na.rm = T)

ggplot(df, aes(E_clean, E_estimate, color = Trt)) +
  ggtitle("All data") + xlab("Observed E") + ylab("Estimated E") +
  geom_point() + xlim(-1, maxi) + ylim(-1, maxi) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
ggplot(df[df$after_dark2light,], aes(E_clean, E_estimate, color = Trt)) +
  ggtitle("Start of day") + xlab("Observed E") + ylab("Estimated E") +
  geom_point() + xlim(-1, maxi) + ylim(-1, maxi) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
ggplot(df[df$light2dark,], aes(E_clean, E_estimate, color = Trt)) +
  ggtitle("End of day") + xlab("Observed E") + ylab("Estimated E") +
  geom_point() + xlim(-1, maxi) + ylim(-1, maxi) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
ggplot(df[df$after_light2dark,], aes(E_clean, E_estimate, color = Trt)) +
  ggtitle("Start of night") + xlab("Observed E") + ylab("Estimated E") +
  geom_point() + xlim(-1, maxi) + ylim(-1, maxi) +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")))
ggplot(df[df$dark2light, ], aes(E_clean, E_estimate, color = Trt)) +
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
pdf(file.path(getwd(), "Figures", paste(idExp, "Step#02e_visualize_corrected_data_all_pots.pdf", sep = "_")), width = 10, height = 5)
for (irr in c("WW", "WS"))
  {
  for (geno in sort(unique(df$idGenotype[df$idWatering == irr])))
    {
    for (pot in sort(unique(df$idPot[df$idWatering == irr & df$idGenotype == geno])))
      {
      # Prepare the plot
      dat <- df[df$idPot == pot, ]
      prepare_kin(dat,
                  Time_var = "Time",
                  E_var = "E_mmol_per_m2_s",
                  main = paste(geno, irr, paste("Pot", pot, sep = " "), sep = " - "),
                  inside = F, irrig_show_mode = "pot", pot = pot)
      
      # Print robust data
      points(E_clean ~ Time, data = dat, pch = 16, cex = 0.5)
      
      # Print outliers identified manually (whole periods)
      points(E_mmol_per_m2_s ~ Time, data = dat[dat$outlier, ], pch = 1, cex = 0.5, col = "grey")
        
      # Print outliers identified automatically
      points(E_mmol_per_m2_s ~ Time, data = dat[!dat$outlier & dat$rule_out, ], pch = 16, col = "red", cex = 0.5)
      
      # Print transition estimates
      points(E_estimate ~ Time, data = dat[is.na(dat$E_clean), ], pch = "*", cex = 1, col = "goldenrod2")
      
      # Add the legend
      legend("topright", ncol = 3, pch = c(1, 16, charToRaw("*")), col = c("grey", "red", "goldenrod2"), cex = 0.5, pt.cex = c(0.5, 0.5, 1), bty = "n",
             legend = c("Outlier (visual)", "Outlier (auto)", "Transition estimate"))
      
      rm(dat)
      }
    }
  }

dev.off()
    

# # Reset time in the column 'Time' so that 0 corresponds to the desired dark-to-light transition (here Time_ON0)
# df$Time <- df$decimalDay - round(Time_ON0, n.digits)

# Create a column 'E_corr' that corresponds to 'E_clean' (where outliers are NA)
# but updated for the missing transition estimates ('E_estimate')
df$E_corr <- df$E_clean
df$E_corr[is.na(df$E_corr)] <- df$E_estimate[is.na(df$E_corr)]
df$E_corr_per_kPa <- df$E_clean_per_kPa
df$E_corr_per_kPa[is.na(df$E_corr_per_kPa)] <- df$E_estimate_per_kPa[is.na(df$E_corr_per_kPa)]

# Save the results in 'E_corr.csv'
df.save <- df[, c("Trt", "idGenotype", "idWatering", "idPot", "Time", "Exact_Time_min", "SWC", "E_mmol_per_m2_s", "E_clean", "E_estimate", "E_corr", "meanVPD", "E_mmol_per_m2_s_kPa", "E_clean_per_kPa", "E_estimate_per_kPa", "E_corr_per_kPa")]
write.csv(df.save, file.path(getwd(), "Corrected_data", "E_corr.csv"), row.names = F)

