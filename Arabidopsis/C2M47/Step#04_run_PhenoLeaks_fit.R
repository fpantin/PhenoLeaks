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
  source(file.path(here::here(), "Arabidopsis", "C2M47", "Step#00_define_experiment.R"))

# Now the species and ID of the experiment can be found by entering:
#c(spcs, idExp)

# so that the directory of the experiment is now explicitly defined as:
dir_Exp <- file.path(here::here(), spcs, idExp)

# To check all constants, enter:
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
                       period_for_detrend = c("Control", "Recovery 1", "Recovery 2"),
                       remove_mode = "day_or_night") # allows including the first night period while discarding the first day period (incomplete)

# Plot the fits
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#04a_fit_all_pots.pdf", sep = "_")), width = 8, height = 5)
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
      diel_trend <- results_fit$diel_trend[results_fit$idPot == pot][1]
      diel_trend_baseline <- results_fit$diel_trend_baseline[results_fit$idPot == pot][1]
      abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
      
      # Add the fit and the data (here, not detrended) that have been used for the fit
      curve_fit(dat, use_VPD = F, Time_var = "Time", E_var = "E_corr", results = results_fit, col.per = ColorsPeriod, remove_mode = "day_or_night")
      }
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
                               period_for_detrend = c("Control", "Recovery 1", "Recovery 2"),
                               remove_mode = "day_or_night") # allows including the first night period while discarding the first day period (incomplete)

# Plot the fits
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#04b_fit_all_pots_use_VPD.pdf", sep = "_")), width = 8, height = 5)
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
      diel_trend <- results_fit_use_VPD$diel_trend[results_fit_use_VPD$idPot == pot][1]
      diel_trend_baseline <- results_fit_use_VPD$diel_trend_baseline[results_fit_use_VPD$idPot == pot][1]
      abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
      
      # Add the fit
      curve_fit(dat, use_VPD = T, Time_var = "Time", E_var = "E_corr_per_kPa", results = results_fit_use_VPD, col.per = ColorsPeriod, remove_mode = "day_or_night")
      }
    }  
  }
dev.off()

# Save the results
write.csv(results_fit_use_VPD, file.path(dir_Exp, "Fitted_data", "results_fit_use_VPD.csv"), row.names = F)


#------------------------------------------------------------------------------#
# Add a 24-h 'acclimation pulse' when the environment changes (non-normalized) #
#------------------------------------------------------------------------------#

# Set a new 'ColorsPeriod' table with acclimation ranges
ColorsPeriod_acclim <- data.frame(idPeriod = c("Control", "High CO2", "Low light", "Recovery"),
                                  Time1 = c(Time_start_exp, 2, 3, 4),
                                  Time2 = c( 2, 3, 4, Time_end_exp),
                                  Time1_acclim = c(NA, 2, 3, 4),
                                  Time2_acclim = c(NA, 3, 4, 5),
                                  col = c("black", "coral3", "burlywood4", "gray30"))

# Run the fit and extract parameters (takes a few seconds)
results_fit_acclim <- run_fit(df,
                              use_VPD = F,
                              Time_var = "Time",
                              E_var = "E_corr",
                              Trt_var = c("idGenotype", "idWatering"),
                              use_acclim = T,
                              col.per = ColorsPeriod_acclim,
                              linear_detrend = T,
                              period_for_detrend = c("Control", "Recovery"),
                              remove_mode = "day_or_night") # allows including the first night period while discarding the first day period (incomplete)

# Plot the fits
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#04c_fit_all_pots_acclim.pdf", sep = "_")), width = 8, height = 5)
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
      diel_trend <- results_fit_acclim$diel_trend[results_fit_acclim$idPot == pot][1]
      diel_trend_baseline <- results_fit_acclim$diel_trend_baseline[results_fit_acclim$idPot == pot][1]
      abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
      
      # Add the fit and the data (here, not detrended) that have been used for the fit
      curve_fit(dat, use_VPD = F, Time_var = "Time", E_var = "E_corr", results = results_fit_acclim, use_acclim = T, col.per = ColorsPeriod_acclim, remove_mode = "day_or_night")
      }
    }  
  }
dev.off()

# Save the results
write.csv(results_fit_acclim, file.path(dir_Exp, "Fitted_data", "results_fit_acclim.csv"), row.names = F)


#------------------------------------------------------------------------------#
# Add a 24-h 'acclimation pulse' when the environment changes (VPD-normalized) #
#------------------------------------------------------------------------------#

# Run the fit and extract parameters (takes a few seconds)
results_fit_acclim_use_VPD <- run_fit(df,
                                      use_VPD = T,
                                      Time_var = "Time",
                                      E_var = "E_corr_per_kPa",
                                      Trt_var = c("idGenotype", "idWatering"),
                                      use_acclim = T,
                                      col.per = ColorsPeriod_acclim,
                                      linear_detrend = T,
                                      period_for_detrend = c("Control", "Recovery"),
                                      remove_mode = "day_or_night") # allows including the first night period while discarding the first day period (incomplete)

# Plot the fits
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#04d_fit_all_pots_acclim_use_VPD.pdf", sep = "_")), width = 8, height = 5)
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
      diel_trend <- results_fit_acclim_use_VPD$diel_trend[results_fit_acclim_use_VPD$idPot == pot][1]
      diel_trend_baseline <- results_fit_acclim_use_VPD$diel_trend_baseline[results_fit_acclim_use_VPD$idPot == pot][1]
      abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
      
      # Add the fit and the data (here, not detrended) that have been used for the fit
      curve_fit(dat, use_VPD = T, Time_var = "Time", E_var = "E_corr_per_kPa", results = results_fit_acclim_use_VPD, use_acclim = T, col.per = ColorsPeriod_acclim, remove_mode = "day_or_night")
      }
    }  
  }
dev.off()

# Save the results
write.csv(results_fit_acclim_use_VPD, file.path(dir_Exp, "Fitted_data", "results_fit_acclim_use_VPD.csv"), row.names = F)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                           (4)  Demo for publication                          #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Retrieve the original 'ColorsPeriod' (otherwise the black rectangle is incomplete)
source(file.path(dir_Exp, "Step#00_define_experiment.R"))

# Set the function
plot_individual_fit <- function (pot,
                                 YLIM = c(0.35, 2),
                                 X_LEGEND = 2.25,
                                 Y_LEGEND = 1.95,
                                 LEGEND_TITLE = pot,
                                 draw_x_axis = T,
                                 draw_rect_env = T)
                                 #res = results_fit,
                                 #use_acclim = any(!is.na(res$acclim_slope)),
                                 #remove_mode = "day_or_night",
                                 #col.per = ColorsPeriod)
  {
  # Prepare the plot
  dat <- df[df$idPot == pot, ]
  prepare_kin(dat,
              use_VPD = F,
              Time_var = "Time",
              E_var = "E_corr",
              ylim = YLIM,
              xlim = c(-0.7, 5.9),
              draw_x_axis = draw_x_axis,
              draw_rect_env = draw_rect_env,
              export_PPTX = T)
  
  # Show the data
  points(E_corr ~ Time, data = dat, pch = 16, cex = 0.5, col = "grey") # includes the transition estimates, if any (also used for the fit)
  
  # # Show the diel trend used for the fit
  # diel_trend <- res$diel_trend[res$idPot == pot][1]
  # diel_trend_baseline <- res$diel_trend_baseline[res$idPot == pot][1]
  # abline(a = diel_trend_baseline, b = diel_trend, lty = 3, col = "grey")
  
  # Add the fit and the data (here, not detrended) that have been used for the fit
  #curve_fit(dat, use_VPD = F, Time_var = "Time", E_var = "E_corr", results = res, use_acclim = use_acclim, col.per = col.per, export_PPTX = T, remove_mode = remove_mode)
  curve_fit(dat, use_VPD = F, Time_var = "Time", E_var = "E_corr", results = results_fit, use_acclim = F, col.per = ColorsPeriod, export_PPTX = T, remove_mode = "day_or_night")
  curve_fit(dat, use_VPD = F, Time_var = "Time", E_var = "E_corr", results = results_fit_acclim, use_acclim = T, col.per = ColorsPeriod_acclim, export_PPTX = T, add_points = F, acclim_only = T, remove_mode = "day_or_night")
  
  # Add legend
  legend(X_LEGEND, Y_LEGEND, bty = "n", cex = 0.7,
         #title = LEGEND_TITLE, title.adj = 0,
         legend = c("fit without acclimation", "fit with 24-h acclimation"),
         lty = c(1, 1), lwd = c(1, 4),
         col = rgb(t(col2rgb("gray30"))/255, alpha = c(1, 0.5)))
  text(4.6, 1.8, LEGEND_TITLE, adj = 0, cex = 0.9)
  }


# Plot in PPTX files

wd <- 6.5 # fig width
ht <- 2 # fig height

# Col-0
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_individual_fit("242", LEGEND_TITLE = "Col-0 - Pot #242", draw_x_axis = F, draw_rect_env = T)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "fit_Col-a.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_individual_fit("275", LEGEND_TITLE = "Col-0 - Pot #275", draw_x_axis = F, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "fit_Col-b.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_individual_fit("368", LEGEND_TITLE = "Col-0 - Pot #368", draw_x_axis = T, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "fit_Col-c.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()

# pgm
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_individual_fit("233", LEGEND_TITLE = expression(paste(italic("pgm"), " - Pot #233", sep = "")), draw_x_axis = F, draw_rect_env = T)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "fit_pgm-a.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_individual_fit("272", LEGEND_TITLE = expression(paste(italic("pgm"), " - Pot #272", sep = "")), draw_x_axis = F, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "fit_pgm-b.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_individual_fit("324", LEGEND_TITLE = expression(paste(italic("pgm"), " - Pot #324", sep = "")), draw_x_axis = T, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "fit_pgm-c.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()


# sex1
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_individual_fit("226", LEGEND_TITLE = expression(paste(italic("sex1"), " - Pot #226", sep = "")), draw_x_axis = F, draw_rect_env = T)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "fit_sex1-a.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_individual_fit("271", LEGEND_TITLE = expression(paste(italic("sex1"), " - Pot #271", sep = "")), draw_x_axis = F, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "fit_sex1-b.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_individual_fit("316", LEGEND_TITLE = expression(paste(italic("sex1"), " - Pot #316", sep = "")), draw_x_axis = T, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "fit_sex1-c.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()


