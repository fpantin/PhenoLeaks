################################################################################
#                                                                              #
#               PhenoLeaks - Step #03 - Plotting average kinetics              #
#                                                                              #
#                   Script to plot the transpiration kinetics                  #
#                averaged per treatment (genotype x environment)               #
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
#                       Import and manage corrected data                       #
#------------------------------------------------------------------------------#


# Import corrected transpiration data
df <- read.csv(file.path(dir_Exp, "Corrected_data", "E_corr.csv"))


#------------------------------------------------------------------------------#
#                     Compute average kinetics by treatment                    #
#------------------------------------------------------------------------------#

## 1 - Non-normalized data

# Get average kinetics and error intervals
Tr.avg <- compute_avg(df, E_var = "E_corr", Trt_var = "Trt") 

# Correct the average (not the error bars!) for artifactual "jumps" due to missing data
Tr.avg <- offset_avg(df, Tr.avg, interval_est = 0.15)

# Show the correction in a PDF file
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03a_plot_missing_data_correction.pdf", sep = "_")), width = 8, height = 5)
par(mfrow = c(2, 1), oma = c(0, 0, 2, 10))
for (trt in sort(unique(df$Trt)))
  {
  Delta <- Tr.avg[Tr.avg$Trt == trt, ]
  Delta$Delta <- Tr.avg$E_corr_avg[Tr.avg$Trt == trt] - Tr.avg$E_corr_offset_avg[Tr.avg$Trt == trt]
  
  prepare_kin(Delta, Time_var = "Time", E_var = "E_corr_avg", mar = c(3, 4, 0, 0), main = "", draw_rect_env = F)
  points(E_corr_avg ~ Time, data = Delta, cex = 0.5)
  points(E_corr_offset_avg ~ Time, data = Delta, col = 2, pch = 16, cex = 0.5)
  
  mtext(trt, side = 3, outer = T, cex = 1.5)
  
  prepare_kin(Delta, Time_var = "Time", E_var = "Delta", mar = c(3, 4, 0, 0), main = "", draw_rect_env = F)
  points(Delta ~ Time, data = Delta, col = 3, cex = 0.5)
  
  legend("topright", inset = c(-0.35, -0.25), legend = c("Mean of available data", "Offset correction for missing data", "Difference"),
         col = 1:3, pch = c(1, 16, 1), pt.cex = 0.5, bty = "n", xpd = NA, cex = 0.5)
  }
dev.off()

# Save the results
write.csv(Tr.avg, file.path(dir_Exp, "Corrected_data", "E_corr_offset_avg.csv"), row.names = F)


## 2 - Data normalized by the VPD

# Get average kinetics and error intervals
Tr.avg.kPa <- compute_avg(df, E_var = "E_corr_per_kPa", Trt_var = "Trt")

# Correct the average (not the error bars!) for artifactual "jumps" due to missing data
Tr.avg.kPa <- offset_avg(df, Tr.avg.kPa, interval_est = 0.15)

# Estimate the average VPD
Tr.avg.kPa <- append_VPD_avg(df, Tr.avg.kPa)

# Show the correction in a PDF file
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03b_plot_missing_data_correction_use_VPD.pdf", sep = "_")), width = 8, height = 5)
par(mfrow = c(2, 1), oma = c(0, 0, 2, 10))
for (trt in sort(unique(df$Trt)))
  {
  Delta <- Tr.avg.kPa[Tr.avg.kPa$Trt == trt, ]
  Delta$Delta <- Tr.avg.kPa$E_corr_per_kPa_avg[Tr.avg.kPa$Trt == trt] - Tr.avg.kPa$E_corr_per_kPa_offset_avg[Tr.avg.kPa$Trt == trt]
  
  prepare_kin(Delta, use_VPD = T, Time_var = "Time", E_var = "E_corr_per_kPa_avg", mar = c(3, 4, 0, 0), main = "", draw_rect_env = F)
  points(E_corr_per_kPa_avg ~ Time, data = Delta, cex = 0.5)
  points(E_corr_per_kPa_offset_avg ~ Time, data = Delta, col = 2, pch = 16, cex = 0.5)
  
  mtext(trt, side = 3, outer = T, cex = 1.5)
  
  prepare_kin(Delta, use_VPD = T, add_VPD = F, Time_var = "Time", E_var = "Delta", mar = c(3, 4, 0, 0), main = "", draw_rect_env = F)
  points(Delta ~ Time, data = Delta, col = 3, cex = 0.5)
  
  legend("topright", inset = c(-0.35, -0.25), legend = c("Mean of available data", "Offset correction for missing data", "Difference"),
         col = 1:3, pch = c(1, 16, 1), pt.cex = 0.5, bty = "n", xpd = NA, cex = 0.5)
  }
dev.off()


# Save the results
write.csv(Tr.avg.kPa, file.path(dir_Exp, "Corrected_data", "E_corr_per_kPa_offset_avg.csv"), row.names = F)


#------------------------------------------------------------------------------#
#                 Plot the average kinetics and error interval                 #
#------------------------------------------------------------------------------#

XLIM <- c(-0.7, 4)
YLIM <- c(0.8, 3.5)
MAR <- c(2.5, 3.5, 0.5, 3.5)
X.LEGEND <- 1
Y.LEGEND <- 3.4
Y.AXIS.BREAK <- 0.8

## 1 - Non-normalized data

common.args <- list(avg_dat = Tr.avg,
                    use_VPD = F,
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR,
                    x.legend = X.LEGEND, y.legend = Y.LEGEND,
                    y.axis.break = Y.AXIS.BREAK)

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03c_plot_average_per_condition.pdf", sep = "_")), width = 8, height = 5)

do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW"),
                                               text.legend = c("Col-0"))))

dev.off()


## 2 - Data normalized by the VPD

# pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03d_plot_average_per_condition_use_VPD.pdf", sep = "_")), width = 8, height = 5)
# plot_avg_kin(Tr.avg.kPa, c("Col-0 - WW"), text.legend = c("Col-0"), idExperiment = idExp, xlim = XLIM, ylim = YLIM, use_VPD = T, x.legend = 1, y.legend = 3.4, y.axis.break = 0.8)
# dev.off()

common.args <- list(avg_dat = Tr.avg.kPa,
                    use_VPD = T,
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR,
                    x.legend = X.LEGEND, y.legend = Y.LEGEND,
                    y.axis.break = Y.AXIS.BREAK)

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03d_plot_average_per_condition_use_VPD.pdf", sep = "_")), width = 8, height = 5)

do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW"),
                                               text.legend = c("Col-0"))))
dev.off()

