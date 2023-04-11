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
#                                     C2M47                                    #
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
#                       Import and manage corrected data                       #
#------------------------------------------------------------------------------#


# Import corrected transpiration data
df <- read.csv(file.path(dir_Exp, "Corrected_data", "E_corr.csv"))

# Keep only the time interval that will be used for averaging
df <- df[df$Time >= Time_start_exp, ]


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
    
  prepare_kin(Delta, Time_var = "Time", E_var = "E_corr_avg", xlim = c(-0.7, 5.9), mar = c(3, 4, 0, 0), main = "", draw_rect_env = F)
  points(E_corr_avg ~ Time, data = Delta, cex = 0.5)
  points(E_corr_offset_avg ~ Time, data = Delta, col = 2, pch = 16, cex = 0.5)
  
  mtext(trt, side = 3, outer = T, cex = 1.5)
  
  prepare_kin(Delta, Time_var = "Time", E_var = "Delta", xlim = c(-0.7, 5.9), mar = c(3, 4, 0, 0), main = "", draw_rect_env = F)
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
  
  prepare_kin(Delta, use_VPD = T, Time_var = "Time", E_var = "E_corr_per_kPa_avg", xlim = c(-0.7, 5.9), mar = c(3, 4, 0, 0), main = "", draw_rect_env = F)
  points(E_corr_per_kPa_avg ~ Time, data = Delta, cex = 0.5)
  points(E_corr_per_kPa_offset_avg ~ Time, data = Delta, col = 2, pch = 16, cex = 0.5)
  
  mtext(trt, side = 3, outer = T, cex = 1.5)
  
  prepare_kin(Delta, use_VPD = T, add_VPD = F, Time_var = "Time", E_var = "Delta", xlim = c(-0.7, 5.9), mar = c(3, 4, 0, 0), main = "", draw_rect_env = F)
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

XLIM <- c(-0.7, 5.9)
YLIM <- c(0.35, 1.85)
MAR <- c(2.5, 3.5, 0.5, 3.5)


## 1 - Non-normalized data

# pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03c_plot_average_per_condition.pdf", sep = "_")), width = 8, height = 5)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "Col-0 ProStarv::LUC - WW"), text.legend = c("Col-0", "Col-0*"), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "pgm-1 - WW"), text.legend = c("Col-0", expression(italic("pgm"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "sex1-3 - WW"), text.legend = c("Col-0", expression(italic("sex1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 ProStarv::LUC - WW", "isa1-1 ProStarv::LUC - WW"), text.legend = c("Col-0*", expression(paste(italic("isa1"), "*", sep = ""))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 ProStarv::LUC - WW", "ss4-3 ProStarv::LUC - WW"), text.legend = c("Col-0*", expression(paste(italic("ss4"), "*", sep = ""))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "mex1-1 - WW"), text.legend = c("Col-0", expression(italic("mex1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "dpe1-2 - WW"), text.legend = c("Col-0", expression(italic("dpe1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "dpe2-5 - WW"), text.legend = c("Col-0", expression(italic("dpe2"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "pgi1-1 - WW"), text.legend = c("Col-0", expression(italic("pgi"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "amy3-2 bam1 - WW"), text.legend = c("Col-0", expression(italic("amy3 bam1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "bam1bam3 - WW"), text.legend = c("Col-0", expression(italic("bam1 bam3"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "amy3-2 - WW"), text.legend = c("Col-0", expression(italic("amy3"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "bam1 - WW"), text.legend = c("Col-0", expression(italic("bam1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "bam3 - WW"), text.legend = c("Col-0", expression(italic("bam3"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "abcb14-1 - WW"), text.legend = c("Col-0", expression(paste(italic("abcb14"), "-1", sep = ""))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "abcb14-2 - WW"), text.legend = c("Col-0", expression(paste(italic("abcb14"), "-2", sep = ""))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("abcb14-1 - WW", "abcb14-2 - WW"), text.legend = c(expression(paste(italic("abcb14"), "-1", sep = "")), expression(paste(italic("abcb14"), "-2", sep = ""))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "Col-0 - WS"), text.legend = c("Col-0 well-watered", "Col-0 water stress"), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("pgm-1 - WW", "pgm-1 - WS"), text.legend = c(expression(paste(italic("pgm"), " well-watered")), expression(paste(italic("pgm"), " water stress"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# plot_avg_kin(Tr.avg, c("Col-0 - WS", "pgm-1 - WS"), text.legend = c("Col-0 water stress", expression(paste(italic("pgm"), " water stress"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM)
# dev.off()

common.args <- list(avg_dat = Tr.avg,
                    use_VPD = F,
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR)

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03c_plot_average_per_condition.pdf", sep = "_")), width = 8, height = 5)

do.call(plot_avg_kin, append(common.args, list(treatments =  c("Col-0 - WW", "Col-0 ProStarv::LUC - WW"),
                                               text.legend = c("Col-0", expression("Col-0"^"†")))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "pgm-1 - WW"),
                                               text.legend = c("Col-0", expression(italic("pgm"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "sex1-3 - WW"),
                                               text.legend = c("Col-0", expression(italic("sex1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 ProStarv::LUC - WW", "isa1-1 ProStarv::LUC - WW"),
                                               text.legend = c(expression("Col-0"^"†"),  expression(italic("isa1")^"†")))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 ProStarv::LUC - WW", "ss4-3 ProStarv::LUC - WW"),
                                               text.legend = c(expression("Col-0"^"†"), expression(italic("ss4")^"†")))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "mex1-1 - WW"),
                                               text.legend = c("Col-0", expression(italic("mex1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "dpe1-2 - WW"),
                                               text.legend = c("Col-0", expression(italic("dpe1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "dpe2-5 - WW"),
                                               text.legend = c("Col-0", expression(italic("dpe2"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "pgi1-1 - WW"),
                                               text.legend = c("Col-0", expression(italic("pgi"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "amy3-2 bam1 - WW"),
                                               text.legend = c("Col-0", expression(italic("amy3 bam1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam1bam3 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam1 bam3"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "amy3-2 - WW"),
                                               text.legend = c("Col-0", expression(italic("amy3"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam1 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam3 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam3"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "abcb14-1 - WW"),
                                               text.legend = c("Col-0", expression(paste(italic("abcb14"), "-1", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "abcb14-2 - WW"),
                                               text.legend = c("Col-0", expression(paste(italic("abcb14"), "-2", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("abcb14-1 - WW", "abcb14-2 - WW"),
                                               text.legend = c(expression(paste(italic("abcb14"), "-1", sep = "")), expression(paste(italic("abcb14"), "-2", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "Col-0 - WS"),
                                               text.legend = c("Col-0 well-watered", "Col-0 water stress"))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("pgm-1 - WW", "pgm-1 - WS"),
                                               text.legend = c(expression(paste(italic("pgm"), " well-watered")), expression(paste(italic("pgm"), " water stress"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WS", "pgm-1 - WS"),
                                               text.legend = c("Col-0 water stress", expression(paste(italic("pgm"), " water stress"))))))

dev.off()


# Plot the averages into separate PPTX files (for publication)
# Here we have set 'write_y_axis_label' to FALSE so that the name of the y-axis should be added manually in PowerPoint.
# This is because as some symbols are not currently supported by 'graph2ppt()' (see https://github.com/tomwenseleers/export/issues/49 for details).

wd <- 6.5 # fig width
ht <- 2 # fig height
if (!dir.exists(file.path(dir_Exp, "Figures", "PPTX"))) { dir.create(file.path(dir_Exp, "Figures", "PPTX")) }

dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "pgm-1 - WW"), text.legend = c("Col-0", expression(italic("pgm"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_x_axis = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_pgm.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "sex1-3 - WW"), text.legend = c("Col-0", expression(italic("sex1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_sex1.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()

dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 ProStarv::LUC - WW", "isa1-1 ProStarv::LUC - WW"), text.legend = c(expression("Col-0"^"†"), expression(italic("isa1")^"†")), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_x_axis = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_isa1.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 ProStarv::LUC - WW", "ss4-3 ProStarv::LUC - WW"), text.legend = c(expression("Col-0"^"†"), expression(italic("ss4")^"†")), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_ss4.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()

dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "mex1-1 - WW"), text.legend = c("Col-0", expression(italic("mex1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_x_axis = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_mex1.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "dpe1-2 - WW"), text.legend = c("Col-0", expression(italic("dpe1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_x_axis = F, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_dpe1.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "dpe2-5 - WW"), text.legend = c("Col-0", expression(italic("dpe2"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_dpe2.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()

dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "pgi1-1 - WW"), text.legend = c("Col-0", expression(italic("pgi"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_x_axis = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_pgi.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "amy3-2 bam1 - WW"), text.legend = c("Col-0", expression(italic("amy3 bam1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_x_axis = F, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_amy3 bam1.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "bam1bam3 - WW"), text.legend = c("Col-0", expression(italic("bam1 bam3"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_bam1 bam3.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()

dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "amy3-2 - WW"), text.legend = c("Col-0", expression(italic("amy3"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_x_axis = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_amy3.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "bam1 - WW"), text.legend = c("Col-0", expression(italic("bam1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_x_axis = F, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_bam1.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "bam3 - WW"), text.legend = c("Col-0", expression(italic("bam3"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_bam3.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()

dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "abcb14-1 - WW"), text.legend = c("Col-0", expression(paste(italic("abcb14"), "-1", sep = ""))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_x_axis = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_abcb14-1.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "abcb14-2 - WW"), text.legend = c("Col-0", expression(paste(italic("abcb14"), "-2", sep = ""))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, export_PPTX = T, draw_rect_env = F)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_abcb14-2.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()


# Zoom into the first complete diel cycle for Col-0, pgm and sex1 (for publication)
dev.new(width = 3, height = 2.7, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "pgm-1 - WW", "sex1-3 - WW"), text.legend = c("Col-0", expression(italic("pgm")), expression(italic("sex1"))), idExperiment = idExp, xlim = c(0, 1), ylim = c(0.65, 1.87), export_PPTX = T, draw_x_axis = T, x_axis_unit = "hours", xlab = "Time (h)", draw_rect_env = T, x.legend = 0.6)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin zoom_pgm_sex1.pptx", sep = "_")), paper = "A4", orient = "portrait", width = 3, height = 2.7)
dev.off()



## 2 - Data normalized by the VPD

YLIM <- c(0.45, 2.05)

common.args <- list(avg_dat = Tr.avg.kPa,
                    use_VPD = T,
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR)

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03d_plot_average_per_condition_use_VPD.pdf", sep = "_")), width = 8, height = 5)

do.call(plot_avg_kin, append(common.args, list(treatments =  c("Col-0 - WW", "Col-0 ProStarv::LUC - WW"),
                                               text.legend = c("Col-0", expression("Col-0"^"†")))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "pgm-1 - WW"),
                                               text.legend = c("Col-0", expression(italic("pgm"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "sex1-3 - WW"),
                                               text.legend = c("Col-0", expression(italic("sex1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 ProStarv::LUC - WW", "isa1-1 ProStarv::LUC - WW"),
                                               text.legend = c(expression("Col-0"^"†"), expression(italic("isa1")^"†")))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 ProStarv::LUC - WW", "ss4-3 ProStarv::LUC - WW"),
                                               text.legend = c(expression("Col-0"^"†"), expression(italic("ss4")^"†")))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "mex1-1 - WW"),
                                               text.legend = c("Col-0", expression(italic("mex1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "dpe1-2 - WW"),
                                               text.legend = c("Col-0", expression(italic("dpe1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "dpe2-5 - WW"),
                                               text.legend = c("Col-0", expression(italic("dpe2"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "pgi1-1 - WW"),
                                               text.legend = c("Col-0", expression(italic("pgi"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "amy3-2 bam1 - WW"),
                                               text.legend = c("Col-0", expression(italic("amy3 bam1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam1bam3 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam1 bam3"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "amy3-2 - WW"),
                                               text.legend = c("Col-0", expression(italic("amy3"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam1 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam3 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam3"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "abcb14-1 - WW"),
                                               text.legend = c("Col-0", expression(paste(italic("abcb14"), "-1", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "abcb14-2 - WW"),
                                               text.legend = c("Col-0", expression(paste(italic("abcb14"), "-2", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("abcb14-1 - WW", "abcb14-2 - WW"),
                                               text.legend = c(expression(paste(italic("abcb14"), "-1", sep = "")), expression(paste(italic("abcb14"), "-2", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "Col-0 - WS"),
                                               text.legend = c("Col-0 well-watered", "Col-0 water stress"))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("pgm-1 - WW", "pgm-1 - WS"),
                                               text.legend = c(expression(paste(italic("pgm"), " well-watered")), expression(paste(italic("pgm"), " water stress"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WS", "pgm-1 - WS"),
                                               text.legend = c("Col-0 water stress", expression(paste(italic("pgm"), " water stress"))))))

dev.off()



#------------------------------------------------------------------------------#
#             Plot E color-coded for SWC (Col-0 and pgm, WW vs.WS)             #
#------------------------------------------------------------------------------#

## 1 - Non-normalized data

# Estimate the average SWC
Tr.avg <- append_SWC_avg(df, Tr.avg)

# How many replicates?
(Nmean <- aggregate(list(N = Tr.avg$N), by = list(Trt = Tr.avg$Trt), FUN = mean))
(Nmax <- aggregate(list(N = Tr.avg$N), by = list(Trt = Tr.avg$Trt), FUN = max))

# Plot the average kinetics and error interval
X.LEGEND <- 2.5
Y.LEGEND <- 2.6
YLIM <- c(0.05, 2.8)
common.args <- list(avg_dat = Tr.avg,
                    use_VPD = F,
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR,
                    x.legend = X.LEGEND, y.legend = Y.LEGEND,
                    color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(1, 3))

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03e_plot_average_per_condition_use_SWC.pdf", sep = "_")), width = 8, height = 5)
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "Col-0 - WS"),
                                               title_by_SWC = "Col-0",
                                               text.legend = c("well-watered (n = 8)", "stable water stress (n = 8)"))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("pgm-1 - WW", "pgm-1 - WS"),
                                               title_by_SWC = expression(italic("pgm")),
                                               text.legend = c("well-watered (n = 8)", "stable water stress (n = 8)"))))
dev.off()


# Plot into separate PPTX files (for publication)
wd <- 6.5 # fig width
ht <- 2 # fig height
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "Col-0 - WS"), text.legend = c("Col-0 - stable well-watered regime (n = 8)",
                                                                    "Col-0 - stable soil water deficit (n = 8)"),
             idExperiment = idExp, xlim = c(-0.7, 5.9), ylim = c(0, 2.4), export_PPTX = T, draw_x_axis = F,
             color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(1, 4), title_by_SWC = "Experiment #4",#"Col-0",
             scale_bar_pos_by_SWC = c(-0.6, 0.3), cex.legend = 0.6, x.legend = 2.75, y.legend = 2.2, irrig_show_mode = "mean")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_WS_Col.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("pgm-1 - WW", "pgm-1 - WS"), text.legend = c(expression(paste(italic("pgm"), " - stable well-watered regime (n = 8)", sep = "")),
                                                                    expression(paste(italic("pgm"), " - stable soil water deficit (n = 8)", sep = ""))),
             idExperiment = idExp, xlim = c(-0.7, 5.9), ylim = c(0, 2.4), export_PPTX = T, draw_x_axis = F,
             color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(1, 4), title_by_SWC = "Experiment #4",#expression(italic("pgm")),
             scale_bar_pos_by_SWC = c(-0.6, 0.3), cex.legend = 0.6, x.legend = 2.75, y.legend = 2.2, irrig_show_mode = "mean")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_WS_pgm.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()


## 2 - Data normalized by the VPD

# Estimate the average SWC
Tr.avg.kPa <- append_SWC_avg(df, Tr.avg.kPa)

# How many replicates?
(Nmean <- aggregate(list(N = Tr.avg.kPa$N), by = list(Trt = Tr.avg.kPa$Trt), FUN = mean))
(Nmax <- aggregate(list(N = Tr.avg.kPa$N), by = list(Trt = Tr.avg.kPa$Trt), FUN = max))

# Plot the average kinetics and error interval
common.args <- list(avg_dat = Tr.avg.kPa,
                    use_VPD = T, add_VPD_mode = "single",
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR,
                    x.legend = X.LEGEND, y.legend = Y.LEGEND,
                    color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(1, 3))

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03f_plot_average_per_condition_use_SWC_use_VPD.pdf", sep = "_")), width = 8, height = 5)
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "Col-0 - WS"),
                                               title_by_SWC = "Col-0",
                                               text.legend = c("well-watered (n = 8)", "stable water stress (n = 8)"))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("pgm-1 - WW", "pgm-1 - WS"),
                                               title_by_SWC = expression(italic("pgm")),
                                               text.legend = c("well-watered (n = 8)", "stable water stress (n = 8)"))))
dev.off()
