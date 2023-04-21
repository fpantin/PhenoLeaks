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


options(repos = "https://cran.rstudio.com/") # RStudio
for (pkg in c("here", "ggtext"))
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

XLIM <- c(-0.75, 4.25)
YLIM <- c(0.05, 2.5)
MAR <- c(2.5, 3.5, 0.5, 3.5)
X.LEGEND <- 2.75
Y.LEGEND <- 2.3


## 1 - Non-normalized data

# pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03c_plot_average_per_condition.pdf", sep = "_")), width = 8, height = 5)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "pgm - WW"), text.legend = c("Col-0", expression(italic("pgm"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, use_VPD = F, x.legend = 2.75, y.legend = 2.3, mar = MAR)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "sex1-3 - WW"), text.legend = c("Col-0", expression(italic("sex1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, use_VPD = F, x.legend = 2.75, y.legend = 2.3, mar = MAR)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "bam1 bam3 - WW"), text.legend = c("Col-0", expression(italic("bam1 bam3"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, use_VPD = F, x.legend = 2.75, y.legend = 2.3, mar = MAR)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "bam1 - WW"), text.legend = c("Col-0", expression(italic("bam1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, use_VPD = F, x.legend = 2.75, y.legend = 2.3, mar = MAR)
# plot_avg_kin(Tr.avg, c("Col-0 - WW", "bam3 - WW"), text.legend = c("Col-0", expression(italic("bam3"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, use_VPD = F, x.legend = 2.75, y.legend = 2.3, mar = MAR)
# dev.off()

common.args <- list(avg_dat = Tr.avg,
                    use_VPD = F,
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR,
                    x.legend = X.LEGEND, y.legend = Y.LEGEND)

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03c_plot_average_per_condition.pdf", sep = "_")), width = 8, height = 5)

do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "pgm - WW"),
                                               text.legend = c("Col-0", expression(italic("pgm"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "sex1-3 - WW"),
                                               text.legend = c("Col-0", expression(italic("sex1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam1 bam3 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam1 bam3"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam1 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam3 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam3"))))))

dev.off()


# Plot the kinetics of all genotypes into separate PPTX files (for publication)
wd <- 6.5 # fig width
ht <- 2 # fig height
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "pgm - WW", "sex1-3 - WW"),
             text.legend = c("Col-0 - moderate to fast incipient soil drying (n = 7)",
                             expression(paste(italic("pgm"), " - slow to moderate incipient soil drying (n = 7)")),
                             expression(paste(italic("sex1"), " - slow incipient soil drying (n = 4)"))),
             idExperiment = idExp, xlim = c(-0.7, 5.9), ylim = c(0, 2.4), export_PPTX = T, draw_x_axis = F,
             cex.legend = 0.6, x.legend = 2.75, y.legend = 2.3, irrig_show_mode = "none")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_Col_pgm_sex1.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - WW", "bam3 - WW", "bam1 bam3 - WW"),
             text.legend = c("Col-0 - moderate to fast incipient soil drying (n = 7)",
                             expression(paste(italic("bam3"), " - slow incipient soil drying (n = 4)")),
                             expression(paste(italic("bam1 bam3"), " - slow to moderate incipient soil drying (n = 5)"))),
             idExperiment = idExp, xlim = c(-0.7, 5.9), ylim = c(0, 2.4), export_PPTX = T, draw_x_axis = T, draw_rect_env = F,
             cex.legend = 0.6, x.legend = 2.75, y.legend = 2.3, irrig_show_mode = "none")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_Col_bam3_bam1bam3.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()



## 2 - Data normalized by the VPD

# pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03d_plot_average_per_condition_use_VPD.pdf", sep = "_")), width = 8, height = 5)
# plot_avg_kin(Tr.avg.kPa, c("Col-0 - WW", "pgm - WW"), text.legend = c("Col-0", expression(italic("pgm"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, use_VPD = T, x.legend = 2.75, y.legend = 2.3, mar = MAR)
# plot_avg_kin(Tr.avg.kPa, c("Col-0 - WW", "sex1-3 - WW"), text.legend = c("Col-0", expression(italic("sex1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, use_VPD = T, x.legend = 2.75, y.legend = 2.3, mar = MAR)
# plot_avg_kin(Tr.avg.kPa, c("Col-0 - WW", "bam1 bam3 - WW"), text.legend = c("Col-0", expression(italic("bam1 bam3"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, use_VPD = T, x.legend = 2.75, y.legend = 2.3, mar = MAR)
# plot_avg_kin(Tr.avg.kPa, c("Col-0 - WW", "bam1 - WW"), text.legend = c("Col-0", expression(italic("bam1"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, use_VPD = T, x.legend = 2.75, y.legend = 2.3, mar = MAR)
# plot_avg_kin(Tr.avg.kPa, c("Col-0 - WW", "bam3 - WW"), text.legend = c("Col-0", expression(italic("bam3"))), idExperiment = idExp, xlim = XLIM, ylim = YLIM, use_VPD = T, x.legend = 2.75, y.legend = 2.3, mar = MAR)
# dev.off()

common.args <- list(avg_dat = Tr.avg.kPa,
                    use_VPD = T,
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR,
                    x.legend = X.LEGEND, y.legend = Y.LEGEND)

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03d_plot_average_per_condition_use_VPD.pdf", sep = "_")), width = 8, height = 5)

do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "pgm - WW"),
                                               text.legend = c("Col-0", expression(italic("pgm"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "sex1-3 - WW"),
                                               text.legend = c("Col-0", expression(italic("sex1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam1 bam3 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam1 bam3"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam1 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam1"))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - WW", "bam3 - WW"),
                                               text.legend = c("Col-0", expression(italic("bam3"))))))

dev.off()


#------------------------------------------------------------------------------#
#            Same procedure but account for the severity of SWC drop           #
#------------------------------------------------------------------------------#

# Assign the severity of SWC drop
meanSWC <- aggregate(list(SWC = df$SWC), by = list(idPot = df$idPot), FUN = mean, na.rm = T)
meanSWC[order(meanSWC$SWC), ]
minSWC <- aggregate(list(SWC = df$SWC), by = list(idPot = df$idPot), FUN = min, na.rm = T)
minSWC[order(minSWC$SWC), ]

df$SWC_drop <- "moderate"
fast <- minSWC$idPot[minSWC$SWC <= 0.8]
df$SWC_drop[df$idPot %in% fast] <- "fast"
slow <- minSWC$idPot[minSWC$SWC > 1.2]
df$SWC_drop[df$idPot %in% slow] <- "slow"

unique(df$idGenotype[df$SWC_drop == "fast"]) # only Col-0 and bam1
unique(df$idGenotype[df$SWC_drop == "moderate"]) # only Col-0, pgm and bam1 bam3
unique(df$idGenotype[df$SWC_drop == "slow"]) # no Col-0

# Check the relationship with rosette surface
surf <- read.csv(file.path(dir_Exp, "Processed_data", "C2M43B_pot_transpiration.csv"))
minSWC_inisurf <- data.frame(idPot = NULL, idGenotype = NULL, SWC_drop = NULL, minSWC = NULL, inisurf = NULL)
for (pot in unique(minSWC$idPot))
  {
  minSWC_inisurf <- rbind(minSWC_inisurf,
                          data.frame(idPot = pot,
                                     idGenotype = df$idGenotype[df$idPot == pot][1],
                                     SWC_drop = df$SWC_drop[df$idPot == pot][1],
                                     minSWC = minSWC$SWC[minSWC$idPot == pot],
                                     inisurf = min(surf$surface[surf$idPot == pot], na.rm = T)/100))
  }
rm(surf)
minSWC_inisurf$idGenotype <- factor(minSWC_inisurf$idGenotype,
                                    levels = c("Col-0", "pgm", "sex1-3", "bam1", "bam3", "bam1 bam3"), ordered = T)
plot(minSWC ~ inisurf, data = minSWC_inisurf, type = "n", log = "y",
     xlab = expression(paste("Initial rosette area (cm"^2, ")", sep = "")),
     ylab = expression(paste("Lowest soil water content (g"["water"], " g"["dry soil"]^"-1", ")", sep = "")))
points(minSWC ~ inisurf, data = minSWC_inisurf)
abline(h = c(0.8, 1.2), lty = 2)

# For publication:
#require(ggplot2)
#require(ggtext)
#require(export)
wd <- 4 # fig width
ht <- 3 # fig height
if (!dir.exists(file.path(dir_Exp, "Figures", "PPTX"))) { dir.create(file.path(dir_Exp, "Figures", "PPTX")) }

dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
ggplot(minSWC_inisurf[minSWC_inisurf$idGenotype != "bam1", ],
       aes(y = minSWC , x = inisurf,
           shape = SWC_drop, colour = idGenotype)) + 
  geom_point(aes(shape = SWC_drop, colour = idGenotype), size = 2) +
  #scale_y_continuous(trans = "log10") +
  scale_shape_manual(values = c(16, 17, 18),
                     name = "Soil drying") +
  scale_color_manual(values = ColorsTrt$col[match(levels(minSWC_inisurf$idGenotype)[!levels(minSWC_inisurf$idGenotype) %in% "bam1"], ColorsTrt$idGenotype)],
                     name = "Genotype",
                     labels = c("Col-0", "*pgm*", "*sex1*", "*bam3*", "*bam1 bam3*")) +
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 3)),
         shape = guide_legend(override.aes = list(colour = "grey"))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_markdown(),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black")) +
  labs(x = expression(paste("Initial rosette area (cm"^2, ")", sep = "")),
       y = expression(paste("Lowest soil water content (g"["water"], " g"["dry soil"]^"-1", ")", sep = ""))) +
  geom_hline(yintercept = c(0.8, 1.2), linetype = "dashed", colour = "grey")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "Soil_drying_rosette_area.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()


## 1a - Non-normalized data

# Get average kinetics and error intervals
df$Trt <- paste(df$idGenotype, df$SWC_drop, sep = " - ")
Tr.avg <- compute_avg(df, E_var = "E_corr", Trt_var = "Trt") 

# Correct the average (not the error bars!) for artifactual "jumps" due to missing data
Tr.avg <- offset_avg(df, Tr.avg, interval_est = 0.15)

# Estimate the average SWC
Tr.avg <- append_SWC_avg(df, Tr.avg)

# How many replicates?
(Nmean <- aggregate(list(N = Tr.avg$N), by = list(Trt = Tr.avg$Trt), FUN = mean))
(Nmax <- aggregate(list(N = Tr.avg$N), by = list(Trt = Tr.avg$Trt), FUN = max))

# Plot the average kinetics and error interval
YLIM <- c(0.05, 2.8)
X.LEGEND <- 2.5
Y.LEGEND <- 2.6
common.args <- list(avg_dat = Tr.avg,
                    use_VPD = F,
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR,
                    x.legend = X.LEGEND, y.legend = Y.LEGEND)

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03e_plot_average_per_condition_use_SWC.pdf", sep = "_")), width = 8, height = 5)
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - moderate", "Col-0 - fast"),
                                               color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(2, 3),
                                               title_by_SWC = "Col-0",
                                               text.legend = c("moderate soil drying (n = 3)", "fast soil drying (n = 4)"))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("pgm - slow", "pgm - moderate"),
                                               color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(1, 2),
                                               title_by_SWC = expression(italic("pgm")),
                                               text.legend = c("slow soil drying (n = 4)", "moderate soil drying (n = 3)"))))
do.call(plot_avg_kin, append(common.args, list(treatments = "sex1-3 - slow",
                                               color_E_by_SWC = T, col_by_SWC = "black", lty_by_SWC = 1,
                                               title_by_SWC = expression(italic("sex1")),
                                               text.legend = "slow soil drying (n = 4)")))
do.call(plot_avg_kin, append(common.args, list(treatments = c("bam1 bam3 - slow", "bam1 bam3 - moderate"),
                                               color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(1, 2),
                                               title_by_SWC = expression(italic("bam1 bam3")),
                                               text.legend = c("slow soil drying (n = 1)", "moderate soil drying (n = 4)"))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("bam1 - slow", "bam1 - fast"),
                                               color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(1, 3),
                                               title_by_SWC = expression(italic("bam1")),
                                               text.legend = c("slow soil drying (n = 1)", "fast soil drying (n = 2)"))))
do.call(plot_avg_kin, append(common.args, list(treatments = "bam3 - slow",
                                               color_E_by_SWC = T, col_by_SWC = "black", lty_by_SWC = 1,
                                               title_by_SWC = expression(italic("bam3")),
                                               text.legend = "slow soil drying (n = 4)")))
dev.off()


# Plot Col-0 and pgm into separate PPTX files (for publication)
wd <- 6.5 # fig width
ht <- 2 # fig height
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("Col-0 - moderate", "Col-0 - fast"), text.legend = c("Col-0 - moderate incipient soil drying (n = 3)",
                                                                            "Col-0 - fast incipient soil drying (n = 4)"),
             idExperiment = idExp, xlim = c(-0.7, 5.9), ylim = c(0, 2.4), export_PPTX = T, draw_x_axis = T,
             color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(2, 3), title_by_SWC = "Experiment #3",#"Col-0",
             scale_bar_pos_by_SWC = c(-0.25, 0.3), cex.legend = 0.6, x.legend = 2.75, y.legend = 2.2)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_WS_Col.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()
dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
plot_avg_kin(Tr.avg, c("pgm - slow", "pgm - moderate"), text.legend =  c(expression(paste(italic("pgm"), " - slow incipient soil drying (n = 4)")),
                                                                         expression(paste(italic("pgm"), " - moderate incipient soil drying (n = 3)"))),
             idExperiment = idExp, xlim = c(-0.7, 5.9), ylim = c(0, 2.4), export_PPTX = T, draw_x_axis = T,
             color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(1, 2), title_by_SWC = "Experiment #3",#expression(italic("pgm")),
             scale_bar_pos_by_SWC = c(-0.25, 0.3), cex.legend = 0.6, x.legend = 2.75, y.legend = 2.2)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_WS_pgm.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
dev.off()



## 1b - Non-normalized data - Pool slow and moderate drying

# Get average kinetics and error intervals
df$SWC_drop_pool <- df$SWC_drop
df$SWC_drop_pool[df$SWC_drop %in% c("slow", "moderate")] <- "slow_moderate"
df$Trt <- paste(df$idGenotype, df$SWC_drop_pool, sep = " - ")
Tr.avg <- compute_avg(df, E_var = "E_corr", Trt_var = "Trt") 

# Correct the average (not the error bars!) for artifactual "jumps" due to missing data
Tr.avg <- offset_avg(df, Tr.avg, interval_est = 0.15)

# How many replicates?
(Nmean <- aggregate(list(N = Tr.avg$N), by = list(Trt = Tr.avg$Trt), FUN = mean))
(Nmax <- aggregate(list(N = Tr.avg$N), by = list(Trt = Tr.avg$Trt), FUN = max))

# Update the color dataframe (note that the "fast" soil drying will not be plotted)
ColorsTrt$Trt <- paste(ColorsTrt$idGenotype, "slow_moderate", sep = " - ")

# Plot the average kinetics and error interval
X.LEGEND <- 2.75
common.args <- list(avg_dat = Tr.avg,
                    use_VPD = F,
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR,
                    x.legend = X.LEGEND, y.legend = Y.LEGEND)

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03f_plot_average_per_condition_use_SWC_poolLowModerate.pdf", sep = "_")), width = 8, height = 5)
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - slow_moderate", "pgm - slow_moderate"),
                                               text.legend = c("Col-0 (n = 3)",
                                                               expression(paste(italic("pgm"), " (n = 7)", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - slow_moderate", "sex1-3 - slow_moderate"),
                                               text.legend = c("Col-0 (n = 3)",
                                                               expression(paste(italic("sex1"), " (n = 4)", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - slow_moderate", "bam1 bam3 - slow_moderate"),
                                               text.legend = c("Col-0 (n = 3)",
                                                               expression(paste(italic("bam1 bam3"), " (n = 5)", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - slow_moderate", "bam1 - slow_moderate"),
                                               text.legend = c("Col-0 (n = 3)",
                                                               expression(paste(italic("bam1"), " (n = 1)", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - slow_moderate", "bam3 - slow_moderate"),
                                               text.legend = c("Col-0 (n = 3)",
                                                               expression(paste(italic("pgm"), " (n = 4)", sep = ""))))))
dev.off()


# # Plot the kinetics into separate PPTX files (for publication)
# wd <- 6.5 # fig width
# ht <- 2 # fig height
# dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
# plot_avg_kin(Tr.avg, c("Col-0 - slow_moderate", "pgm - slow_moderate", "sex1-3 - slow_moderate"),
#              text.legend = c("Col-0 - moderate incipient soil drying (n = 3)",
#                              expression(paste(italic("pgm"), " - slow to moderate incipient soil drying (n = 7)")),
#                              expression(paste(italic("sex1"), " - slow incipient soil drying (n = 4)"))),
#              idExperiment = idExp, xlim = c(-0.7, 5.9), ylim = c(0, 2.4), export_PPTX = T, draw_x_axis = F,
#              cex.legend = 0.6, x.legend = 2.75, y.legend = 2.3, irrig_show_mode = "none")
# graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_Col_pgm_sex1.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
# dev.off()
# dev.new(width = wd, height = ht, unit = "in", noRStudioGD = T)
# plot_avg_kin(Tr.avg, c("Col-0 - slow_moderate", "bam3 - slow_moderate", "bam1 bam3 - slow_moderate"),
#              text.legend = c("Col-0 - moderate incipient soil drying (n = 3)",
#                              expression(paste(italic("bam3"), " - slow incipient soil drying (n = 4)")),
#                              expression(paste(italic("bam1 bam3"), " - slow to moderate incipient soil drying (n = 5)"))),
#              idExperiment = idExp, xlim = c(-0.7, 5.9), ylim = c(0, 2.4), export_PPTX = T, draw_x_axis = T, draw_rect_env = F,
#              cex.legend = 0.6, x.legend = 2.75, y.legend = 2.3, irrig_show_mode = "none")
# graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin_Col_bam3_bam1bam3.pptx", sep = "_")), paper = "A4", orient = "portrait", width = wd, height = ht)
# dev.off()





## 2a - Data normalized by the VPD

# Get average kinetics and error intervals
df$Trt <- paste(df$idGenotype, df$SWC_drop, sep = " - ")
Tr.avg.kPa <- compute_avg(df, E_var = "E_corr_per_kPa", Trt_var = "Trt")

# Correct the average (not the error bars!) for artifactual "jumps" due to missing data
Tr.avg.kPa <- offset_avg(df, Tr.avg.kPa, interval_est = 0.15)

# Estimate the average VPD
Tr.avg.kPa <- append_VPD_avg(df, Tr.avg.kPa)

# Estimate the average SWC
Tr.avg.kPa <- append_SWC_avg(df, Tr.avg.kPa)

# How many replicates?
(Nmean <- aggregate(list(N = Tr.avg.kPa$N), by = list(Trt = Tr.avg.kPa$Trt), FUN = mean))
(Nmax <- aggregate(list(N = Tr.avg.kPa$N), by = list(Trt = Tr.avg.kPa$Trt), FUN = max))

# Plot the average kinetics and error interval
YLIM <- c(0.05, 2.8)
X.LEGEND <- 2.5
Y.LEGEND <- 2.6
common.args <- list(avg_dat = Tr.avg.kPa,
                    use_VPD = T, add_VPD_mode = "single",
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR,
                    x.legend = X.LEGEND, y.legend = Y.LEGEND)

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03g_plot_average_per_condition_use_SWC_use_VPD.pdf", sep = "_")), width = 8, height = 5)
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - moderate", "Col-0 - fast"),
                                               color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(2, 3),
                                               title_by_SWC = "Col-0",
                                               text.legend = c("moderate soil drying (n = 3)", "fast soil drying (n = 4)"))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("pgm - slow", "pgm - moderate"),
                                               color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(1, 2),
                                               title_by_SWC = expression(italic("pgm")),
                                               text.legend = c("slow soil drying (n = 4)", "moderate soil drying (n = 3)"))))
do.call(plot_avg_kin, append(common.args, list(treatments = "sex1-3 - slow",
                                               color_E_by_SWC = T, col_by_SWC = "black", lty_by_SWC = 1,
                                               title_by_SWC = expression(italic("sex1")),
                                               text.legend = "slow soil drying (n = 4)")))
do.call(plot_avg_kin, append(common.args, list(treatments = c("bam1 bam3 - slow", "bam1 bam3 - moderate"),
                                               color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(1, 2),
                                               title_by_SWC = expression(italic("bam1 bam3")),
                                               text.legend = c("slow soil drying (n = 1)", "moderate soil drying (n = 4)"))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("bam1 - slow", "bam1 - fast"),
                                               color_E_by_SWC = T, col_by_SWC = c("black", "black"), lty_by_SWC = c(1, 3),
                                               title_by_SWC = expression(italic("bam1")),
                                               text.legend = c("slow soil drying (n = 1)", "fast soil drying (n = 2)"))))
do.call(plot_avg_kin, append(common.args, list(treatments = "bam3 - slow",
                                               color_E_by_SWC = T, col_by_SWC = "black", lty_by_SWC = 1,
                                               title_by_SWC = expression(italic("bam3")),
                                               text.legend = "slow soil drying (n = 4)")))
dev.off()


## 2b - Data normalized by the VPD - Pool slow and moderate drying

# Get average kinetics and error intervals
df$SWC_drop_pool <- df$SWC_drop
df$SWC_drop_pool[df$SWC_drop %in% c("slow", "moderate")] <- "slow_moderate"
df$Trt <- paste(df$idGenotype, df$SWC_drop_pool, sep = " - ")
Tr.avg.kPa <- compute_avg(df, E_var = "E_corr_per_kPa", Trt_var = "Trt") 

# Correct the average (not the error bars!) for artifactual "jumps" due to missing data
Tr.avg.kPa <- offset_avg(df, Tr.avg.kPa, interval_est = 0.15)

# Estimate the average VPD
Tr.avg.kPa <- append_VPD_avg(df, Tr.avg.kPa)

# How many replicates?
(Nmean <- aggregate(list(N = Tr.avg.kPa$N), by = list(Trt = Tr.avg.kPa$Trt), FUN = mean))
(Nmax <- aggregate(list(N = Tr.avg.kPa$N), by = list(Trt = Tr.avg.kPa$Trt), FUN = max))

# Update the color dataframe (note that the "fast" soil drying will not be plotted)
ColorsTrt$Trt <- paste(ColorsTrt$idGenotype, "slow_moderate", sep = " - ")

# Plot the average kinetics and error interval
X.LEGEND <- 2.75
common.args <- list(avg_dat = Tr.avg.kPa,
                    use_VPD = T, add_VPD_mode = "single",
                    idExperiment = idExp,
                    xlim = XLIM, ylim = YLIM, mar = MAR,
                    x.legend = X.LEGEND, y.legend = Y.LEGEND)

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03h_plot_average_per_condition_use_SWC_poolLowModerate_use_VPD.pdf", sep = "_")), width = 8, height = 5)
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - slow_moderate", "pgm - slow_moderate"),
                                               text.legend = c("Col-0 (n = 3)",
                                                               expression(paste(italic("pgm"), " (n = 7)", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - slow_moderate", "sex1-3 - slow_moderate"),
                                               text.legend = c("Col-0 (n = 3)",
                                                               expression(paste(italic("sex1"), " (n = 4)", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - slow_moderate", "bam1 bam3 - slow_moderate"),
                                               text.legend = c("Col-0 (n = 3)",
                                                               expression(paste(italic("bam1 bam3"), " (n = 5)", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - slow_moderate", "bam1 - slow_moderate"),
                                               text.legend = c("Col-0 (n = 3)",
                                                               expression(paste(italic("bam1"), " (n = 1)", sep = ""))))))
do.call(plot_avg_kin, append(common.args, list(treatments = c("Col-0 - slow_moderate", "bam3 - slow_moderate"),
                                               text.legend = c("Col-0 (n = 3)",
                                                               expression(paste(italic("pgm"), " (n = 4)", sep = ""))))))
dev.off()
