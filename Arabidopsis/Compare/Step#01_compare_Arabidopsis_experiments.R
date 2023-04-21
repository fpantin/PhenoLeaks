################################################################################
#                                                                              #
#                PhenoLeaks - Comparing Arabidopsis experiments                #
#                                                                              #
#               Script to compare several Arabidopsis experiments              #
#                         (C3M31, C2M43A and B, C2M47)                         #
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
#                                    Compare                                   #
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
  source(file.path(here::here(), "Arabidopsis", "Compare", "Step#00_define_experiment.R"))

# Now the species and ID of the experiment can be found by entering:
#c(spcs, idExp)

# so that the directory of the experiment is now explicitly defined as:
dir_Exp <- file.path(here::here(), spcs, idExp)

# To check all constants, enter:
#set_constants_Compare()

# To check the colors, enter:
#set_colors_Compare()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                             (4)  Run the script                              #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#------------------------------------------------------------------------------#
#                                Set directories                               #
#------------------------------------------------------------------------------#

dir_C3M31 <- file.path(here::here(), spcs, "C3M31")
dir_C2M43A <- file.path(here::here(), spcs, "C2M43A")
dir_C2M43B <- file.path(here::here(), spcs, "C2M43B")
dir_C2M47 <- file.path(here::here(), spcs, "C2M47")



#------------------------------------------------------------------------------#
#         Plot the average kinetics and error interval - Non-normalized        #
#------------------------------------------------------------------------------#

# Import transpiration average
df_C3M31 <- read.csv(file.path(dir_C3M31, "Corrected_data", "E_corr_offset_avg.csv"))
df_C2M43A <- read.csv(file.path(dir_C2M43A, "Corrected_data", "E_corr_offset_avg.csv"))
df_C2M43B <- read.csv(file.path(dir_C2M43B, "Corrected_data", "E_corr_offset_avg.csv"))
df_C2M47 <- read.csv(file.path(dir_C2M47, "Corrected_data", "E_corr_offset_avg.csv"))

# Homogenize genotype names (C2M43B)
df_C2M43B$Trt[df_C2M43B$Trt == "bam1 bam3 - WW"] <- "bam1bam3 - WW"
df_C2M43B$Trt[df_C2M43B$Trt == "pgm - WW"] <- "pgm-1 - WW"

# Add idExperiment for merging
df_C3M31 <- cbind(data.frame(idExperiment = "Experiment #1"), df_C3M31)
df_C2M43A <- cbind(data.frame(idExperiment = "Experiment #2"), df_C2M43A)
df_C2M43B <- cbind(data.frame(idExperiment = "Experiment #3"), df_C2M43B)
df_C2M47 <- cbind(data.frame(idExperiment = "Experiment #4"), df_C2M47)

# Merge datasets
df <- rbind(df_C3M31, df_C2M43A, df_C2M43B, df_C2M47)

# Distinguish between 'idGenotype' and 'idWatering'
df$idGenotype <- sub("(.*) - .*", "\\1", df$Trt)
df$idWatering <- sub(".* - ", "", df$Trt)

# Create the full 'Trt' column
df$Trt <- paste(df$idExperiment, df$idGenotype, df$idWatering, sep = " - ")

# Set the colors for the experiments
ColorsTrt$col <- ColorsTrt$col_Expe

# Select the data of interest
dat <- df[df$idWatering == "WW" & df$Time <= 2, ]

# Plot the averages of the control period into a single PDF file (for visualization)
if (!dir.exists(file.path(dir_Exp, "Figures"))) { dir.create(file.path(dir_Exp, "Figures")) }
pdf(file.path(dir_Exp, "Figures", "Compare_Step#01a_plot_average_per_condition.pdf"), width = 8, height = 5)
plot_avg_kin(dat,
             c("Experiment #1 - Col-0 - WW", "Experiment #2 - Col-0 - WW", "Experiment #3 - Col-0 - WW", "Experiment #4 - Col-0 - WW"),
             use_VPD = F,
             text.legend = c("Col-0 – Experiment #1 (0.75 kPa)", "Col-0 – Experiment #2 (1 kPa)", "Col-0 – Experiment #3 (0.9 kPa)", "Col-0 – Experiment #4 (0.9 kPa)"),
             idExperiment = idExp, x.legend = -1)
plot_avg_kin(dat,
             c("Experiment #3 - pgm-1 - WW", "Experiment #4 - pgm-1 - WW"),
             use_VPD = F,
             text.legend = c(expression(paste(italic("pgm"), " – Experiment #3 (0.9 kPa)", sep = "")), expression(paste(italic("pgm"), " – Experiment #4 (0.9 kPa)", sep = ""))),
             idExperiment = idExp, x.legend = -1)
plot_avg_kin(dat,
             c("Experiment #3 - sex1-3 - WW", "Experiment #4 - sex1-3 - WW"),
             use_VPD = F,
             text.legend = c(expression(paste(italic("sex1"), " – Experiment #3 (0.9 kPa)", sep = "")), expression(paste(italic("sex1"), " – Experiment #4 (0.9 kPa)", sep = ""))),
             idExperiment = idExp, x.legend = -1)
plot_avg_kin(dat,
             c("Experiment #3 - bam1 - WW", "Experiment #4 - bam1 - WW"),
             use_VPD = F,
             text.legend = c(expression(paste(italic("bam1"), " – Experiment #3 (0.9 kPa)", sep = "")), expression(paste(italic("bam1"), " – Experiment #4 (0.9 kPa)", sep = ""))),
             idExperiment = idExp, x.legend = -1)
plot_avg_kin(dat,
             c("Experiment #3 - bam3 - WW", "Experiment #4 - bam3 - WW"),
             use_VPD = F,
             text.legend = c(expression(paste(italic("bam3"), " – Experiment #3 (0.9 kPa)", sep = "")), expression(paste(italic("bam3"), " – Experiment #4 (0.9 kPa)", sep = ""))),
             idExperiment = idExp, x.legend = -1)
plot_avg_kin(dat,
             c("Experiment #3 - bam1bam3 - WW", "Experiment #4 - bam1bam3 - WW"),
             use_VPD = F,
             text.legend = c(expression(paste(italic("bam1 bam3"), " – Experiment #3 (0.9 kPa)", sep = "")), expression(paste(italic("bam1 bam3"), " – Experiment #4 (0.9 kPa)", sep = ""))),
             idExperiment = idExp, x.legend = -1)
dev.off()

# Zoom into the first complete diel cycle of Col-0 (for publication)
if (!dir.exists(file.path(dir_Exp, "Figures", "PPTX"))) { dir.create(file.path(dir_Exp, "Figures", "PPTX")) }
dev.new(width = 3.47, height = 2.7, unit = "in", noRStudioGD = T)
plot_avg_kin(dat[dat$Time > -0.2 & dat$Time < 1.2, ],
             c("Experiment #1 - Col-0 - WW", "Experiment #2 - Col-0 - WW", "Experiment #3 - Col-0 - WW", "Experiment #4 - Col-0 - WW"),
             use_VPD = F,
             text.legend = c("#1 (0.75 kPa, n = 9)", "#2 (1 kPa, n = 6)", "#3 (0.9 kPa, n = 7)", "#4 (0.9 kPa, n = 8)"),
             xlim = c(0, 1), ylim = c(0.4, 2.8),
             idExperiment = idExp,
             draw_x_axis = T, x_axis_unit = "hours", xlab = "Time (h)", draw_rect_env = T,
             x.legend = 0.5,# y.legend = 2,
             export_PPTX = T)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin zoom_Col_Expe1_2_3_4.pptx", sep = "_")), paper = "A4", orient = "portrait", width = 3.47, height = 2.7)
dev.off()


#------------------------------------------------------------------------------#
#     Plot the average kinetics and error interval - Normalized by the VPD     #
#------------------------------------------------------------------------------#

# Import transpiration average
df_C3M31 <- read.csv(file.path(dir_C3M31, "Corrected_data", "E_corr_per_kPa_offset_avg.csv"))
df_C2M43A <- read.csv(file.path(dir_C2M43A, "Corrected_data", "E_corr_per_kPa_offset_avg.csv"))
df_C2M43B <- read.csv(file.path(dir_C2M43B, "Corrected_data", "E_corr_per_kPa_offset_avg.csv"))
df_C2M47 <- read.csv(file.path(dir_C2M47, "Corrected_data", "E_corr_per_kPa_offset_avg.csv"))

# Homogenize genotype names (C2M43B)
df_C2M43B$Trt[df_C2M43B$Trt == "bam1 bam3 - WW"] <- "bam1bam3 - WW"
df_C2M43B$Trt[df_C2M43B$Trt == "pgm - WW"] <- "pgm-1 - WW"

# Add idExperiment for merging
df_C3M31 <- cbind(data.frame(idExperiment = "Experiment #1"), df_C3M31)
df_C2M43A <- cbind(data.frame(idExperiment = "Experiment #2"), df_C2M43A)
df_C2M43B <- cbind(data.frame(idExperiment = "Experiment #3"), df_C2M43B)
df_C2M47 <- cbind(data.frame(idExperiment = "Experiment #4"), df_C2M47)

# Merge datasets
df <- rbind(df_C3M31, df_C2M43A, df_C2M43B, df_C2M47)

# Distinguish between 'idGenotype' and 'idWatering'
df$idGenotype <- sub("(.*) - .*", "\\1", df$Trt)
df$idWatering <- sub(".* - ", "", df$Trt)

# Create the full 'Trt' column
df$Trt <- paste(df$idExperiment, df$idGenotype, df$idWatering, sep = " - ")

# Set the colors for the experiments
ColorsTrt$col <- ColorsTrt$col_Expe

# Select the data of interest
dat <- df[df$idWatering == "WW" & df$Time <= 2, ]

# Plot the averages of the control period into a single PDF file (for visualization)
pdf(file.path(dir_Exp, "Figures", "Compare_Step#01b_plot_average_per_condition_use_VPD.pdf"), width = 8, height = 5)
plot_avg_kin(dat,
             c("Experiment #1 - Col-0 - WW", "Experiment #2 - Col-0 - WW", "Experiment #3 - Col-0 - WW", "Experiment #4 - Col-0 - WW"),
             use_VPD = T, add_VPD = F, add_VPD_mode = "multiple", # need to specify 'add_VPD_mode' otherwise 'add_VPD' is not passed to 'prepare_kin()' (in fact R interprets 'add_VPD' as a shortcut for 'add_VPD_mode' if only the former is specified)
             text.legend = c("Col-0 – Experiment #1 (0.75 kPa)", "Col-0 – Experiment #2 (1 kPa)", "Col-0 – Experiment #3 (0.9 kPa)", "Col-0 – Experiment #4 (0.9 kPa)"),
             idExperiment = idExp, x.legend = -1)
plot_avg_kin(dat,
             c("Experiment #3 - pgm-1 - WW", "Experiment #4 - pgm-1 - WW"),
             use_VPD = T, add_VPD = F, add_VPD_mode = "multiple", # need to specify 'add_VPD_mode' otherwise 'add_VPD' is not passed to 'prepare_kin()' (in fact R interprets 'add_VPD' as a shortcut for 'add_VPD_mode' if only the former is specified)
             text.legend = c(expression(paste(italic("pgm"), " – Experiment #3 (0.9 kPa)", sep = "")), expression(paste(italic("pgm"), " – Experiment #4 (0.9 kPa)", sep = ""))),
             idExperiment = idExp, x.legend = -1)
plot_avg_kin(dat,
             c("Experiment #3 - sex1-3 - WW", "Experiment #4 - sex1-3 - WW"),
             use_VPD = T, add_VPD = F, add_VPD_mode = "multiple", # need to specify 'add_VPD_mode' otherwise 'add_VPD' is not passed to 'prepare_kin()' (in fact R interprets 'add_VPD' as a shortcut for 'add_VPD_mode' if only the former is specified)
             text.legend = c(expression(paste(italic("sex1"), " – Experiment #3 (0.9 kPa)", sep = "")), expression(paste(italic("sex1"), " – Experiment #4 (0.9 kPa)", sep = ""))),
             idExperiment = idExp, x.legend = -1)
plot_avg_kin(dat,
             c("Experiment #3 - bam1 - WW", "Experiment #4 - bam1 - WW"),
             use_VPD = T, add_VPD = F, add_VPD_mode = "multiple", # need to specify 'add_VPD_mode' otherwise 'add_VPD' is not passed to 'prepare_kin()' (in fact R interprets 'add_VPD' as a shortcut for 'add_VPD_mode' if only the former is specified)
             text.legend = c(expression(paste(italic("bam1"), " – Experiment #3 (0.9 kPa)", sep = "")), expression(paste(italic("bam1"), " – Experiment #4 (0.9 kPa)", sep = ""))),
             idExperiment = idExp, x.legend = -1)
plot_avg_kin(dat,
             c("Experiment #3 - bam3 - WW", "Experiment #4 - bam3 - WW"),
             use_VPD = T, add_VPD = F, add_VPD_mode = "multiple", # need to specify 'add_VPD_mode' otherwise 'add_VPD' is not passed to 'prepare_kin()' (in fact R interprets 'add_VPD' as a shortcut for 'add_VPD_mode' if only the former is specified)
             text.legend = c(expression(paste(italic("bam3"), " – Experiment #3 (0.9 kPa)", sep = "")), expression(paste(italic("bam3"), " – Experiment #4 (0.9 kPa)", sep = ""))),
             idExperiment = idExp, x.legend = -1)
plot_avg_kin(dat,
             c("Experiment #3 - bam1bam3 - WW", "Experiment #4 - bam1bam3 - WW"),
             use_VPD = T, add_VPD = F, add_VPD_mode = "multiple", # need to specify 'add_VPD_mode' otherwise 'add_VPD' is not passed to 'prepare_kin()' (in fact R interprets 'add_VPD' as a shortcut for 'add_VPD_mode' if only the former is specified)
             text.legend = c(expression(paste(italic("bam1 bam3"), " – Experiment #3 (0.9 kPa)", sep = "")), expression(paste(italic("bam1 bam3"), " – Experiment #4 (0.9 kPa)", sep = ""))),
             idExperiment = idExp, x.legend = -1)
dev.off()

# Zoom into the first complete diel cycle of Col-0 (for publication)
dev.new(width = 3.47, height = 2.7, unit = "in", noRStudioGD = T)
plot_avg_kin(dat[dat$Time > -0.2 & dat$Time < 1.2, ],
             c("Experiment #1 - Col-0 - WW", "Experiment #2 - Col-0 - WW", "Experiment #3 - Col-0 - WW", "Experiment #4 - Col-0 - WW"),
             use_VPD = T, add_VPD = F, add_VPD_mode = "multiple", # need to specify 'add_VPD_mode' otherwise 'add_VPD' is not passed to 'prepare_kin()' (in fact R interprets 'add_VPD' as a shortcut for 'add_VPD_mode' if only the former is specified)
             text.legend = c("#1 (0.75 kPa, n = 9)", "#2 (1 kPa, n = 6)", "#3 (0.9 kPa, n = 7)", "#4 (0.9 kPa, n = 8)"),
             xlim = c(0, 1),# ylim = c(0, 2.5),
             idExperiment = idExp,
             draw_x_axis = T, x_axis_unit = "hours", xlab = "Time (h)", draw_rect_env = T,
             x.legend = 0.5,# y.legend = 2,
             export_PPTX = T)
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "kin zoom_Col_Expe1_2_3_4_use_VPD.pptx", sep = "_")), paper = "A4", orient = "portrait", width = 3.47, height = 2.7)
dev.off()



# Reset the colors
ColorsTrt$col <- ColorsTrt$col_Geno


#------------------------------------------------------------------------------#
#                           List the useful variables                          #
#------------------------------------------------------------------------------#

# Variables available at all periods
VAR_all <- c("E_diel_obs", "E_day_mean_obs", "E_night_mean_obs", "E_EON_obs", "daily_amp_obs",
             "A_rapid_op_obs", "A_rapid_clo_obs",
             "E_day_max_obs", "t_day_max_obs", "sigma_day_obs", "Delta_day_obs", "Sigma_preclo_obs",
             "E_night_min_obs", "t_night_min_obs", "sigma_night_obs", "Delta_night_obs", "Sigma_preop_obs",
             "E_day_mean_mod", "E_night_mean_mod", "E_EON_mod", "daily_amp_mod",
             "E_day_max_mod", "t_day_max_mod", "sigma_day_mod", "Delta_day_mod", "Sigma_preclo_mod",
             "E_night_min_mod", "t_night_min_mod", "sigma_night_mod", "Delta_night_mod", "Sigma_preop_mod")

# Variables common for acclimating and acclimated periods
VAR_common_fit <- c("E_mean", "A_SQW", "A1", "phi1", "A2", "phi2", "A1_plus_A2", "phi1_minus_phi2", "adjR2", "RMSE")

# Variables specific to  acclimating periods
VAR_acclim <- c("acclim_slope")

# Variables specific to one pot
VAR_pot <- c("diel_trend")


#------------------------------------------------------------------------------#
#           Statistical analysis of the parameters - C2M43B and C2M47          #
#------------------------------------------------------------------------------#

# Import data - no need to use VPD because both experiments have been performed at 0.9 kPa
res_C2M43B <- read.csv(file.path(dir_C2M43B, "Fitted_data", "results_fit_acclim_custom.csv"))
res_C2M47 <- read.csv(file.path(dir_C2M47, "Fitted_data", "results_fit_acclim.csv"))

# Remove one aberrant phase fit (C2M47, mex1 #349 low light)
res_C2M47[which.min(res_C2M47$phi1), c("phi1", "phi2")] <- NA

# Homogenize genotype names (C2M43B)
res_C2M43B$idGenotype[res_C2M43B$idGenotype == "bam1 bam3"] <- "bam1bam3"
res_C2M43B$idGenotype[res_C2M43B$idGenotype == "pgm"] <- "pgm-1"

# Remove low light data for pots that showed fast soil drying (C2M43B)
res_C2M43B[res_C2M43B$SWC_drop == "fast" & res_C2M43B$idPeriod %in% c("Low light acclimating", "Low light acclimated"), which(names(res_C2M43B) == "E_diel_obs"):ncol(res_C2M43B)] <- NA

# Remove 'SWC_drop' (C2M43B)
res_C2M43B <- res_C2M43B[, -which(names(res_C2M43B) == "SWC_drop")]

# # Add 'SWC_drop' (C2M47)
# res_C2M47 <- cbind(res_C2M47[1:3], data.frame(SWC_drop = "slow"), res_C2M47[4:ncol(res_C2M47)])

# Add idExperiment for merging
res_C2M43B <- cbind(data.frame(idExperiment = "Experiment #3"), res_C2M43B)
res_C2M47 <- cbind(data.frame(idExperiment = "Experiment #4"), res_C2M47)

# Merge datasets
res <- rbind(res_C2M43B, res_C2M47)

# Compute new variables
res$phi1_minus_phi2 <- res$phi1-res$phi2
res$A1_plus_A2 <- res$A1+res$A2
res$daily_amp_mod <- res$E_day_max_mod - res$E_night_min_mod
res$daily_amp_obs <- res$E_day_max_obs - res$E_night_min_obs

# Switch to hours (since the start of the day or night) where relevant
res[, c("phi1", "phi2", "phi1_minus_phi2", "t_day_max_mod", "t_day_max_obs")] <- res[, c("phi1", "phi2", "phi1_minus_phi2", "t_day_max_mod", "t_day_max_obs")]*24
res[, c("t_night_min_mod", "t_night_min_obs")] <- (res[, c("t_night_min_mod", "t_night_min_obs")] - 0.5) *24


#------------------------------------------------------------------------------#
# Exp. 1 & 2, control & low light, acclim_custom

# Plot all variables in a PDF file

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#02_stat_Expe4_3_acclim_custom.pdf", sep = "_")), width = 8, height = 5)

c(myplots1, VAR_all_common_fit) := anova_batch(res,
                                               anova_jitter_function = anova3_jitter,
                                               VAR_SELECT = c(VAR_all, VAR_common_fit),
                                               EXPE_SELECT = c("Experiment #4", "Experiment #3"),
                                               PER_SELECT = c("Control", "Low light acclimating"),
                                               GEN_SELECT = c("Col-0", "pgm-1", "sex1-3", "bam3", "bam1bam3"),
                                               IRR_SELECT = "WW",
                                               filter_phi_for_low_A = F,
                                               use_within = F,
                                               MAIN_FACTOR = "idPeriod",
                                               GROUP = "idGenotype",
                                               FACET = "idExperiment",
                                               X_LABELS = c("Control", "Low light"),
                                               LEGEND_LABELS = c("Col-0",
                                                                 expression(italic("pgm")),
                                                                 expression(italic("sex1")),
                                                                 expression(italic("bam3")),
                                                                 expression(italic("bam1 bam3"))),
                                               export_PPTX = F)
for (v in VAR_all_common_fit) { print(myplots1[[which(VAR_all_common_fit == v)]]) }


c(myplots3, VAR_acclim_pot) := anova_batch(res,
                                           anova_jitter_function = anova2_jitter,
                                           VAR_SELECT = c(VAR_acclim, VAR_pot),
                                           EXPE_SELECT = c("Experiment #4", "Experiment #3"),
                                           PER_SELECT = "Low light acclimating",
                                           GEN_SELECT = c("Col-0", "pgm-1", "sex1-3", "bam3", "bam1bam3"),
                                           IRR_SELECT = "WW",
                                           filter_phi_for_low_A = F,
                                           use_within = F,
                                           MAIN_FACTOR = "idGenotype",
                                           GROUP = "idExperiment",
                                           X_LABELS = c("Col-0",
                                                        expression(italic("pgm")),
                                                        expression(italic("sex1")),
                                                        expression(italic("bam3")),
                                                        expression(italic("bam1 bam3"))),
                                           export_PPTX = F)
for (v in VAR_acclim_pot) { print(myplots3[[which(VAR_acclim_pot == v)]]) }

dev.off()


# Export some variables as PPTX for publication

c(myplots1, VAR_all_common_fit) := anova_batch(res,
                                               anova_jitter_function = anova3_jitter,
                                               VAR_SELECT = c(VAR_all, VAR_common_fit),
                                               EXPE_SELECT = c("Experiment #4", "Experiment #3"),
                                               PER_SELECT = c("Control", "Low light acclimating"),
                                               GEN_SELECT = c("Col-0", "pgm-1", "sex1-3", "bam3", "bam1bam3"),
                                               IRR_SELECT = "WW",
                                               filter_phi_for_low_A = F,
                                               use_within = F,
                                               MAIN_FACTOR = "idPeriod",
                                               GROUP = "idGenotype",
                                               FACET = "idExperiment",
                                               LEGEND_LABELS = c("Col-0",
                                                                 expression(italic("pgm")),
                                                                 expression(italic("sex1")),
                                                                 expression(italic("bam3")),
                                                                 expression(italic("bam1 bam3"))),
                                               export_PPTX = T)

## Measured
dev.new(width = 6.5, height = 6*2/3, unit = "in", noRStudioGD = T)
ggarrange(myplots1[[which(VAR_all_common_fit == "sigma_day_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all_common_fit == "Sigma_preclo_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all_common_fit == "sigma_night_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all_common_fit == "Sigma_preop_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          ncol = 2, nrow = 2, common.legend = T, align = "v")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "_stat_Expe4_3_acclim_custom_obs.pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = 6*2/3)
dev.off()

## Fitted
dev.new(width = 6.5, height = 6, unit = "in", noRStudioGD = T)
ggarrange(myplots1[[which(VAR_all_common_fit == "E_mean")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all_common_fit == "A_SQW")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all_common_fit == "A1")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all_common_fit == "phi1")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all_common_fit == "A2")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all_common_fit == "phi2")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          ncol = 2, nrow = 3, common.legend = T, align = "v")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "_stat_Expe4_3_acclim_custom_fit.pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = 6)
dev.off()



#------------------------------------------------------------------------------#
#            Statistical analysis of the parameters - 4 experiments            #
#------------------------------------------------------------------------------#

# Import data
res_C2M43B <- read.csv(file.path(dir_C2M43B, "Fitted_data", "results_fit_use_VPD_custom.csv"))
res_C2M47 <- read.csv(file.path(dir_C2M47, "Fitted_data", "results_fit_use_VPD.csv"))
res_C2M43A <- read.csv(file.path(dir_C2M43A, "Fitted_data", "results_fit_use_VPD_custom.csv"))
res_C3M31 <- read.csv(file.path(dir_C3M31, "Fitted_data", "results_fit_use_VPD_custom.csv"))

# Remove one aberrant phase fit (C2M47, mex1 #349 low light)
res_C2M47[which.min(res_C2M47$phi1), c("phi1", "phi2")] <- NA

# Homogenize genotype names (C2M43B)
res_C2M43B$idGenotype[res_C2M43B$idGenotype == "bam1 bam3"] <- "bam1bam3"
res_C2M43B$idGenotype[res_C2M43B$idGenotype == "pgm"] <- "pgm-1"

# Remove low light data for pots that showed fast soil drying (C2M43B)
res_C2M43B[res_C2M43B$SWC_drop == "fast" & res_C2M43B$idPeriod %in% c("Low light acclimating", "Low light acclimated"), which(names(res_C2M43B) == "E_diel_obs"):ncol(res_C2M43B)] <- NA

# Remove 'SWC_drop' (C2M43B)
res_C2M43B <- res_C2M43B[, -which(names(res_C2M43B) == "SWC_drop")]

# # Remove 'acclim_slope' (C2M43B, C2M47)
# res_C2M43B <- res_C2M43B[, -which(names(res_C2M43B) == "acclim_slope")]
# res_C2M47 <- res_C2M47[, -which(names(res_C2M47) == "acclim_slope")]

# Add idExperiment for merging
res_C3M31 <- cbind(data.frame(idExperiment = "Experiment #1"), res_C3M31)
res_C2M43A <- cbind(data.frame(idExperiment = "Experiment #2"), res_C2M43A)
res_C2M43B <- cbind(data.frame(idExperiment = "Experiment #3"), res_C2M43B)
res_C2M47 <- cbind(data.frame(idExperiment = "Experiment #4"), res_C2M47)

# Merge datasets
res <- rbind(res_C2M43B, res_C2M47, res_C2M43A, res_C3M31)

# Compute new variables
res$phi1_minus_phi2 <- res$phi1-res$phi2
res$A1_plus_A2 <- res$A1+res$A2
res$daily_amp_mod <- res$E_day_max_mod - res$E_night_min_mod
res$daily_amp_obs <- res$E_day_max_obs - res$E_night_min_obs

# Switch to hours (since the start of the day or night) where relevant
res[, c("phi1", "phi2", "phi1_minus_phi2", "t_day_max_mod", "t_day_max_obs")] <- res[, c("phi1", "phi2", "phi1_minus_phi2", "t_day_max_mod", "t_day_max_obs")]*24
res[, c("t_night_min_mod", "t_night_min_obs")] <- (res[, c("t_night_min_mod", "t_night_min_obs")] - 0.5) *24


#------------------------------------------------------------------------------#
# Exp. 1 to 4, Col-0, control, use_VPD_custom

# Plot all variables in a PDF file

pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#03_stat_Expe1_2_3_4_use_VPD_custom.pdf", sep = "_")), width = 8, height = 5)

c(myplots1, VAR_all_common_fit_pot) := anova_batch(res,
                                                   anova_jitter_function = anova1_jitter,
                                                   VAR_SELECT = c(VAR_all, VAR_common_fit, VAR_pot),
                                                   EXPE_SELECT = c("Experiment #1", "Experiment #2", "Experiment #3", "Experiment #4"),
                                                   PER_SELECT = "Control",
                                                   GEN_SELECT = "Col-0",
                                                   IRR_SELECT = "WW",
                                                   filter_phi_for_low_A = F,
                                                   use_within = F,
                                                   MAIN_FACTOR = "idExperiment",
                                                   export_PPTX = F)
for (v in VAR_all_common_fit_pot) { print(myplots1[[which(VAR_all_common_fit_pot == v)]]) }

dev.off()


# Export some variables as PPTX for publication

c(myplots1, VAR_all_common_fit) := anova_batch(res,
                                               anova_jitter_function = anova1_jitter,
                                               VAR_SELECT = c(VAR_all, VAR_common_fit),
                                               EXPE_SELECT = c("Experiment #1", "Experiment #2", "Experiment #3", "Experiment #4"),
                                               PER_SELECT = "Control",
                                               GEN_SELECT = "Col-0",
                                               IRR_SELECT = "WW",
                                               filter_phi_for_low_A = F,
                                               use_within = F,
                                               MAIN_FACTOR = "idExperiment",
                                               export_PPTX = T)

dev.new(width = 6.5, height = (4-0.3)*5/3, unit = "in", noRStudioGD = T)
ggarrange(myplots1[[which(VAR_all_common_fit == "E_mean")]],
          myplots1[[which(VAR_all_common_fit == "A1")]],
          myplots1[[which(VAR_all_common_fit == "phi1")]],
          myplots1[[which(VAR_all_common_fit == "Sigma_preclo_obs")]],
          myplots1[[which(VAR_all_common_fit == "A_SQW")]],
          myplots1[[which(VAR_all_common_fit == "A2")]],
          myplots1[[which(VAR_all_common_fit == "phi2")]],
          myplots1[[which(VAR_all_common_fit == "Sigma_preop_obs")]],
          myplots1[[which(VAR_all_common_fit == "E_diel_obs")]],
          myplots1[[which(VAR_all_common_fit == "daily_amp_obs")]],
          myplots1[[which(VAR_all_common_fit == "A_rapid_op_obs")]],
          myplots1[[which(VAR_all_common_fit == "A_rapid_clo_obs")]],
          myplots1[[which(VAR_all_common_fit == "E_day_mean_obs")]],
          myplots1[[which(VAR_all_common_fit == "t_day_max_obs")]],
          myplots1[[which(VAR_all_common_fit == "sigma_day_obs")]],
          myplots1[[which(VAR_all_common_fit == "Delta_day_obs")]],
          myplots1[[which(VAR_all_common_fit == "E_night_mean_obs")]],
          myplots1[[which(VAR_all_common_fit == "t_night_min_obs")]],
          myplots1[[which(VAR_all_common_fit == "sigma_night_obs")]],
          myplots1[[which(VAR_all_common_fit == "Delta_night_obs")]],
          ncol = 4, nrow = 5, common.legend = T, align = "v")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "_stat_Expe1_2_3_4_use_VPD_custom.pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = (4-0.3)*5/3)
dev.off()
