################################################################################
#                                                                              #
#                 PhenoLeaks - Step #05 - Statistical analyses                 #
#                                                                              #
#                    Script to perform statistical analyses                    #
#                on fitted and observed transpiration parameters               #
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
for (pkg in c("gplots", "here"))
  {
  if (!pkg %in% installed.packages()[, "Package"]) { install.packages(pkg) }
  #update.packages(pkg)
  library(pkg, character.only = T)
  }

if (!"ComplexHeatmap" %in% installed.packages()[, "Package"])
  {
  if (!"BiocManager" %in% installed.packages()[, "Package"]) { install.packages("BiocManager") }
  library("BiocManager")
  BiocManager::install("ComplexHeatmap")
  }
library("ComplexHeatmap")



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
#                         Import and manage fitted data                        #
#------------------------------------------------------------------------------#

# Import data
res <- read.csv(file.path(dir_Exp, "Fitted_data", "results_fit_acclim.csv"))
res <- cbind(data.frame(idExperiment = idExp), res)

# Back to 'idPeriod' former name (for compatibility)
res$idPeriod[res$idPeriod == "High CO2 acclimating"] <- "High CO2"
res$idPeriod[res$idPeriod == "Low light acclimating"] <- "Low light"
res$idPeriod[res$idPeriod == "Recovery acclimating"] <- "Recovery 1"
res$idPeriod[res$idPeriod == "Recovery acclimated"] <- "Recovery 2"

# Compute new variables
res$phi1_minus_phi2 <- res$phi1-res$phi2
res$A1_plus_A2 <- res$A1+res$A2
res$daily_amp_mod <- res$E_day_max_mod - res$E_night_min_mod
res$daily_amp_obs <- res$E_day_max_obs - res$E_night_min_obs

# Switch to hours (since the start of the day or night) where relevant
res[, c("phi1", "phi2", "phi1_minus_phi2", "t_day_max_mod", "t_day_max_obs")] <- res[, c("phi1", "phi2", "phi1_minus_phi2", "t_day_max_mod", "t_day_max_obs")]*24
res[, c("t_night_min_mod", "t_night_min_obs")] <- (res[, c("t_night_min_mod", "t_night_min_obs")] - 0.5) *24

# Remove one aberrant phase fit (mex1 #349 low light)
hist(res$phi1)
hist(res$phi1[res$phi1 <= (-1)])
hist(res$phi1[res$phi1 <= (-2)]) # only one observation below -4h for phi1, reaching -11h...
res[which.min(res$phi1), ]
res[which.min(res$phi1), c("phi1", "phi2", "phi1_minus_phi2")] <- NA


#------------------------------------------------------------------------------#
#               Effect of genotypes across periods (well-watered)              #
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
#0# Set generic functions

# Plot all variables in a PDF file
plot_geno_PDF <- function (GEN_SELECT, GENO_LABELS, gen_file,
                           IRR_SELECT = "WW",
                           MAIN_FACTOR = "idPeriod",
                           GROUP = "idGenotype",
                           PERIOD_LABELS = c("Control", "High CO2", "Low light", "Recovery 1", "Recovery 2"),
                           PERIOD_FIT_LABELS = c("Control", "High CO2", "Low light", "Recovery"),
                           PERIOD_ACCLIM_LABELS = c("High CO2", "Low light", "Recovery 1"))
  {
  pdf(file.path(dir_Exp, "Figures", paste(idExp, "_Step#05a_stat_", gen_file, ".pdf", sep = "")), width = 8, height = 5)
  
  c(myplots1, VAR_all) := anova_batch(res,
                                      anova_jitter_function = anova2_jitter,
                                      VAR_SELECT = VAR_all,
                                      EXPE_SELECT = idExp,
                                      PER_SELECT = c("Control", "High CO2", "Low light", "Recovery 1", "Recovery 2"),
                                      GEN_SELECT = GEN_SELECT,
                                      IRR_SELECT = IRR_SELECT,
                                      filter_phi_for_low_A = F,
                                      use_within = F,
                                      MAIN_FACTOR = MAIN_FACTOR,
                                      GROUP = GROUP,
                                      LEGEND_LABELS = if (MAIN_FACTOR == "idPeriod") GENO_LABELS else PERIOD_LABELS,
                                      X_LABELS = if (MAIN_FACTOR == "idPeriod") PERIOD_LABELS else GENO_LABELS, 
                                      export_PPTX = F)
  for (v in VAR_all) { print(myplots1[[which(VAR_all == v)]]) }
  
  c(myplots2, VAR_common_fit) := anova_batch(res,
                                             anova_jitter_function = anova2_jitter,
                                             VAR_SELECT = VAR_common_fit,
                                             EXPE_SELECT = idExp,
                                             PER_SELECT = c("Control", "High CO2", "Low light", "Recovery 1"),
                                             GEN_SELECT = GEN_SELECT,
                                             IRR_SELECT = IRR_SELECT,
                                             filter_phi_for_low_A = F,
                                             use_within = F,
                                             MAIN_FACTOR = MAIN_FACTOR,
                                             GROUP = GROUP,
                                             LEGEND_LABELS = if (MAIN_FACTOR == "idPeriod") GENO_LABELS else PERIOD_FIT_LABELS,
                                             X_LABELS = if (MAIN_FACTOR == "idPeriod") PERIOD_FIT_LABELS else GENO_LABELS, 
                                             export_PPTX = F)
  for (v in VAR_common_fit) { print(myplots2[[which(VAR_common_fit == v)]]) }

  c(myplots3, VAR_acclim) := anova_batch(res,
                                         anova_jitter_function = anova2_jitter,
                                         VAR_SELECT = VAR_acclim,
                                         EXPE_SELECT = idExp,
                                         PER_SELECT = c("High CO2", "Low light", "Recovery 1"),
                                         GEN_SELECT = GEN_SELECT,
                                         IRR_SELECT = IRR_SELECT,
                                         filter_phi_for_low_A = F,
                                         use_within = F,
                                         MAIN_FACTOR = MAIN_FACTOR,
                                         GROUP = GROUP,
                                         LEGEND_LABELS = if (MAIN_FACTOR == "idPeriod") GENO_LABELS else PERIOD_ACCLIM_LABELS,
                                         X_LABELS = if (MAIN_FACTOR == "idPeriod") PERIOD_ACCLIM_LABELS else GENO_LABELS, 
                                         export_PPTX = F)
  for (v in VAR_acclim) { print(myplots3[[which(VAR_acclim == v)]]) }
  
  if (MAIN_FACTOR == "idPeriod")
    {
    c(myplots4, VAR_pot) := anova_batch(res,
                                        anova_jitter_function = anova1_jitter,
                                        VAR_SELECT = VAR_pot,
                                        EXPE_SELECT = idExp,
                                        PER_SELECT = "Control",
                                        GEN_SELECT = GEN_SELECT,
                                        IRR_SELECT = IRR_SELECT,
                                        filter_phi_for_low_A = F,
                                        use_within = F,
                                        MAIN_FACTOR = "idGenotype",
                                        LEGEND_LABELS = GENO_LABELS,
                                        X_LABELS = GENO_LABELS,
                                        export_PPTX = F)
    for (v in VAR_pot) { print(myplots4[[which(VAR_pot == v)]]) }
    }
  
  dev.off()
  }


# Plot all variables before export in a PPTX file
plot_geno_PPTX <- function (GEN_SELECT, GENO_LABELS)
  {
  c(myplots1, VAR_all) := anova_batch(res,
                                      anova_jitter_function = anova2_jitter,
                                      VAR_SELECT = VAR_all,
                                      EXPE_SELECT = idExp,
                                      PER_SELECT = c("Control", "High CO2", "Low light", "Recovery 1", "Recovery 2"),
                                      GEN_SELECT = GEN_SELECT,
                                      IRR_SELECT = "WW",
                                      filter_phi_for_low_A = F,
                                      use_within = F,
                                      MAIN_FACTOR = "idPeriod",
                                      GROUP = "idGenotype",
                                      LEGEND_LABELS = GENO_LABELS,
                                      export_PPTX = T)
  
  c(myplots2, VAR_common_fit) := anova_batch(res,
                                             anova_jitter_function = anova2_jitter,
                                             VAR_SELECT = VAR_common_fit,
                                             EXPE_SELECT = idExp,
                                             PER_SELECT = c("Control", "High CO2", "Low light", "Recovery 1"),
                                             GEN_SELECT = GEN_SELECT,
                                             IRR_SELECT = "WW",
                                             filter_phi_for_low_A = F,
                                             use_within = F,
                                             MAIN_FACTOR = "idPeriod",
                                             GROUP = "idGenotype",
                                             LEGEND_LABELS = GENO_LABELS,
                                             export_PPTX = T)

  c(myplots3, VAR_acclim) := anova_batch(res,
                                         anova_jitter_function = anova2_jitter,
                                         VAR_SELECT = VAR_acclim,
                                         EXPE_SELECT = idExp,
                                         PER_SELECT = c("High CO2", "Low light", "Recovery 1"),
                                         GEN_SELECT = GEN_SELECT,
                                         IRR_SELECT = "WW",
                                         filter_phi_for_low_A = F,
                                         use_within = F,
                                         MAIN_FACTOR = "idPeriod",
                                         GROUP = "idGenotype",
                                         LEGEND_LABELS = GENO_LABELS,
                                         export_PPTX = T)
  
  c(myplots4, VAR_pot) := anova_batch(res,
                                      anova_jitter_function = anova1_jitter,
                                      VAR_SELECT = VAR_pot,
                                      EXPE_SELECT = idExp,
                                      PER_SELECT = "Recovery 1",
                                      GEN_SELECT = GEN_SELECT,
                                      IRR_SELECT = "WW",
                                      filter_phi_for_low_A = F,
                                      use_within = F,
                                      MAIN_FACTOR = "idGenotype",
                                      LEGEND_LABELS = GENO_LABELS,
                                      export_PPTX = T)

  return (list(myplots1, myplots2, myplots3, myplots4))
  }


# Grids of plots for PPTX export (not used for pgm and sex1, which have a slightly different layout)

PPTX_stat_main <- function (gen_file)
  {
  dev.new(width = 6.5, height = (4-0.3)*2/3, unit = "in", noRStudioGD = T)
  myfig <- ggarrange(myplots2[[which(VAR_common_fit == "E_mean")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
                     myplots2[[which(VAR_common_fit == "A1")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
                     myplots2[[which(VAR_common_fit == "phi1")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
                     myplots1[[which(VAR_all == "Sigma_preclo_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
                     myplots2[[which(VAR_common_fit == "A_SQW")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
                     myplots2[[which(VAR_common_fit == "A2")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
                     myplots2[[which(VAR_common_fit == "phi2")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
                     myplots1[[which(VAR_all == "Sigma_preop_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
                     ncol = 4, nrow = 2, common.legend = T, align = "v")
  print(myfig)
  graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "_stat_main_", gen_file, ".pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = (4-0.3)*2/3)
  dev.off()
  }

PPTX_stat_suppl <- function (gen_file)
  {
  dev.new(width = 6.5, height = 4-0.3, unit = "in", noRStudioGD = T)
  myfig <- ggarrange(myplots1[[which(VAR_all == "E_diel_obs")]],
                     myplots1[[which(VAR_all == "daily_amp_obs")]],
                     myplots1[[which(VAR_all == "A_rapid_op_obs")]],
                     myplots1[[which(VAR_all == "A_rapid_clo_obs")]],
                     myplots1[[which(VAR_all == "E_day_mean_obs")]],
                     myplots1[[which(VAR_all == "t_day_max_obs")]],
                     myplots1[[which(VAR_all == "sigma_day_obs")]],
                     myplots1[[which(VAR_all == "Delta_day_obs")]],
                     myplots1[[which(VAR_all == "E_night_mean_obs")]],
                     myplots1[[which(VAR_all == "t_night_min_obs")]],
                     myplots1[[which(VAR_all == "sigma_night_obs")]],
                     myplots1[[which(VAR_all == "Delta_night_obs")]],
                     ncol = 4, nrow = 3, common.legend = T, align = "v")
  print(myfig)
  graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "_stat_suppl_", gen_file, ".pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = 4-0.3)
  dev.off()
  }


#------------------------------------------------------------------------------#
#1# pgm and sex1

GEN_SELECT <- c("Col-0",
                "pgm-1",
                "sex1-3")
GENO_LABELS <- c("Col-0",
                   expression(italic("pgm")),
                   expression(italic("sex1")))
gen_file <- "geno_pgm_sex1"

# Get the mean and SD
#dat <- res[res$idGenotype %in% GEN_SELECT & res$idWatering == "WW", ]
#VAR_SELECT <- "phi1" #VAR_all_common_fit
#aggregate(list(dat[, VAR_SELECT]), by = list(idGenotype = dat$idGenotype, idPeriod = dat$idPeriod), FUN = mean, na.rm = T)
#aggregate(list(dat[, VAR_SELECT]), by = list(idGenotype = dat$idGenotype, idPeriod = dat$idPeriod), FUN = function (x) { sd(x, na.rm = T) / sqrt(length(na.omit(x)) - 1) } )

# Plot all variables in a PDF file - STATS PER GENOTYPE (default)
plot_geno_PDF(GEN_SELECT, GENO_LABELS, gen_file)

# Plot all variables in a PDF file - STATS PER PERIOD
# The effect sizes are the same as before, but the pairwise t-tests are performed across periods
plot_geno_PDF(GEN_SELECT, GENO_LABELS,
              gen_file = "period_pgm_sex1", MAIN_FACTOR = "idGenotype", GROUP = "idPeriod")

# Export some variables as PPTX for publication
c(myplots1, myplots2, myplots3, myplots4) := plot_geno_PPTX(GEN_SELECT, GENO_LABELS)

dev.new(width = 6.5/2, height = 4-0.3, unit = "in", noRStudioGD = T)
ggarrange(myplots1[[which(VAR_all == "A_rapid_op_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all == "A_rapid_clo_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all == "t_day_max_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all == "t_night_min_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all == "Sigma_preclo_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all == "Sigma_preop_obs")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          ncol = 2, nrow = 3, common.legend = T, align = "v")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "_stat_main_", gen_file, "_obs.pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5/2, height = 4-0.3)
dev.off()

dev.new(width = 6.5, height = (4-0.3)*2/3, unit = "in", noRStudioGD = T)
ggarrange(myplots1[[which(VAR_all == "E_diel_obs")]],
          myplots1[[which(VAR_all == "E_day_mean_obs")]],
          myplots1[[which(VAR_all == "sigma_day_obs")]],
          myplots1[[which(VAR_all == "Delta_day_obs")]],
          myplots1[[which(VAR_all == "daily_amp_obs")]],
          myplots1[[which(VAR_all == "E_night_mean_obs")]],
          myplots1[[which(VAR_all == "sigma_night_obs")]],
          myplots1[[which(VAR_all == "Delta_night_obs")]],
          ncol = 4, nrow = 2, common.legend = T, align = "v")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "_stat_suppl_", gen_file, ".pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = (4-0.3)*2/3)
dev.off()

dev.new(width = 6.5, height = (4-0.3)*2/3, unit = "in", noRStudioGD = T)
ggarrange(myplots2[[which(VAR_common_fit == "E_mean")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots2[[which(VAR_common_fit == "A1")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots2[[which(VAR_common_fit == "phi1")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all == "Sigma_preclo_mod")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots2[[which(VAR_common_fit == "A_SQW")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots2[[which(VAR_common_fit == "A2")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots2[[which(VAR_common_fit == "phi2")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots1[[which(VAR_all == "Sigma_preop_mod")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          ncol = 4, nrow = 2, common.legend = T, align = "v")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "_stat_main_", gen_file, "_fit.pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = (4-0.3)*2/3)
dev.off()


#------------------------------------------------------------------------------#
#2# isa1 and ss4

GEN_SELECT <- c("Col-0 ProStarv::LUC",
                "isa1-1 ProStarv::LUC",
                "ss4-3 ProStarv::LUC")
GENO_LABELS <- c(expression("Col-0"^"†"),
                   expression(italic("isa1")^"†"),
                   expression(italic("ss4")^"†"))
gen_file <- "geno_isa1_ss4"


# Plot all variables in a PDF file
plot_geno_PDF(GEN_SELECT, GENO_LABELS, gen_file)

# Export some variables as PPTX for publication
c(myplots1, myplots2, myplots3, myplots4) := plot_geno_PPTX(GEN_SELECT, GENO_LABELS)
PPTX_stat_main(gen_file)
PPTX_stat_suppl(gen_file)  


#------------------------------------------------------------------------------#
#3# mex1, dpe1 and dpe2

GEN_SELECT <- c("Col-0",
                "mex1-1",
                "dpe1-2",
                "dpe2-5")
GENO_LABELS <- c("Col-0",
                   expression(italic("mex1")),
                   expression(italic("dpe1")),
                   expression(italic("dpe2")))
gen_file <- "geno_mex1_dpe1_dpe2"

# Plot all variables in a PDF file
plot_geno_PDF(GEN_SELECT, GENO_LABELS, gen_file)

# Export some variables as PPTX for publication
c(myplots1, myplots2, myplots3, myplots4) := plot_geno_PPTX(GEN_SELECT, GENO_LABELS)
PPTX_stat_main(gen_file)
PPTX_stat_suppl(gen_file)  


#------------------------------------------------------------------------------#
#4# pgi, amy3 bam1 and bam1 bam3

GEN_SELECT <- c("Col-0",
                "pgi1-1",
                "amy3-2 bam1",
                "bam1bam3")
GENO_LABELS <- c("Col-0",
                   expression(italic("pgi")),
                   expression(italic("amy3 bam1")),
                   expression(italic("bam1 bam3")))
gen_file <- "geno_pgi_amy3bam1_bam1bam3"

# Plot all variables in a PDF file
plot_geno_PDF(GEN_SELECT, GENO_LABELS, gen_file)

# Export some variables as PPTX for publication
c(myplots1, myplots2, myplots3, myplots4) := plot_geno_PPTX(GEN_SELECT, GENO_LABELS)
PPTX_stat_main(gen_file)
PPTX_stat_suppl(gen_file)  


#------------------------------------------------------------------------------#
#5# amy3, bam1 and bam3

GEN_SELECT <- c("Col-0",
                "amy3-2",
                "bam1",
                "bam3")
GENO_LABELS <- c("Col-0",
                   expression(italic("amy3")),
                   expression(italic("bam1")),
                   expression(italic("bam3")))
gen_file <- "geno_amy3_bam1_bam3"

# Plot all variables in a PDF file
plot_geno_PDF(GEN_SELECT, GENO_LABELS, gen_file)

# Export some variables as PPTX for publication
c(myplots1, myplots2, myplots3, myplots4) := plot_geno_PPTX(GEN_SELECT, GENO_LABELS)
PPTX_stat_main(gen_file)
PPTX_stat_suppl(gen_file)  


#------------------------------------------------------------------------------#
#6# abcb14-1 and -2

GEN_SELECT <- c("Col-0",
                "abcb14-1",
                "abcb14-2")
GENO_LABELS <- c("Col-0",
                   expression(paste(italic("abcb14"), "-1", sep = "")),
                   expression(paste(italic("abcb14"), "-2", sep = "")))
gen_file <- "geno_abcb14"

# Plot all variables in a PDF file
plot_geno_PDF(GEN_SELECT, GENO_LABELS, gen_file)

# Export some variables as PPTX for publication
c(myplots1, myplots2, myplots3, myplots4) := plot_geno_PPTX(GEN_SELECT, GENO_LABELS)
PPTX_stat_main(gen_file)
PPTX_stat_suppl(gen_file)


#------------------------------------------------------------------------------#
#7# Acclimation slope - All genotypes and periods

GEN_SELECT <-c("Col-0",
               "Col-0 ProStarv::LUC",
               "pgi1-1",
               "pgm-1",
               "isa1-1 ProStarv::LUC",
               "ss4-3 ProStarv::LUC",
               "sex1-3",
               "amy3-2",
               "amy3-2 bam1",
               "bam1",
               "bam3",
               "bam1bam3",
               "mex1-1",
               "dpe1-2",
               "dpe2-5",
               "abcb14-1",
               "abcb14-2")

GENO_LABELS <- c(expression("Col-0"),
                 expression("Col-0"^"†"),
                 expression(italic("pgi")),
                 expression(italic("pgm")),               
                 expression(italic("isa1")^"†"),
                 expression(italic("ss4")^"†"),                
                 expression(italic("sex1")),
                 expression(italic("amy3")),
                 expression(italic("amy3 bam1")),
                 expression(italic("bam1")),
                 expression(italic("bam3")),
                 expression(italic("bam1 bam3")),
                 expression(italic("mex1")),
                 expression(italic("dpe1")),
                 expression(italic("dpe2")),                
                 expression(paste(italic("abcb14"), "-1", sep = "")),
                 expression(paste(italic("abcb14"), "-2", sep = "")))     

# Note: manually change step.increase = 0.001 in add_xy_position() of anova2_jitter() before running the command below, and then revert
c(myplots_acclim, VAR_acclim) := anova_batch(res,
                                       anova_jitter_function = anova2_jitter,
                                       VAR_SELECT = VAR_acclim,
                                       EXPE_SELECT = idExp,
                                       PER_SELECT = c("High CO2", "Low light", "Recovery 1"),
                                       GEN_SELECT = GEN_SELECT,
                                       IRR_SELECT = "WW",
                                       filter_phi_for_low_A = F,
                                       use_within = F,
                                       MAIN_FACTOR = "idPeriod",
                                       GROUP = "idGenotype",
                                       LEGEND_LABELS = GENO_LABELS,
                                       export_PPTX = T)

# The code is commented to secure the pptx file (manual step above)
#dev.new(width = 6.5, height = 4.5, unit = "in", noRStudioGD = T)
#myplots_acclim 
#graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "_stat_acclim.pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = 4.5)
#dev.off()


#------------------------------------------------------------------------------#
#         Effect of irrigation regimes across periods (Col-0 and pgm)          #
#------------------------------------------------------------------------------#

# Plot all variables in a PDF file
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#05b_stat_irrigation_Col_pgm.pdf", sep = "_")), width = 8, height = 5)
  
  c(myplots1, VAR_all) := anova_batch(res,
                                      anova_jitter_function = anova3_jitter,
                                      VAR_SELECT = VAR_all,
                                      EXPE_SELECT = idExp,
                                      PER_SELECT = c("Control", "High CO2", "Low light", "Recovery 1", "Recovery 2"),
                                      GEN_SELECT = c("Col-0", "pgm-1"),
                                      IRR_SELECT = c("WW", "WS"),
                                      filter_phi_for_low_A = F,
                                      use_within = F,
                                      MAIN_FACTOR = "idPeriod",
                                      GROUP = "idGenotype",
                                      FACET = "idWatering",
                                      LEGEND_LABELS = c("Col-0", expression(italic("pgm"))),
                                      export_PPTX = F)
  for (v in VAR_all) { print(myplots1[[which(VAR_all == v)]]) }
  
  c(myplots2, VAR_common_fit) := anova_batch(res,
                                             anova_jitter_function = anova3_jitter,
                                             VAR_SELECT = VAR_common_fit,
                                             EXPE_SELECT = idExp,
                                             PER_SELECT = c("Control", "High CO2", "Low light", "Recovery 1"),
                                             GEN_SELECT = c("Col-0", "pgm-1"),
                                             IRR_SELECT = c("WW", "WS"),
                                             filter_phi_for_low_A = F,
                                             use_within = F,
                                             MAIN_FACTOR = "idPeriod",
                                             GROUP = "idGenotype",
                                             FACET = "idWatering",
                                             LEGEND_LABELS = c("Col-0", expression(italic("pgm"))),
                                             X_LABELS = c("Control", "High CO2", "Low light", "Recovery"),
                                             export_PPTX = F)
  for (v in VAR_common_fit) { print(myplots2[[which(VAR_common_fit == v)]]) }
  
  c(myplots3, VAR_acclim) := anova_batch(res,
                                         anova_jitter_function = anova3_jitter,
                                         VAR_SELECT = VAR_acclim,
                                         EXPE_SELECT = idExp,
                                         PER_SELECT = c("High CO2", "Low light", "Recovery 1"),
                                         GEN_SELECT = c("Col-0", "pgm-1"),
                                         IRR_SELECT = c("WW", "WS"),
                                         filter_phi_for_low_A = F,
                                         use_within = F,
                                         MAIN_FACTOR = "idPeriod",
                                         GROUP = "idGenotype",
                                         FACET = "idWatering",
                                         LEGEND_LABELS = c("Col-0", expression(italic("pgm"))),
                                         export_PPTX = F)
  for (v in VAR_acclim) { print(myplots3[[which(VAR_acclim == v)]]) }
  
  c(myplots4, VAR_pot) := anova_batch(res,
                                      anova_jitter_function = anova2_jitter,
                                      VAR_SELECT = VAR_pot,
                                      EXPE_SELECT = idExp,
                                      PER_SELECT = "Control",
                                      GEN_SELECT = c("Col-0", "pgm-1"),
                                      IRR_SELECT = c("WW", "WS"),
                                      filter_phi_for_low_A = F,
                                      use_within = F,
                                      MAIN_FACTOR = "idGenotype",
                                      GROUP = "idWatering",
                                      X_LABELS = c("Col-0", expression(italic("pgm"))),
                                      export_PPTX = F)
  for (v in VAR_pot) { print(myplots4[[which(VAR_pot == v)]]) }
  
dev.off()


# Export some variables as PPTX for publication
res$idWatering[res$idWatering == "WW"] <- "Well-watered"
res$idWatering[res$idWatering == "WS"] <- "Water stress"

c(myplots2, VAR_common_fit) := anova_batch(res,
                                           anova_jitter_function = anova3_jitter,
                                           VAR_SELECT = VAR_common_fit,
                                           EXPE_SELECT = idExp,
                                           PER_SELECT = c("Control", "High CO2", "Low light", "Recovery 1"),
                                           GEN_SELECT = c("Col-0", "pgm-1"),
                                           IRR_SELECT = c("Well-watered", "Water stress"),
                                           filter_phi_for_low_A = F,
                                           use_within = F,
                                           MAIN_FACTOR = "idPeriod",
                                           GROUP = "idGenotype",
                                           FACET = "idWatering",
                                           LEGEND_LABELS = c("Col-0", expression(italic("pgm"))),
                                           export_PPTX = T)

dev.new(width = 6.5, height = 6, unit = "in", noRStudioGD = T)
ggarrange(myplots2[[which(VAR_common_fit == "E_mean")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots2[[which(VAR_common_fit == "A_SQW")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots2[[which(VAR_common_fit == "A1")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots2[[which(VAR_common_fit == "phi1")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots2[[which(VAR_common_fit == "A2")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          myplots2[[which(VAR_common_fit == "phi2")]] + theme(plot.margin = unit(c(0,1,0,13), "points")),
          ncol = 2, nrow = 3, common.legend = T, align = "v")
graph2ppt(file = file.path(dir_Exp, "Figures", "PPTX", paste(idExp, "_stat_irrigation_Col_pgm.pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = 6)
dev.off()

res$idWatering[res$idWatering == "Well-watered"] <- "WW"
res$idWatering[res$idWatering == "Water stress"] <- "WS"




#------------------------------------------------------------------------------#
#                                    Heatmap                                   #
#------------------------------------------------------------------------------#

#require(ComplexHeatmap)
#require(gplots)

#------------------------------------------------------------------------------#
#1# With acclimation

# Filter the phases for low amplitudes (uncomment if required)
res_dupl <- res
#res_dupl[res_dupl$A1 < 0.05 & !is.na(res_dupl$A1), c("phi1", "phi1_minus_phi2", "t_day_max_mod", "t_day_max_obs", "t_night_min_mod", "t_night_min_obs")] <- NA
#res_dupl[res_dupl$A2 < 0.02 & !is.na(res_dupl$A2), c("phi2", "phi1_minus_phi2")] <- NA


# Note that "Recovery 2" (acclimated) has the same fitted parameters as "Recovery 1" (acclimating) except for the acclimation slope,
# which has been set to 0 for "Recovery 2" by the function 'run_fit()' of 'PhenoLeaks_fit.R'.
# We change it to NA:
res_dupl$acclim_slope[res_dupl$idPeriod == "Recovery 2"] <- NA

# Create the average matrix
avg <- aggregate(res_dupl[, !colnames(res_dupl) %in% c("idExperiment", "idGenotype", "idWatering", "idPot", "idPeriod", "diel_trend", "diel_trend_baseline", "t1_fit", "t2_fit", "adjR2", "RMSE", "E_EON_before_obs", "E_EON_before_mod", "A_rapid_op_transition_obs", "A_rapid_op_stable_obs",
                                                       "E_day_mean_mod", "E_day_mean_obs", "E_night_mean_mod", "E_night_mean_obs", "E_EON_mod", "E_EON_obs", "E_day_max_obs", "E_day_max_mod", "E_night_min_obs", "E_night_min_mod", "phi1_minus_phi2", "A1_plus_A2")],
                 by = list (idWatering = res_dupl$idWatering, idGenotype = res_dupl$idGenotype, idPeriod = res_dupl$idPeriod),
                 FUN = mean, na.rm = T)
mat <- as.matrix(avg[-(1:3)])
rownames(mat) <- paste(avg[,1], avg[,2], avg[,3])


# Set vector of row names
LABROW <- c(expression(paste(italic("abcb14"), "-1", sep = "")),
            expression(paste(italic("abcb14"), "-2", sep = "")),
            expression(italic("amy3")),
            expression(italic("amy3 bam1")),
            expression(italic("bam1")),
            expression(italic("bam1 bam3")),
            expression(italic("bam3")),
            expression("Col-0 (WS)"),
            expression("Col-0"),
            expression("Col-0"^"†"),
            expression(italic("dpe1")),
            expression(italic("dpe2")),
            expression(italic("isa1")^"†"),
            expression(italic("mex1")),
            expression(italic("pgi")),
            expression(paste(italic("pgm"), " (WS)", sep = "")),
            expression(italic("pgm")),
            expression(italic("sex1")),
            expression(italic("ss4")^"†"))

# Set vector of column names
LABCOL <- set_names_units("names")[match(colnames(mat), set_names_units("par"))] 

# Set the colors for the traits
trait_colors <- c(palette("default")[1], palette("default")[3], palette("default")[3],
                  "darkorange", "darkorange3", "darkorange3", "darkorange",
                  "deepskyblue2", "dodgerblue4", "dodgerblue4", "deepskyblue2",
                  "aquamarine4", palette("default")[3], "coral3", "coral3", "coral1", "coral1", "coral4",
                  "darkorange", "darkorange3", "darkorange3", "darkorange",
                  "deepskyblue2", "dodgerblue4", "dodgerblue4", "deepskyblue2",
                  "palegreen3", "palegreen3")

# Set the colors for the periods
vector_colperiod <- ColorsPeriod$col
names(vector_colperiod) <- ColorsPeriod$idPeriod


#-----------------------------#
#1.1# Well-watered only

mat_WW <- mat[avg$idWatering == "WW", ]
mat_WW <- scale(mat_WW)

LABROW_WW <- LABROW[!grepl("WS", LABROW)]

row_ha <- rowAnnotation(Environment = avg$idPeriod[avg$idWatering == "WW"],
                        col = list(Environment = vector_colperiod),
                        show_annotation_name = F,
                        annotation_legend_param = list(Environment = list(labels = c("Control", expression("High CO"[2]), "Low light", "Recov.1", "Recov.2"),
                                                                          nrow = 1,
                                                                          title_position = "topleft",
                                                                          labels_gp = gpar(fontsize = 8),
                                                                          title_gp = gpar(fontsize = 10, fontface = "plain"))))

hmp_acclim_WW <- Heatmap(mat_WW, name = "Trait value (scaled and centered)",
                         col = colorpanel(256, "purple", "black", "orange"),
                         #rect_gp = gpar(col = "gray", lwd = 0.01),
                         row_labels = rep(LABROW_WW, 5),
                         column_labels = LABCOL,
                         row_names_gp = gpar(fontsize = 7, col = ColorsTrt$col[match(avg$idGenotype[avg$idWatering == "WW"], ColorsTrt$idGenotype)]),
                         column_names_gp = gpar(fontsize = 8, col = trait_colors),
                         column_names_rot = 45,
                         left_annotation = row_ha,
                         heatmap_legend_param = list(direction = "horizontal",
                                                     title_position = "topcenter",
                                                     legend_width = unit(6, "cm"),
                                                     title_gp = gpar(fontsize = 10, fontface = "plain"),
                                                     labels_gp = gpar(fontsize = 8)),
                         row_split = 4, column_split = 4,
                         #row_gap = unit(c(2, 4, 2), "pt"), column_gap = unit(c(2, 4, 2), "pt"),
                         row_title = NULL, column_title = "Model with acclimation - well-watered only")

#dev.new(width = 6.5, height = 7.5, unit = "in", noRStudioGD = T)
#draw(hmp_acclim_WW, heatmap_legend_side = "bottom")


#-----------------------------#
#1.2# Including water stress

mat_all <- scale(mat)

row_ha <- rowAnnotation(Environment = avg$idPeriod,
                        col = list(Environment = vector_colperiod),
                        show_annotation_name = F,
                        annotation_legend_param = list(Environment = list(labels = c("Control", expression("High CO"[2]), "Low light", "Recov.1", "Recov.2"),
                                                                          nrow = 1,
                                                                          title_position = "topleft",
                                                                          labels_gp = gpar(fontsize = 8),
                                                                          title_gp = gpar(fontsize = 10, fontface = "plain"))))

hmp_acclim <- Heatmap(mat_all, name = "Trait value (scaled and centered)",
                      col = colorpanel(256, "purple", "black", "orange"),
                      #rect_gp = gpar(col = "gray", lwd = 0.01),
                      row_labels = rep(LABROW, 5),
                      column_labels = LABCOL,
                      row_names_gp = gpar(fontsize = 7, col = ColorsTrt$col[match(avg$idGenotype, ColorsTrt$idGenotype)]),
                      column_names_gp = gpar(fontsize = 8, col = trait_colors),
                      column_names_rot = 45,
                      left_annotation = row_ha,
                      heatmap_legend_param = list(direction = "horizontal",
                                                  title_position = "topcenter",
                                                  legend_width = unit(6, "cm"),
                                                  title_gp = gpar(fontsize = 10, fontface = "plain"),
                                                  labels_gp = gpar(fontsize = 8)),
                      row_split = 4, column_split = 4,
                      #row_gap = unit(c(2, 4, 2), "pt"), column_gap = unit(c(2, 4, 2), "pt"),
                      row_title = NULL, column_title = "Model with acclimation - including water stress (WS)")

#dev.new(width = 6.5, height = 7.5, unit = "in", noRStudioGD = T)
#draw(hmp_acclim, heatmap_legend_side = "bottom")


#------------------------------------------------------------------------------#
#2# Without acclimation

# Import data without acclimation
res <- read.csv(file.path(dir_Exp, "Fitted_data", "results_fit.csv"))
res <- cbind(data.frame(idExperiment = idExp), res)

# Back to 'idPeriod' former name (for compatibility)
res$idPeriod[res$idPeriod == "High CO2 acclimating"] <- "High CO2"
res$idPeriod[res$idPeriod == "Low light acclimating"] <- "Low light"
res$idPeriod[res$idPeriod == "Recovery acclimating"] <- "Recovery 1"
res$idPeriod[res$idPeriod == "Recovery acclimated"] <- "Recovery 2"

# Compute new variables
res$phi1_minus_phi2 <- res$phi1-res$phi2
res$A1_plus_A2 <- res$A1+res$A2
res$daily_amp_mod <- res$E_day_max_mod - res$E_night_min_mod
res$daily_amp_obs <- res$E_day_max_obs - res$E_night_min_obs

# Switch to hours (since the start of the day or night) where relevant
res[, c("phi1", "phi2", "phi1_minus_phi2", "t_day_max_mod", "t_day_max_obs")] <- res[, c("phi1", "phi2", "phi1_minus_phi2", "t_day_max_mod", "t_day_max_obs")]*24
res[, c("t_night_min_mod", "t_night_min_obs")] <- (res[, c("t_night_min_mod", "t_night_min_obs")] - 0.5) *24

# Remove one aberrant phase fit (mex1 #349 low light)
# hist(res$phi1)
# hist(res$phi1[res$phi1 <= (-1)])
# hist(res$phi1[res$phi1 <= (-2)]) # only one observation below -4h for phi1, reaching -11h...
# res[which.min(res$phi1), ]
res[which.min(res$phi1), c("phi1", "phi2", "phi1_minus_phi2")] <- NA


# Filter the phases for low amplitudes (uncomment if required)
res_dupl <- res
#res_dupl[res_dupl$A1 < 0.05 & !is.na(res_dupl$A1), c("phi1", "phi1_minus_phi2", "t_day_max_mod", "t_day_max_obs", "t_night_min_mod", "t_night_min_obs")] <- NA
#res_dupl[res_dupl$A2 < 0.02 & !is.na(res_dupl$A2), c("phi2", "phi1_minus_phi2")] <- NA

# Create the average matrix
avg <- aggregate(res_dupl[, !colnames(res_dupl) %in% c("idExperiment", "idGenotype", "idWatering", "idPot", "idPeriod", "diel_trend", "diel_trend_baseline", "t1_fit", "t2_fit", "adjR2", "RMSE", "E_EON_before_obs", "E_EON_before_mod", "A_rapid_op_transition_obs", "A_rapid_op_stable_obs",
                                                       "E_day_mean_mod", "E_day_mean_obs", "E_night_mean_mod", "E_night_mean_obs", "E_EON_mod", "E_EON_obs", "E_day_max_obs", "E_day_max_mod", "E_night_min_obs", "E_night_min_mod", "phi1_minus_phi2", "A1_plus_A2",
                                                       "acclim_slope")],
                 by = list (idWatering = res_dupl$idWatering, idGenotype = res_dupl$idGenotype, idPeriod = res_dupl$idPeriod),
                 FUN = mean, na.rm = T)
mat <- as.matrix(avg[-(1:3)])
rownames(mat) <- paste(avg[,1], avg[,2], avg[,3])


# Set vector of column names
LABCOL <- set_names_units("names")[match(colnames(mat), set_names_units("par"))] 

# Set the colors for the traits
trait_colors <- c(palette("default")[1], palette("default")[3], palette("default")[3],
                  "darkorange", "darkorange3", "darkorange3", "darkorange",
                  "deepskyblue2", "dodgerblue4", "dodgerblue4", "deepskyblue2",
                  "aquamarine4", palette("default")[3], "coral3", "coral3", "coral1", "coral1",
                  "darkorange", "darkorange3", "darkorange3", "darkorange",
                  "deepskyblue2", "dodgerblue4", "dodgerblue4", "deepskyblue2",
                  "palegreen3", "palegreen3")



#-----------------------------#
#2.1# Well-watered only

# This heatmap is similar to the one presented in the first bioRxiv submission:
# https://www.biorxiv.org/content/10.1101/2022.10.07.511256v1
# except that:
#   - more data are available
#       * retrieval of plants that were initially discarded thanks to new outlier detection process
#       * no filtering of phase for low amplitude (which influences the position of mex1 low light and high CO2, but not the clustering of parameters)
#   - the rows and columns are not manually reordered

mat_WW <- mat[avg$idWatering == "WW", ]
mat_WW <- scale(mat_WW)

LABROW_WW <- LABROW[!grepl("WS", LABROW)]

row_ha <- rowAnnotation(Environment = avg$idPeriod[avg$idWatering == "WW"],
                        col = list(Environment = vector_colperiod),
                        show_annotation_name = F,
                        annotation_legend_param = list(Environment = list(labels = c("Control", expression("High CO"[2]), "Low light", "Recov.1", "Recov.2"),
                                                                          nrow = 1,
                                                                          title_position = "topleft",
                                                                          labels_gp = gpar(fontsize = 8),
                                                                          title_gp = gpar(fontsize = 10, fontface = "plain"))))

hmp_WW <- Heatmap(mat_WW, name = "Trait value (scaled and centered)",
                  col = colorpanel(256, "purple", "black", "orange"),
                  #rect_gp = gpar(col = "gray", lwd = 0.01),
                  row_labels = rep(LABROW_WW, 5),
                  column_labels = LABCOL,
                  row_names_gp = gpar(fontsize = 7, col = ColorsTrt$col[match(avg$idGenotype[avg$idWatering == "WW"], ColorsTrt$idGenotype)]),
                  column_names_gp = gpar(fontsize = 8, col = trait_colors),
                  column_names_rot = 45,
                  left_annotation = row_ha,
                  heatmap_legend_param = list(direction = "horizontal",
                                              title_position = "topcenter",
                                              legend_width = unit(6, "cm"),
                                              title_gp = gpar(fontsize = 10, fontface = "plain"),
                                              labels_gp = gpar(fontsize = 8)),
                  row_split = 4, column_split = 4,
                  #row_gap = unit(c(2, 4, 2), "pt"), column_gap = unit(c(2, 4, 2), "pt"),
                  row_title = NULL, column_title = "Model without acclimation - well-watered only")

#dev.new(width = 6.5, height = 7.5, unit = "in", noRStudioGD = T)
#draw(hmp_WW, heatmap_legend_side = "bottom")


#-----------------------------#
#2.2# Including water stress

mat_all <- scale(mat)

row_ha <- rowAnnotation(Environment = avg$idPeriod,
                        col = list(Environment = vector_colperiod),
                        show_annotation_name = F,
                        annotation_legend_param = list(Environment = list(labels = c("Control", expression("High CO"[2]), "Low light", "Recov.1", "Recov.2"),
                                                                          nrow = 1,
                                                                          title_position = "topleft",
                                                                          labels_gp = gpar(fontsize = 8),
                                                                          title_gp = gpar(fontsize = 10, fontface = "plain"))))

hmp <- Heatmap(mat_all, name = "Trait value (scaled and centered)",
               col = colorpanel(256, "purple", "black", "orange"),
               #rect_gp = gpar(col = "gray", lwd = 0.01),
               row_labels = rep(LABROW, 5),
               column_labels = LABCOL,
               row_names_gp = gpar(fontsize = 7, col = ColorsTrt$col[match(avg$idGenotype, ColorsTrt$idGenotype)]),
               column_names_gp = gpar(fontsize = 8, col = trait_colors),
               column_names_rot = 45,
               left_annotation = row_ha,
               heatmap_legend_param = list(direction = "horizontal",
                                           title_position = "topcenter",
                                           legend_width = unit(6, "cm"),
                                           title_gp = gpar(fontsize = 10, fontface = "plain"),
                                           labels_gp = gpar(fontsize = 8)),
               row_split = 4, column_split = 4,
               #row_gap = unit(c(2, 4, 2), "pt"), column_gap = unit(c(2, 4, 2), "pt"),
               row_title = NULL, column_title = "Model without acclimation - including water stress (WS)")

#dev.new(width = 6.5, height = 7.5, unit = "in", noRStudioGD = T)
#draw(hmp, heatmap_legend_side = "bottom")



#------------------------------------------------------------------------------#
# Plot the heatmaps in a pdf file
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#05c_heatmaps.pdf", sep = "_")), width = 6.5, height = 7.5)
draw(hmp_acclim_WW, heatmap_legend_side = "bottom")
draw(hmp_acclim, heatmap_legend_side = "bottom")
draw(hmp_WW, heatmap_legend_side = "bottom")
draw(hmp, heatmap_legend_side = "bottom")
dev.off()
