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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                              (1)  Load libraries                             #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#options(repos = "https://cran.rstudio.com/") # RStudio
for (pkg in c("gplots"))
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


dir_PhenoLeaks <- file.path(dirname(getwd()), "PhenoLeaks")
source(file.path(dir_PhenoLeaks, "PhenoLeaks_generic.R"))
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
#                         Import and manage fitted data                        #
#------------------------------------------------------------------------------#

# Import data
#res <- read.csv(file.path(getwd(), "Fitted_data", "results_fit_use_VPD.csv"))
res <- read.csv(file.path(getwd(), "Fitted_data", "results_fit_acclim.csv"))
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
  pdf(file.path(getwd(), "Figures", paste(idExp, "_Step#05a_stat_", gen_file, ".pdf", sep = "")), width = 8, height = 5)
  
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
  graph2ppt(file = file.path(getwd(), "Figures", "PPTX", paste(idExp, "_stat_main_", gen_file, ".pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = (4-0.3)*2/3)
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
  graph2ppt(file = file.path(getwd(), "Figures", "PPTX", paste(idExp, "_stat_suppl_", gen_file, ".pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = 4-0.3)
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
graph2ppt(file = file.path(getwd(), "Figures", "PPTX", paste(idExp, "_stat_main_", gen_file, "_obs.pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5/2, height = 4-0.3)
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
graph2ppt(file = file.path(getwd(), "Figures", "PPTX", paste(idExp, "_stat_suppl_", gen_file, ".pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = (4-0.3)*2/3)
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
graph2ppt(file = file.path(getwd(), "Figures", "PPTX", paste(idExp, "_stat_main_", gen_file, "_fit.pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = (4-0.3)*2/3)
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
#graph2ppt(file = file.path(getwd(), "Figures", "PPTX", paste(idExp, "_stat_acclim.pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = 4.5)
#dev.off()


#------------------------------------------------------------------------------#
#         Effect of irrigation regimes across periods (Col-0 and pgm)          #
#------------------------------------------------------------------------------#

# Plot all variables in a PDF file
pdf(file.path(getwd(), "Figures", paste(idExp, "Step#05b_stat_irrigation_Col_pgm.pdf", sep = "_")), width = 8, height = 5)
  
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
graph2ppt(file = file.path(getwd(), "Figures", "PPTX", paste(idExp, "_stat_irrigation_Col_pgm.pptx", sep = "")), paper = "A4", orient = "portrait", width = 6.5, height = 6)
dev.off()

res$idWatering[res$idWatering == "Well-watered"] <- "WW"
res$idWatering[res$idWatering == "Water stress"] <- "WS"









#------------------------------------------------------------------------------#
#                                    Heatmap                                   #
#------------------------------------------------------------------------------#


#require(ComplexHeatmap)
#require(gplots)

# Filter the phases for low amplitudes
res_dupl <- res
res_dupl[res_dupl$A1 < 0.05 & !is.na(res_dupl$A1), c("phi1", "phi1_minus_phi2", "t_day_max_mod", "t_day_max_obs", "t_night_min_mod", "t_night_min_obs")] <- NA
res_dupl[res_dupl$A2 < 0.02 & !is.na(res_dupl$A2), c("phi2", "phi1_minus_phi2")] <- NA

# Create the average matrix
avg <- aggregate(res_dupl[, -which(colnames(res_dupl) %in% c("idExperiment", "idGenotype", "idWatering", "idPot", "idPeriod", "diel_trend", "diel_trend_baseline", "t1_fit", "t2_fit", "adjR2", "RMSE", "E_EON_before_obs", "E_EON_before_mod", "A_rapid_op_transition_obs", "A_rapid_op_stable_obs"))],
                 by = list (idWatering = res_dupl$idWatering, idGenotype = res_dupl$idGenotype, idPeriod = res_dupl$idPeriod),
                 FUN = mean, na.rm = T)
mat <- as.matrix(avg[-(1:3)])
mat <- mat[, !colnames(mat) %in% c("E_day_mean_mod", "E_day_mean_obs", "E_night_mean_mod", "E_night_mean_obs", "E_EON_mod", "E_EON_obs", "E_day_max_obs", "E_day_max_mod", "E_night_min_obs", "E_night_min_mod", "phi1_minus_phi2", "A1_plus_A2")]
row.names(mat) <- paste(avg[,1], avg[,2], avg[,3])

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
LABCOL <- c(expression("E"["diel"]),
            expression("A"["rapid op"]),
            expression("A"["rapid clo"]),
            expression("t"["day max"]),
            expression(sigma["day"]),
            expression(Delta["day"]),
            expression(Sigma["preclo"]),
            expression("t"["night min"]),
            expression(sigma["night"]),
            expression(Delta["night"]),
            expression(Sigma["preop"]),
            expression("E"["mean"]),
            expression("A"["SQW"]),
            expression("A"["1"]),
            expression(varphi["1"]),
            expression("A"["2"]),
            expression(varphi["2"]),
            expression(hat("t")["day max"]),
            expression(hat(sigma)["day"]),
            expression(hat(Delta)["day"]),
            expression(hat(Sigma)["preclo"]),
            expression(hat("t")["night min"]),
            expression(hat(sigma)["night"]),
            expression(hat(Delta)["night"]),
            expression(hat(Sigma)["preop"]),
            expression(hat("A")["diel"]),
            expression("A"["diel"]))

# Set vector of column names for PPTX export
# (several symbols are not supported by 'graph2ppt()', need to add them manually on the PPTX file)
LABCOL_PPTX <- c(expression("E"["diel"]),
                 expression("A"["rapid op"]),
                 expression("A"["rapid clo"]),
                 expression("t"["day max"]),
                 expression("s"["day"]),#expression(sigma["day"]),
                 expression("D"["day"]),#expression(Delta["day"]),
                 expression("S"["preclo"]),#expression(Sigma["preclo"]),
                 expression("t"["night min"]),
                 expression("s"["night"]),#expression(sigma["night"]),
                 expression("D"["night"]),#expression(Delta["night"]),
                 expression("S"["preop"]),#expression(Sigma["preop"]),
                 expression("E"["mean"]),
                 expression("A"["SQW"]),
                 expression("A"["1"]),
                 expression("j"["1"]),#expression(varphi["1"]),
                 expression("A"["2"]),
                 expression("j"["2"]),#expression(varphi["2"]),
                 expression("t"["day max"]),#expression(hat("t")["day max"]),
                 expression("s"["day"]),#expression(hat(sigma)["day"]),
                 expression("D"["day"]),#expression(hat(Delta)["day"]),
                 expression("S"["preclo"]),#expression(hat(Sigma)["preclo"]),
                 expression("t"["night min"]),#expression(hat("t")["night min"]),
                 expression("s"["night"]),#expression(hat(sigma)["night"]),
                 expression("D"["night"]),#expression(hat(Delta)["night"]),
                 expression("S"["preop"]),#expression(hat(Sigma)["preop"]),
                 expression("A"["diel"]),#expression(hat("A")["diel"]),
                 expression("A"["diel"]))


# Set the colors for the traits
trait_colors <- c(palette("default")[1], palette("default")[3], palette("default")[3],
                  "darkorange", "darkorange3", "darkorange3", "darkorange",
                  "deepskyblue2", "dodgerblue4", "dodgerblue4", "deepskyblue2",
                  "aquamarine4", palette("default")[3], "coral3", "coral3", "coral1", "coral1",
                  "darkorange", "darkorange3", "darkorange3", "darkorange",
                  "deepskyblue2", "dodgerblue4", "dodgerblue4", "deepskyblue2",
                  "palegreen3", "palegreen3")

# Set the colors for the periods
vector_colperiod <- ColorsPeriod$col
names(vector_colperiod) <- ColorsPeriod$idPeriod


#1# Heatmap without water stress
mat_WW <- mat[avg$idWatering == "WW", ]
mat_WW <- scale(mat_WW)

LABROW_WW <- LABROW[!grepl("WS", LABROW)]

row_ha <- rowAnnotation(Environment = avg[avg$idWatering == "WW", 3],
                        col = list(Environment = vector_colperiod),
                        show_annotation_name = F,
                        annotation_legend_param = list(Environment = list(labels = c("Control", expression("High CO"[2]), "Low light", "Recov.1", "Recov.2"),
                                                                          nrow = 1,
                                                                          title_position = "topleft",
                                                                          labels_gp = gpar(fontsize = 8),
                                                                          title_gp = gpar(fontsize = 10, fontface = "plain"))))

hmp <- Heatmap(mat_WW, name = "Trait value (scaled and centered)",
               col = colorpanel(256, "purple", "black", "orange"),
               #rect_gp = gpar(col = "gray", lwd = 0.01),
               row_labels = rep(LABROW_WW, 5), 
               column_labels = LABCOL,
               row_names_gp = gpar(fontsize = 7, col = ColorsTrt$col[match(avg[avg$idWatering == "WW", 2], ColorsTrt$idGenotype)]),
               column_names_gp = gpar(fontsize = 8, col = trait_colors),
               column_names_rot = 45,
               left_annotation = row_ha,
               heatmap_legend_param = list(direction = "horizontal",
                                           title_position = "topcenter",
                                           legend_width = unit(6, "cm"),
                                           title_gp = gpar(fontsize = 10, fontface = "plain"),
                                           labels_gp = gpar(fontsize = 8)),
               # cluster_rows = reorder(as.dendrogram(hclust(dist(mat))),
               #                        wts = ifelse(rownames(mat) == "pgm-1 0-1-2", 0, 1) + ifelse(rownames(mat) == "dpe2-5 0-1-2", 1, 0),
               #                        agglo.FUN = mean),
               # cluster_columns = reorder(as.dendrogram(hclust(dist(t(mat)))),
               #                           wts = ifelse(colnames(mat) == "A_SQW", 0, 2) + ifelse(colnames(mat) == "E_diel_obs", 0, 1) + ifelse(colnames(mat) == "sigma_day_obs", 1, 0) + ifelse(colnames(mat) == "A2", 2, 0) + ifelse(colnames(mat) == "Sigma_preop_mod", 3, 0),
               #                           agglo.FUN = mean),
               row_split = 5, column_split = 4,
               #row_gap = unit(c(2, 4, 2), "pt"), column_gap = unit(c(2, 4, 2), "pt"),
               row_title = NULL, column_title = NULL)

pdf(file.path(getwd(), "Figures", paste(idExp, "Step#05d_heatmap_WW.pdf", sep = "_")), width = 6.5, height = 7.5)
draw(hmp, heatmap_legend_side = "top")
dev.off()

dev.new(width = 6.5, height = 7.5, unit = "in", noRStudioGD = T)
draw(hmp, heatmap_legend_side = "top")

graph2ppt(file = "Output/Heatmap.pptx", paper = "A4", orient = "portrait", width = 6.5, height = 7.5)






#2# With water stress

mat_all <- scale(mat)

row_ha <- rowAnnotation(Environment = avg[,3],
                        col = list(Environment = vector_colperiod),
                        show_annotation_name = F,
                        annotation_legend_param = list(Environment = list(labels = c("Control", expression("High CO"[2]), "Low light", "Recov.1", "Recov.2"),
                                                                          nrow = 1,
                                                                          title_position = "topleft",
                                                                          labels_gp = gpar(fontsize = 8),
                                                                          title_gp = gpar(fontsize = 10, fontface = "plain"))))
















# export_anova2_posthoc <- function (file,
#                                    dat, # the dataframe
#                                    Y, # the name of the dependent variable (character)
#                                    MAIN_FACTOR = "idPeriod",
#                                    GROUP = "idGenotype")
#   {
#   #require(rstatix)
#   
#   write.table(paste(rep("~", nchar(Y)), sep = "", collapse = ""), file, row.names = F, col.names = F, sep = ",", append = T)
#   write.table(Y, file, row.names = F, col.names = F, sep = ",", append = T)
#   write.table(paste(rep("~", nchar(Y)), sep = "", collapse = ""), file, row.names = F, col.names = F, sep = ",", append = T)
#   write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
#   
#   # Two-way (mixed) ANOVA
#   if ("idPeriod" %in% c(MAIN_FACTOR, GROUP))
#     {
#     between_var <- c(MAIN_FACTOR, GROUP)[!c(MAIN_FACTOR, GROUP) %in% "idPeriod"]
#     res.aov <- anova_test(data = dat,
#                           dv = all_of(Y),
#                           between = all_of(between_var),
#                           within = idPeriod,
#                           wid = idPot)
#     }
#   else
#     {
#     res.aov <- anova_test(data = dat,
#                           dv = all_of(Y),
#                           between = c(all_of(MAIN_FACTOR), all_of(GROUP)))
#     }
#   table.aov <- get_anova_table(res.aov)
#   write.table("ANOVA table", file, row.names = F, col.names = F, sep = ",", append = T)
#   write.table(table.aov, file, row.names = F, sep = ",", append = T)
#   write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
#   
#   # Pairwise t-tests
#   write.table(paste(GROUP, "pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
#   if (table.aov[3, "p"] < 0.05) write.table("WARNING - The significant interaction may make these comparisons irrelevant", file, row.names = F, col.names = F, sep = ",", append = T)
#   pwc <- dat %>%
#     pairwise_t_test(., reformulate(GROUP, Y), p.adjust.method = "bonferroni")
#   write.table(pwc, file, row.names = F, sep = ",", append = T)
#   write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
#   
#   write.table(paste(MAIN_FACTOR, "pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
#   if (table.aov[3, "p"] < 0.05) write.table("WARNING - The significant interaction may make these comparisons irrelevant", file, row.names = F, col.names = F, sep = ",", append = T)
#   pwc <- dat %>%
#       pairwise_t_test(., reformulate(MAIN_FACTOR, Y), p.adjust.method = "bonferroni")
#   write.table(pwc, file, row.names = F, sep = ",", append = T)
#   write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
#   
#   write.table(paste(GROUP, "conditional pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
#   pwc <- dat %>%
#     group_by(get(MAIN_FACTOR)) %>%
#     pairwise_t_test(., reformulate(GROUP, Y), p.adjust.method = "bonferroni")
#   colnames(pwc)[1] <- MAIN_FACTOR
#   write.table(pwc, file, row.names = F, sep = ",", append = T)
#   write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
#   
#   write.table(paste(MAIN_FACTOR, "conditional pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
#   pwc <- dat %>%
#     group_by(get(GROUP)) %>%
#     pairwise_t_test(., reformulate(MAIN_FACTOR, Y), p.adjust.method = "bonferroni")
#   colnames(pwc)[1] <- GROUP
#   write.table(pwc, file, row.names = F, sep = ",", append = T)
#   write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
#   
#   write.table(paste(rep("_", 150), sep = "", collapse = ""), file, row.names = F, col.names = F, sep = ",", append = T)  
#   write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
#   }
# 
# 
# 
#                             
# 
#          
# 
# 
# 
# 
# require(rstatix)
# gen <- c("Col-0", "pgm-1", "sex1-3")
# var <- "phi2"#"phi1"#"E_diel_obs" #
# res.aov <- anova_test(data = res[res$idGenotype %in% gen & res$idWatering=="WW", ],
#                       dv = all_of(var),
#                       wid = idPot,
#                       between = idGenotype,
#                       within = idPeriod)
# a=get_anova_table(res.aov)
# 
# b=res[res$idGenotype %in% gen & res$idWatering=="WW", ] %>%
#   group_by(idPeriod) %>%
#   pairwise_t_test(., reformulate("idGenotype", var), p.adjust.method = "bonferroni")
# c=res[res$idGenotype %in% gen & res$idWatering=="WW", ] %>%
#   pairwise_t_test(., reformulate("idGenotype", var), p.adjust.method = "bonferroni")
# 
# d=res[res$idGenotype %in% gen & res$idWatering=="WW", ] %>%
#   group_by(idGenotype) %>%
#   pairwise_t_test(., reformulate("idPeriod", var), p.adjust.method = "bonferroni")
# e=res[res$idGenotype %in% gen & res$idWatering=="WW", ] %>%
#   pairwise_t_test(., reformulate("idPeriod", var), p.adjust.method = "bonferroni")
# 
# write.table(a, "test.csv", row.names = F, sep = ",", col.names = T)
# write.table(b, "test.csv", row.names = F, sep = ",", col.names = T, append = T)
# 
# 
# 
# 
# 
# 
# 
# 
# dat_paired <- res[!res$outlier & res$geno %in% gen, c("geno", "pot", "idperiod", var)]
# names(dat_paired)[4] <- "DV" # dependent variable
# dat_paired_full <- data.frame(geno = NULL, pot = NULL, idperiod = NULL)
# for (g in gen)
# {
#   dat_paired_full <- rbind(dat_paired_full,
#                            expand.grid(geno = g, pot = unique(dat_paired$pot[dat_paired$geno == g]), idperiod = unique(dat_paired$idperiod)))
# }
# trt_paired <- paste(dat_paired$geno, dat_paired$pot, dat_paired$idperiod)
# trt_paired_full <- paste(dat_paired_full$geno, dat_paired_full$pot, dat_paired_full$idperiod)
# dat_paired_full$DV <- NA
# dat_paired_full$DV[match(trt_paired, trt_paired_full)] <- dat_paired$DV
# 
# ## WARNING - NOT SURE HOW MISSING DATA ARE HANDLED!!!
# # If GxE interaction
# dat_paired_full %>%
#   group_by(geno) %>%
#   pairwise_t_test(DV ~ idperiod, paired = T, p.adjust.method = "bonferroni")
# # If no GxE interaction
# dat_paired_full %>%
#   pairwise_t_test(DV ~ idperiod, paired = T, p.adjust.method = "bonferroni")
# 
# aggregate(dat_paired_full$DV, by = list(geno = dat_paired_full$geno, idperiod=dat_paired_full$idperiod), FUN = mean, na.rm = T)
# aggregate(dat_paired_full$DV, by = list(geno = dat_paired_full$geno, idperiod=dat_paired_full$idperiod), FUN = function (x) { sd(x, na.rm = T) / sqrt(length(na.omit(x)) - 1) })
# 
# 
# 
# 
# 
