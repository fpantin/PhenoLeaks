################################################################################
#                                                                              #
#         PhenoLeaks - Step #01 - Estimation of transpiration                  #
#                                                                              #
#              Script to clean raw data and estimate transpiration values      #
#                                                                              #
#                             Adriaan Westgeest, 2023                          #
#                                                                              #
################################################################################



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                              (1)  Load libraries                             #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#options(repos = "https://cran.rstudio.com/") # RStudio
for (pkg in c("splitstackshape"))
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

# Note that the working directory is expected to be the one of the PhenoLeaks project directory,
# e.g. in RStudio: Session > Set Working Directory > To Project Directory
dir_PhenoLeaks <- file.path(getwd(), "_core")
source(file.path(dir_PhenoLeaks, "PhenoLeaks_generic.R"))
source(file.path(dir_PhenoLeaks, "PhenoLeaks_preparation.R"))
source(file.path(dir_PhenoLeaks, "PhenoLeaks_transpiration.R"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                 (3)  Retrieve the features of the experiment                 #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

idExp = "C2M47"
source(file.path(getwd(), "Arabidopsis", idExp, "Step#00_define_experiment.R"))

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
#                       Import raw data                                        #
#------------------------------------------------------------------------------#

# Meteo data
meteo <- read.csv(file.path(getwd(), "Arabidopsis", idExp, "Raw_data", paste0(idExp,"_starch_meteo.csv")),stringsAsFactors=F, sep=',')

# Genotype data
genolist <- read.csv(file.path(getwd(), "Arabidopsis", idExp, "Raw_data", paste0(idExp,"_starch_genotype_list.csv")),stringsAsFactors=F, sep=';')

# Leaf surface data
surface <- read.csv(file.path(getwd(), "Arabidopsis", idExp, "Raw_data", paste0(idExp,"_starch_leafsurface.csv")),stringsAsFactors=F, sep=',')

# Gravimetric data
grv <- read.csv(file.path(getwd(), "Arabidopsis", idExp, "Raw_data", paste0(idExp,"_starch_gravimetric.csv")),header = T,sep=",")


#------------------------------------------------------------------------------#
#                       Processing                                             #
#------------------------------------------------------------------------------#

# Meteo data
meteo$decimalDay <- decimalDay(column=meteo$date) # add decimal day
colnames(meteo) <- c("idChamber","date","Air.humidity","Air.temperature","Light","VPD","decimalDay") # rename colum names

# Genotype data
genolist$Sowing_decimalDay <- decimalDay(column=genolist$Sowing) # add decimal day

# Leaf surface data

# Gravimetric data
grv = grv[,c("idPotManip","date","After.watering.weight")]
colnames(grv)[3] = "weight"