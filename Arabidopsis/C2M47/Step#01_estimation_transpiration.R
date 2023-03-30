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
#                 (2)  Retrieve the features of the experiment                 #
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
#                             (3)  Source functions                            #
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
#                             (4)  Run the script                              #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#------------------------------------------------------------------------------#
#                       Import raw data                                        #
#------------------------------------------------------------------------------#

# Meteo data 
# source: Phenopsis DB
meteo <- read.csv(file.path(getwd(), "Arabidopsis", idExp, "Raw_data", paste0(idExp,"_starch_meteo.csv")),stringsAsFactors=F, sep=',')

# Genotype data
# source: own file
genolist <- read.csv(file.path(getwd(), "Arabidopsis", idExp, "Raw_data", paste0(idExp,"_starch_genotype_list.csv")),stringsAsFactors=F, sep=';')

# Leaf surface data
# source: output from a Python pipeline to extract rosette surface area
surface <- read.csv(file.path(getwd(), "Arabidopsis", idExp, "Raw_data", paste0(idExp,"_starch_leafsurface.csv")),stringsAsFactors=F, sep=',')

# Gravimetric data
# source: Phenopsis DB
grv <- read.csv(file.path(getwd(), "Arabidopsis", idExp, "Raw_data", paste0(idExp,"_starch_gravimetric.csv")),header = T,sep=",")


#------------------------------------------------------------------------------#
#                       Input file preparation                                 #
#------------------------------------------------------------------------------#
# prepare each inputfile: renaming columns, adding new vectors as a decimalDay

# Meteo data
meteo$decimalDay <- decimalDay(column=meteo$date) # add decimal day
colnames(meteo) <- c("idChamber","date","Air.humidity","Air.temperature","Light","VPD","decimalDay") # rename colum names

# Genotype data
genolist$Sowing_decimalDay <- decimalDay(column=genolist$Sowing) # add decimal day

# Leaf surface data
# need to update this one using raw input file (as for the other experiments)

# Gravimetric data
grv <- grv[,c("idPotManip","date","After.watering.weight")]
colnames(grv)[match("After.watering.weight",colnames(grv))] <- "weight" # rename column
grv$decimalDay <- decimalDay(column=grv$date) # add decimal day
# grv$dayofyear = as.integer(grv$decimalDay) # add day of the year
grv$hourday <- as.numeric(sub(".*\\.","0.",sprintf("%.10f", round(grv$decimalDay,10)) )) # add decimal hour on day
grv <- grv[which(grv$date > from & grv$date < to),] #select the data from the start and end of gravimetric experiment
grv <- grv[!is.na(grv$weight),]# delete all NA's
if(anyDuplicated(grv) > 0){ grv <- grv[!duplicated(grv),] }# remove identical lines
grv <- grv[order(grv$idPotManip,grv$decimalDay),] # order per pot and decimalDay

#------------------------------------------------------------------------------#
#                       Cleaning                                               #
#------------------------------------------------------------------------------#
# data processing, only gravimetric data.

# Gravimetric data

#add experimental time, essential for rehydrations script

grv$ID <- rep(NA,nrow(grv))
grv$Exp_time <- rep(NA,nrow(grv))

for (i in unique(grv$idPotManip)){
  
  grv[which(grv$idPotManip == i),]$ID <- seq.int(nrow(grv[which(grv$idPotManip == i),]))
  
  decday <- as.integer(grv$decimalDay[which(grv$idPotManip == i)][1] )
  
  grv$Exp_time[grv$idPotManip == i] <- grv$decimalDay[grv$idPotManip == i] - decday
}
head(grv)
range(grv$date)