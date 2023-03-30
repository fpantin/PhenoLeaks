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
#                              (0)  Prepare environment                        #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rm(list=ls()) # empty environment (if not potential error )

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
meteo <- read.csv(file.path(getwd(), spcs, idExp, "Raw_data", paste0(idExp,"_starch_meteo.csv")),stringsAsFactors=F, sep=',')

# Genotype data
# source: own file
genolist <- read.csv(file.path(getwd(), spcs, idExp, "Raw_data", paste0(idExp,"_starch_genotype_list.csv")),stringsAsFactors=F, sep=';')

# Leaf surface data
# source: output from a Python pipeline to extract rosette surface area
surface <- read.csv(file.path(getwd(), spcs, idExp, "Raw_data", paste0(idExp,"_starch_leafsurface.csv")),stringsAsFactors=F, sep=',')

# Gravimetric data
# source: Phenopsis DB
grv <- read.csv(file.path(getwd(), spcs, idExp, "Raw_data", paste0(idExp,"_starch_gravimetric.csv")),header = T,sep=",")


#------------------------------------------------------------------------------#
#                       Input file preparation                                 #
#------------------------------------------------------------------------------#
# prepare each inputfile: renaming columns, adding new vectors as a decimalDay

# Meteo data
meteo$decimalDay <- decimalDay(column=meteo$date) # add decimal day
colnames(meteo) <- c("idChamber","date","Air.humidity","Air.temperature","Light","VPD","decimalDay") # rename colum names

# Genotype data
genolist <- genolist[,c("idPotManip","idGenotype","Treatment","Analysis","Sowing")] # select column names (deselect idPot)
genolist$Sowing_decimalDay <- decimalDay(column=genolist$Sowing) # add decimal day
inpots <- genolist$idPotManip[genolist$Analysis == "TR"]
colnames(genolist)[match("idPotManip",colnames(genolist))] <- "idPot" # rename column

# Leaf surface data
# need to update this one using raw input file (as for the other experiments)
surface <- surface[,c("idPotManip","decimalDay","Areamm2","outlier")] # select column names (deselect idPot)
colnames(surface)[match("idPotManip",colnames(surface))] <- "idPot" # rename column
surf_coef <- surf_fit(surface) # calculate the statistical model that fits best the evolution of rosette growth data

# Gravimetric data
grv <- grv[grv$idPotManip %in% inpots,]
grv <- grv[,c("idPotManip","date","After.watering.weight")]
colnames(grv)[match("After.watering.weight",colnames(grv))] <- "weight" # rename column
colnames(grv)[match("idPotManip",colnames(grv))] <- "idPot" # rename column
grv$decimalDay <- decimalDay(column=grv$date) # add decimal day
grv$hourday <- as.numeric(sub(".*\\.","0.",sprintf("%.10f", round(grv$decimalDay,10)) )) # add decimal hour on day
grv <- grv[which(grv$date > from & grv$date < to),] #select the data from the start and end of gravimetric experiment
grv <- grv[!is.na(grv$weight),]# delete all NA's
if(anyDuplicated(grv) > 0){ grv <- grv[!duplicated(grv),] }# remove identical lines
grv <- grv[order(grv$idPot,grv$decimalDay),] # order per pot and decimalDay
grv$ID <- rep(NA,nrow(grv)) # unique chronlogical identifier per measurement per pot
grv$Exp_time <- rep(NA,nrow(grv)) # exp time = decimalDay - day number of first day = continuous exp time from 0,.. to 8,.. days 
for (i in unique(grv$idPot)){
  grv[which(grv$idPot == i),]$ID <- seq.int(nrow(grv[which(grv$idPot == i),])) # add pot unique identifiers
  decday <- as.integer(grv$decimalDay[which(grv$idPot == i)][1] ) # first day number of start experiment
  grv$Exp_time[grv$idPot == i] <- grv$decimalDay[grv$idPot == i] - decday # Exp time is from onset of the experiment for each pot
}
# add lightperiod dark night and transition period
darklight = 5 # minutes up or down the dark light transition
if (night == "inversed"){
  grv$lightPeriod <- "light"
  grv$lightPeriod[grv$hourday > startdark & grv$hourday <= startdark + darkperiod] <- "dark"
  
  grv$lightPeriod[grv$hourday > startdark - darklight/60/24 & grv$hourday <= startdark + darklight/60/24] <- "darklight"
  
  grv$lightPeriod[grv$hourday > enddark - darklight/60/24 & grv$hourday <= enddark + darklight/60/24] <- "darklight"
}

if (night == "normal"){
  grv$lightPeriod <- "dark"
  grv$lightPeriod[grv$hourday > enddark & grv$hourday <= startdark ] <- "light"
  
  grv$lightPeriod[grv$hourday > startdark - darklight/60/24 & grv$hourday <= startdark + darklight/60/24] <- "darklight"
  
  grv$lightPeriod[grv$hourday > enddark - darklight/60/24 & grv$hourday <= enddark + darklight/60/24] <- "darklight"
  
}
grv <- surface_add(input=grv) # add leaf surface per timepoint per pot. 
grv <- merge(grv,genolist,by="idPot")

#------------------------------------------------------------------------------#
#                       Check                                                  #
#------------------------------------------------------------------------------#

# check plant growth data for outliers and fit
if (!dir.exists(file.path(getwd(), spcs, idExp, "Figures"))) { dir.create(file.path(getwd(), spcs, idExp, "Figures")) }
pdf(file.path(getwd(), spcs, idExp, "Figures", paste(idExp, "Step#01_plant_growth.pdf", sep = "_")), width = 10, height = 10)

for (geno in ColorsTrt$idGenotype){
  # geno = ColorsTrt$idGenotype[1]
  idpots <- genolist$idPot[genolist$idGenotype == geno & genolist$Analysis == "TR"]
  
  input <- surface[surface$idPot %in% idpots,]
  plot(Areamm2~decimalDay,input,type='n',ylab=expression(paste("Rosette surface", " (mm"^2, ")")),ylim=range(surface$Areamm2))
  title(paste0(geno,", filled = outlier"))
  rectangle(x = input$decimalDay,y = c(surface$Areamm2))
  
  for (pot in idpots){
    # pot = idpots[1]
    points(Areamm2~decimalDay,input[input$idPot == pot,],type="p",col = ColorsTrt$col[match(pot,idpots)])
    points(Areamm2~decimalDay,input[input$idPot == pot & input$outlier,],type="p",col=ColorsTrt$col[match(pot,idpots)],pch=16)
    
    if(pot %in% unique(surf_coef$idPot)){
      
      mod <- surf_coef$model[surf_coef$idPot == pot & surf_coef$growth == "continu"]
      slope = surf_coef$slope[surf_coef$idPot == pot & surf_coef$growth == "continu"]
      intercept = surf_coef$intercept[surf_coef$idPot == pot & surf_coef$growth == "continu"]
      
      if (mod == "lm"){
        curve(intercept+slope*x ,from = min(surface$decimalDay), to = max(surface$decimalDay),add = T,col = ColorsTrt$col[match(pot,idpots)])
      }
      if (mod == "log"){
        curve(exp(intercept+slope*x), from = min(surface$decimalDay), to = max(surface$decimalDay),add = T,col = ColorsTrt$col[match(pot,idpots)])
      }
      
    }
  }
  legend("topleft",legend = idpots, fill=ColorsTrt$col[1:length(idpots)])
  
  
  
}
dev.off()

#------------------------------------------------------------------------------#
#                       Cleaning                                               #
#------------------------------------------------------------------------------#




#------------------------------------------------------------------------------#
#                       Cleaning                                               #
#------------------------------------------------------------------------------#
