################################################################################
#                                                                              #
#       PhenoLeaks - Step #01 - Computing transpiration from weight data       #
#                                                                              #
#                       Script to clean raw weight data                        #
#                      and estimate transpiration values                       #
#                                                                              #
#                             Adriaan Westgeest, 2023                          #
#                                                                              #
################################################################################

#>                                                                            <#
#     >                                                                  <     #
#           >                                                      <           #
#                 >                                          <                 #
#                       >                              <                       #
#                             >                  <                             #
#                                     C3M31                                    #
#                             >                  <                             #
#                       >                              <                       #
#                 >                                          <                 #
#           >                                                      <           #
#     >                                                                  <     #
#>                                                                            <#

                                                                     
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
# package splitstackshape contains a function Csplit to split columns based on a seperator.
for (pkg in c("splitstackshape", "here", "scales"))
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

# Here the file "Step#00_define_experiment.R" should be sourced.

  ## OPTION 1: open the file and source it manually

  ## OPTION 2: set the correct file path and source it from this script, e.g.:
  source(file.path(here::here(), "Arabidopsis", "C3M31", "Step#00_define_experiment.R"))

# Now the species and ID of the experiment can be found by entering:
#c(spcs, idExp)

# so that the directory of the experiment is now explicitly defined as:
dir_Exp <- file.path(here::here(), spcs, idExp)

# To check all constants, enter:
#set_constants_C3M31()

# To check the colors, enter:
#set_colors_C3M31()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                             (3)  Source functions                            #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Note that paths are built relative to the root directory of the PhenoLeaks project.
dir_PhenoLeaks <- here::here("_core")
source(file.path(dir_PhenoLeaks, "PhenoLeaks_generic.R"))
source(file.path(dir_PhenoLeaks, "PhenoLeaks_preparation.R"))
source(file.path(dir_PhenoLeaks, "PhenoLeaks_weight_to_transpiration.R"))
source(file.path(dir_PhenoLeaks, "PhenoLeaks_graphics.R"))



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
meteo <- read.csv(file.path(dir_Exp, "Raw_data", paste0(idExp,"_meteo.csv")),stringsAsFactors=F, sep=',')
if(dim(meteo)[2] == 1){meteo <- read.csv(file.path(dir_Exp, "Raw_data", paste0(idExp,"_meteo.csv")),stringsAsFactors=F, sep=';')}
# head(meteo)

# Genotype data
# source: own file
genolist <- read.csv(file.path(dir_Exp, "Raw_data", paste0(idExp,"_genotype_list.csv")),stringsAsFactors=F, sep=';')
# head(genolist)

# Leaf surface data
# source: output from a Python pipeline to extract rosette surface area
surface <- read.csv(file.path(dir_Exp, "Raw_data", paste0(idExp,"_leafsurface.csv")),stringsAsFactors=F, sep=',')
if(dim(surface)[2] == 1){surface <- read.csv(file.path(dir_Exp, "Raw_data", paste0(idExp,"_leafsurface.csv")),stringsAsFactors=F, sep=';')}
# head(surface)

# Gravimetric data
# source: Phenopsis DB
grv <- read.csv(file.path(dir_Exp, "Raw_data", paste0(idExp,"_gravimetric.csv")),header = T,sep=",")
if(dim(surface)[2] == 1){grv <- read.csv(file.path(dir_Exp, "Raw_data", paste0(idExp,"_gravimetric.csv")),stringsAsFactors=F, sep=';')}
# head(grv)

# Soil water content data
swc <- read.csv(file.path(dir_Exp, "Raw_data", paste0(idExp,"_soilwatercontent.csv")),header = T,sep=";")

#------------------------------------------------------------------------------#
#                       Input file preparation                                 #
#------------------------------------------------------------------------------#
# prepare each inputfile: renaming columns, adding new vectors as a decimalDay

#------------------------------------------------------------------------------#
# Meteo data
meteo$date <- as.POSIXct(strptime(meteo$date,format= "%Y-%m-%d %H:%M:%S", tz = "UTC")) # now the date/hour is not changing
meteo$decimalDay <- decimalDay(column=meteo$date) # add decimal day
meteo$VPD[meteo$VPD<0] <- NA # negative values
vpd_periods <- c(53+ 11*(1/24),53+ 13*(1/24)) # start-end of vpd period
#head(meteo)

#------------------------------------------------------------------------------#
# Genotype data
genolist$Sowing <- as.POSIXct(genolist$Sowing,format= "%d/%m/%Y", tz = "UTC") # transform to posixct before calculation decimalday
genolist$Sowing_decimalDay <- decimalDay(column=genolist$Sowing) # add decimal day

#------------------------------------------------------------------------------#
# Leaf surface data
# need to update this one using raw input file (as for the other experiments)
surface <- as.data.frame(splitstackshape::cSplit(surface, "Label", sep="-", type.convert=F)) # the file name column
colnames(surface)[match(c("Label_1", "Label_2"),colnames(surface))] = c("Experiment","idPot")
surface$Date1 <- as.POSIXct(strptime(surface$Date, format= "%d/%m/%Y %H:%M:%S",tz = "UTC"))
surface$decimalDay <- decimalDay(column=surface$Date1) # add decimal day
surface$Areamm2 <- surface$Area..cm2. * 100 # recalculation in mm2 units
surface$outlier <- FALSE # no outliers yet identified
surface <- surface[,c("idPot","decimalDay","Areamm2","outlier")] # select column names 
surf_coef <- surf_fit(surface) # calculate the statistical model that fits best the evolution of rosette growth data
surf_coef <- surf_coef[order(as.numeric(surf_coef$idPot)),]

#------------------------------------------------------------------------------#
# Soil water content data
if ("idPotManip" %in% colnames(swc)){
  # when french version was downloaded
  swc <- as.data.frame(splitstackshape::cSplit(swc, "idPotManip", sep="-",type.convert = F)) # split colum with idPotManip to get idPot, need to add as.data.frame to bring back from table to dataframe.
  colnames(swc)[match("idPotManip_2",colnames(swc))] <- "idPot"
  colnames(swc)[match("poidsPotNonTroue",colnames(swc))] <- "nonPerforatedPotWeight"
  colnames(swc)[match("poidsPotTroue",colnames(swc))] <- "perforatedPotWeight"
  colnames(swc)[match("poidsSolSec",colnames(swc))] <- "drySoilWeight"
}else{
  # when english version was downloaded
  swc <- as.data.frame(splitstackshape::cSplit(swc, "idPot", sep="-",type.convert = F)) # split colum with idPotManip to get idPot, need to add as.data.frame to bring back from table to dataframe.
  colnames(swc)[match("idPot_2",colnames(swc))] <- "idPot"
}
swc <- swc[,c("idPot","nonPerforatedPotWeight","perforatedPotWeight","drySoilWeight")] # select column names 

#------------------------------------------------------------------------------#
# Gravimetric data
# preparation  gravimetric data with function
grv <- prep_gravi()
grv$lightPeriod[grv$decimalDay >= vpd_periods[1] & grv$decimalDay <= vpd_periods[2]] <- "VPD"

# just add sowing data to genolist file to add at the end to the transpiration file
genolist$Measuring_decimalDay <- min(grv$dayofyear)
genolist$Days_Sowing_to_measuring <- genolist$Measuring_decimalDay - genolist$Sowing_decimalDay

#------------------------------------------------------------------------------#
#                       Check                                                  #
#------------------------------------------------------------------------------#

# check plant growth data for outliers and fit

if (!dir.exists(file.path(dir_Exp, "Figures"))) { dir.create(file.path(dir_Exp, "Figures")) }

## AJW old version
# pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#01_plant_growth.pdf", sep = "_")), width = 10, height = 10)
# for (geno in unique(ColorsTrt$idGenotype)){
#   # geno = ColorsTrt$idGenotype[1]
#   idpots <- genolist$idPot[genolist$idGenotype == geno]
#   
#   input <- surface[surface$idPot %in% idpots,]
#   plot(Areamm2~decimalDay,input,type='n',ylab=expression(paste("Rosette surface", " (mm"^2, ")")),ylim=range(surface$Areamm2))
#   title(paste0(geno,", filled = outlier"))
#   rectangle(x = input$decimalDay,y = c(surface$Areamm2))
#   
#   for (pot in idpots){
#     # pot = idpots[2]
#     points(Areamm2~decimalDay,input[input$idPot == pot,],type="p",col = Colors[match(pot,idpots)])
#     points(Areamm2~decimalDay,input[input$idPot == pot & input$outlier,],type="p",col=Colors[match(pot,idpots)],pch=16)
#     
#     if(pot %in% unique(surf_coef$idPot)){
#       
#       mod <- surf_coef$model[surf_coef$idPot == pot & surf_coef$growth == "continu"]
#       slope = surf_coef$slope[surf_coef$idPot == pot & surf_coef$growth == "continu"]
#       intercept = surf_coef$intercept[surf_coef$idPot == pot & surf_coef$growth == "continu"]
#       
#       if (mod == "lm"){
#         curve(intercept+slope*x ,from = min(surface$decimalDay), to = max(surface$decimalDay),add = T,col = Colors[match(pot,idpots)])
#       }
#       if (mod == "log"){
#         curve(exp(intercept+slope*x), from = min(surface$decimalDay), to = max(surface$decimalDay),add = T,col = Colors[match(pot,idpots)])
#       }
#       
#     }
#   }
#   legend("topleft",legend = idpots, fill=Colors[1:length(idpots)])
# }
# dev.off()

## New version
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#01a_plant_growth_per_condition.pdf", sep = "_")), width = 10, height = 5)
n.max <- max(aggregate(genolist$idPot[!duplicated(genolist$idPot)], by = list(genolist$idGenotype[!duplicated(genolist$idPot)]), FUN = length)$x) # maximum number of replicates per treatment
for (geno in sort(unique(ColorsTrt$idGenotype)))
  {
  idpots <- genolist$idPot[genolist$idGenotype == geno]
  input <- surface[surface$idPot %in% idpots, ]
  prepare_kin(input, Time_var = "decimalDay", E_var = "Areamm2",
              ylab = expression(paste("Rosette area", " (mm"^2, ")")), 
              main = paste(geno, " - filled = outlier", sep = ""),
              inside = F, irrig_show_mode = "none")
  color = 0
  for (pot in sort(idpots))
    {
    #require(scales)
    color <- color + 1
    points(Areamm2 ~ decimalDay, input[input$idPot == pot, ], type = "p", col = hue_pal()(n.max)[color], cex = 0.5)
    points(Areamm2 ~ decimalDay, input[input$idPot == pot & input$outlier,], type = "p", col = hue_pal()(n.max)[color], cex = 0.5, pch=16)
    
    if (pot %in% unique(surf_coef$idPot))
      {
      mod <- surf_coef$model[surf_coef$idPot == pot & surf_coef$growth == "continu"]
      slope = surf_coef$slope[surf_coef$idPot == pot & surf_coef$growth == "continu"]
      intercept = surf_coef$intercept[surf_coef$idPot == pot & surf_coef$growth == "continu"]
      
      if (mod == "lm")
        {
        curve(intercept+slope*x, from = min(surface$decimalDay), to = max(surface$decimalDay), add = T, col = hue_pal()(n.max)[color])
        }
      if (mod == "log")
        {
        curve(exp(intercept+slope*x), from = min(surface$decimalDay), to = max(surface$decimalDay), add = T, col = hue_pal()(n.max)[color])
        }
      }
    }
  legend("top", ncol = 4, bty = "n",
         legend = paste("Pot ", sort(idpots), sep = ""),
         col = hue_pal()(n.max)[1:color], lty = 1, pch = 21, pt.cex = 0.5)
  rm(input)
  }

dev.off()



#------------------------------------------------------------------------------#
#                               Cleaning                                       #
#------------------------------------------------------------------------------#

#---------------------------STEP 1---------------------------------------------#
output <- Rehy_Corr_v4(input = grv,gap = 1,jumps="positive")
df_rehy <- output$output # corrected output
corr1 <- output$corr1 # 1 run: detected rehydrations

x = c(grv$decimalDay,df_rehy$decimalDay)
y = c(grv$weight,df_rehy$Weight_corr)
xlim = range(x)
ylim = range(y)

## AJW old version
# pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#01_correct_irrigation.pdf", sep = "_")), width = 10, height = 10)
# for (i in genolist$idPot){
#   # i = idpots[3]
#   input <- df_rehy[df_rehy$idPot == i,]
#   if(nrow(input)>1){
#     xlim = range(input$decimalDay)
#     ylim = range(input$weight,input$Weight_corr)
#     plot(weight~decimalDay,input[input$idPot == i,],type='n',ylim = ylim, xlim=xlim, ylab="weight (g)",xlab="decimalDay",main=i)
#     rectangle(x,y)
#     points(weight~decimalDay,input[input$idPot == i,],type="p")
#     if(i %in% unique(corr1$idPot)){
#       # for each unique was an irrigation
#       strt =  corr1$from_id[corr1$idPot == i][1] # take the first one!
#       abline(v=c(corr1$decimalDay[corr1$idPot == i]),col="black")
#       points(Weight_corr~decimalDay,df_rehy[df_rehy$idPot == i & df_rehy$ID > strt,],col="red")
#     }
#     legend("bottomleft",legend = c("corrected data"),fill="red")
#   }
# }
# dev.off()

## New version
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#01b_correct_irrigation.pdf", sep = "_")), width = 10, height = 5)
for (geno in sort(unique(ColorsTrt$idGenotype)))
  {
  idpots <- genolist$idPot[genolist$idGenotype == geno]
  for (pot in sort(idpots))
    {
    input <- df_rehy[df_rehy$idPot == pot, ]
    if (nrow(input) > 1)
      {
      xlim <- range(input$decimalDay)
      ylim <- range(input$weight, input$Weight_corr)
      prepare_kin(input, Time_var = "decimalDay", E_var = "weight",
                  ylab = "Pot weight (g)", ylim = ylim, xlim = xlim,
                  main = paste(geno, " - Pot ", pot, sep = ""),
                  inside = F, irrig_show_mode = "none")
      points(weight ~ decimalDay, data = input, type = "p", cex = 0.5)
      if (pot %in% unique(corr1$idPot))
        {
        # for each unique was an irrigation
        strt <- corr1$from_id[corr1$idPot == pot][1] # take the first one!
        points(Weight_corr ~ decimalDay, df_rehy[df_rehy$idPot == pot & df_rehy$ID > strt,], cex = 0.5, col = "red")
        }
      legend("topright", legend = "corrected data", pch = 1, cex = 0.5, col = "red",
             bty = "n")
      }
    rm(input)
    }
  }
dev.off()

df = df_rehy

# save irrigation information:
if (!dir.exists(file.path(dir_Exp, "Processed_data"))) { dir.create(file.path(dir_Exp, "Processed_data")) }
p2f <- file.path(dir_Exp, "Processed_data", paste(idExp, "info_irrigation.csv", sep = "_"))
write.csv(corr1,p2f,row.names = F)

#---------------------------STEP 2---------------------------------------------#
# Gravimetric outlier detection
# 3 types of outliers:
df$side_outlier <- FALSE # unstable at the start or end
df$obvious_outlier <- FALSE # passing the biological limits
df$hotspot_outlier <- FALSE # determined based on abnormal high influence on the linear regression
df$outlier <- FALSE # if one of the before is true!

# biological limits (depends on day, plant, conditions, PAR or CO2 for example)
# how to define biological limits? average slopes?
min_dark = -1
max_dark = 3.5
min_light = -0.4
max_light = 4.5

for (i in unique(df$idPot)){
  # i = unique(df$idPot)[1]
  # i = "C2M47-217"
  # i =genolist$idPot[1]
  options(warn = 1)
  hotspots <- Outliers_v4_2(time = df$decimalDay[df$idPot == as.character(i)], # time input
                            weight = df$Weight_corr[df$idPot == as.character(i)], # weight input
                            startdark = startdark,
                            darkperiod = darkperiod,
                            min_dark = min_dark,
                            max_dark = max_dark,
                            min_light = min_light,
                            max_light = max_light,
                            surf_i = df$surface[df$idPot == as.character(i)])
  
  outliers <- c()
  out_sides <- unlist(hotspots[2])
  out_obvious <- unlist(hotspots[3])
  hotspots <- unlist(hotspots[1])
  
  # continue with hotspot analysis
  df$side_outlier[df$idPot == i & df$ID %in% c(out_sides) ] <- T # already take this information for the hotspot analysis
  df$obvious_outlier[df$idPot == i & df$ID %in% c(out_obvious) ] <- T # already take this information for the hotspot analysis
  df$outlier[df$side_outlier | df$obvious_outlier] <- T
  
  
  df$hotspot_outlier[df$idPot == i] <- F
  
}


df$outlier[df$side_outlier  | df$obvious_outlier | df$hotspot_outlier] <- T

## AJW old version
# pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#01_correct_gravimetric_outliers.pdf", sep = "_")), width = 10, height = 10)
# for (i in unique(df$idPot)){
#   
#   input <- df[df$idPot == i,]
#   plot(Weight_corr~decimalDay,input,type= "n")
#   points(Weight_corr~decimalDay,input[!input$outlier,])
#   x = input$decimalDay
#   y = input$Weight_corr
#   rectangle(x,y)
#   title(paste0("idPot : ",i))
#   
#   # outliers
#   points(Weight_corr~decimalDay,input[input$side_outlier,],col= "red",pch = 16)
#   points(Weight_corr~decimalDay,input[input$obvious_outlier,],col= "orange",pch = 16)
#   points(Weight_corr~decimalDay,input[input$hotspot_outlier,],col= "violet",pch = 16)
#   
#   legend("topright",legend = c("side","obvious","hotspot"),fill = c("red","orange","violet"),title = "outliers")
#   legend("bottomleft",legend = c("irrigation"),fill = c("black"))
#   
#   if(i %in% unique(corr1$idPot)){
#     abline(v = corr1$decimalDay[corr1$idPot == i])
#   }
# }
# dev.off()

## New version
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#01c_correct_gravimetric_outliers.pdf", sep = "_")), width = 10, height = 5)
for (geno in sort(unique(ColorsTrt$idGenotype)))
  {
  idpots <- genolist$idPot[genolist$idGenotype == geno]
  for (pot in sort(idpots))
    {
    input <- df[df$idPot == pot, ]
    prepare_kin(input, Time_var = "decimalDay", E_var = "Weight_corr",
                ylab = "Corrected pot weight (g)",
                main = paste(geno, " - Pot ", pot, sep = ""),
                inside = F, irrig_show_mode = "pot", pot = pot)
    points(Weight_corr ~ decimalDay, data = input[!input$outlier,], type = "p", cex = 0.5)
    
    # outliers
    points(Weight_corr ~ decimalDay, input[input$side_outlier, ], col = "red", pch = 16, cex = 0.5)
    points(Weight_corr ~ decimalDay, input[input$obvious_outlier, ], col = "orange", pch = 16, cex = 0.5)
    points(Weight_corr ~ decimalDay, input[input$hotspot_outlier, ], col = "violet", pch = 16, cex = 0.5)
    
    legend("topright", legend = c("side", "obvious", "hotspot"), 
           col = c("red", "orange", "violet"), pch = 16, cex = 0.5,
           title = "Outliers", bty = "n")
    
    rm(input)
    }
  }
dev.off()


#------------------------------------------------------------------------------#
#                               Transpiration                                  #
#------------------------------------------------------------------------------#

colnames(df)[match("weight",colnames(df))] <- "initial_weight" # renaming because transpiration function takes "weight" as input.
df$weight <- df$Weight_corr
input = df[!df$outlier,]
# This line takes a few minutes:
dfE_v7 <- Transpi_calc_v7(input = input, freq = 30, min_around = 90, lightsOFF = lightsOFF, nightperiod = Sko_Per, method = "lm",max = 180,nop=2, max_end = 120) # calculate transpiration with function

# add surface
dfE_v7$surface = NA # need to check in the future: 217 start is NA for the surface.. 
for (i in 1:nrow(dfE_v7)){
  dfE_v7$surface[i] = mean(grv$surface[grv$idPot == dfE_v7$idPot[i] & grv$decimalDay>= dfE_v7$min_decimalDay[i] & grv$decimalDay<= dfE_v7$max_decimalDay[i]])
}

dfE_v7$E_mmol_per_m2_s  = dfE_v7$E / (dfE_v7$surface * 10^-6)
# head(dfE_v7)

transpi = merge(dfE_v7,genolist,by="idPot")

# which(diff(aggregate(transpi$idPotManip, by = list(transpi$idPotManip), FUN = length)$x) != 0) 

transpi$meanVPD <- NA

for (i in 1:nrow(transpi)){
  #i=1
  min = transpi$min_decimalDay[i]
  mx = transpi$max_decimalDay[i]
  transpi$meanVPD[i] <- mean(meteo$VPD[meteo$decimalDay > min & meteo$decimalDay < mx],na.rm=T)
}
# head(transpi)
# no VPD data for some parts because of missing values! 

transpi <- transpi[order(transpi$idPot, transpi$decimalDay), ]
transpi$E_mmol_per_m2_s_kPa  <- transpi$E_mmol_per_m2_s  / transpi$meanVPD

## AJW old version
# graphs
# pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#01_transpiration.pdf", sep = "_")), width = 10, height = 10)
# x = transpi$decimalDay
# y = transpi$E_mmol_per_m2_s
# 
# for (i in unique(transpi$idGenotype)){
#   # i = unique(genolist$idGenotype)[1]
#   input <- transpi[transpi$idGenotype == i,]
#   plot(E_mmol_per_m2_s~decimalDay,input,type='n',ylab=ylab_tr,ylim=c(0,4),xlim=range(x,na.rm = T))
#   title(i)
#   rectangle(x,y)
#   
#   for (id in unique(input$idPot)){
#     points(E_mmol_per_m2_s~decimalDay,input[input$idPot == id,],type="l",col=Colors[match(id,unique(input$idPot))])
#     points(E_mmol_per_m2_s~decimalDay,input[input$idPot == id,],type="p",col=Colors[match(id,unique(input$idPot))],cex=0.25)
#     
#   }
#   
#   legend("topright",legend = unique(input$idPot),fill=Colors)
# }
# dev.off()

#------------------------------------------------------------------------------#
#                               Soil Water Content SWC                         #
#------------------------------------------------------------------------------#

# loop on every pot 
# take the SWCs for every pot between the decimaldays used to calculate the transpiration
df <- df[order(df$idPot,df$decimalDay),]

df$SWC <- 0
for (i in unique(df$idPot)){
  # i = "217"
  potweight <- swc$nonPerforatedPotWeight[swc$idPot == i]+swc$perforatedPotWeight[swc$idPot == i]
  drysoilweight <- swc$drySoilWeight[swc$idPot == i]
  df$SWC[df$idPot == i] <- (df$initial_weight[df$idPot == i] - drysoilweight - potweight)/drysoilweight
}

# add the mean SWC to the transpiration file
transpi <- transpi[order(transpi$idPot,transpi$decimalDay),]
transpi$SWC <- 0
for (i in 1:nrow(transpi)){
  # i = 2
  id = transpi$idPot[i]
  transpi$SWC[i] <- mean(df$SWC[df$idPot == id & df$decimalDay > transpi$min_decimalDay[i] & df$decimalDay < transpi$max_decimalDay[i]]) 
}
# head(transpi)

colnames(transpi)[match("Treatment",colnames(transpi))] <- "idWatering" 



#------------------------------------------------------------------------------#
#                          Format time, plot and save                          #
#------------------------------------------------------------------------------#

# Format the dataframe with a uniform time sequence for all pots that starts at 0 at Time_ON0
df <- format_time(transpi, Time_var = "decimalDay", Trt_var = c("idGenotype"), time_step = 30)

## New version of the graphs (with soil water content and VPD)
pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#01d_transpiration_all_pots.pdf", sep = "_")), width = 10, height = 5)
for (geno in sort(unique(df$idGenotype)))
  {
  for (pot in sort(unique(df$idPot[df$idGenotype == geno])))
    {
    dat <- df[df$idPot == pot, ]
    prepare_kin(dat, use_VPD = T, Time_var = "Time", E_var = "E_mmol_per_m2_s_kPa",
                main = paste(geno, " - Pot ", pot, sep = ""),
                inside = F, irrig_show_mode = "pot", pot = pot,
                add_SWC = T, mar = c(2.5, 3.5, 2.5, 7.5))
    points(E_mmol_per_m2_s_kPa ~ Time, data = dat, type = "o", cex = 0.5)
    rm(dat)
    }
  }

dev.off()


pdf(file.path(dir_Exp, "Figures", paste(idExp, "Step#01e_transpiration_all_pots_per_condition.pdf", sep = "_")), width = 10, height = 5)
n.max <- max(aggregate(df$idPot[!duplicated(df$idPot)], by = list(df$idGenotype[!duplicated(df$idPot)]), FUN = length)$x) # maximum number of replicates per treatment
for (geno in sort(unique(df$idGenotype)))
  {
  df.geno <- df[df$idGenotype == geno, ]
  prepare_kin(df.geno, use_VPD = T, add_VPD = T, Time_var = "Time", E_var = "E_mmol_per_m2_s_kPa",
              ylim = c(-0.5, 4),
              main = geno,
              inside = F, irrig_show_mode = "mean")
  color = 0
  for (pot in sort(unique(df.geno$idPot)))
    {
    #require(scales)
    color <- color + 1
    dat <- df[df$idPot == pot, ]
    points(E_mmol_per_m2_s_kPa ~ Time, data = dat, type = "o", col = hue_pal()(n.max)[color], cex = 0.5)
    rm(dat)
    }
  legend("top", ncol = 4, bty = "n",
         legend = paste("Pot ", sort(unique(df.geno$idPot[df.geno$idGenotype == geno]))),
         col = hue_pal()(n.max)[1:color], lty = 1, pch = 21, pt.cex = 0.5)
  rm(df.geno)
  }

dev.off()


p2f <- file.path(dir_Exp, "Processed_data", paste(idExp, "pot_transpiration.csv", sep = "_"))
write.csv(df,p2f,row.names = F)

# cleaned gravimetrical data
p2f <- file.path(dir_Exp, "Processed_data", paste(idExp, "cleaned_gravi_data.csv", sep = "_"))
write.csv(df,p2f,row.names = F)
