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
meteo$date <- as.POSIXct(strptime(meteo$date,format= "%Y-%m-%d %H:%M:%S", tz = "UTC")) # now the date/hour is not changing
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
grv$date <- as.POSIXct(strptime(grv$date,format= "%Y-%m-%d %H:%M:%S", tz = "UTC")) 
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

for (geno in unique(ColorsTrt$idGenotype)){
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
#                               Cleaning                                       #
#------------------------------------------------------------------------------#

#---------------------------STEP 1---------------------------------------------#
# detect and correct for manual water irrigations
# check for outliers if still in 
if(nrow(grv[c(grv$idPot == "C2M47-321" & grv$ID %in% c(220,221)),])>0){grv <- grv[-c(34508,34509),] }# eliminate two points of one pot because two points are outliers

output <- Rehy_Corr_v4(input = grv,gap = 1)

df_rehy <- output$output # corrected output
corr1 <- output$corr1 # 1 run: detected rehydrations

pdf(file.path(getwd(), spcs, idExp, "Figures", paste(idExp, "Step#01_correct_irrigation.pdf", sep = "_")), width = 10, height = 10)
x = c(grv$decimalDay,df_rehy$decimalDay)
y = c(grv$weight,df_rehy$Weight_corr)
xlim = range(x)
ylim = range(y)

for (i in inpots){
  # i = idpots[3]
  input <- df_rehy[df_rehy$idPot == i,]
  xlim = range(input$decimalDay)
  ylim = range(input$weight,input$Weight_corr)
  plot(weight~decimalDay,input[input$idPot == i,],type='n',ylim = ylim, xlim=xlim, ylab="weight (g)",xlab="decimalDay",main=i)
  rectangle(x,y)
  points(weight~decimalDay,input[input$idPot == i,],type="p")
  
  if(i %in% unique(corr1$idPot)){
    # for each unique was an irrigation
    strt =  corr1$from_id[corr1$idPot == i][1] # take the first one!
    
    points(Weight_corr~decimalDay,df_rehy[df_rehy$idPot == i & df_rehy$ID > strt,],col="red")
    
  }
  legend("bottomleft",legend = c("corrected data"),fill="red")
}
dev.off()

df = df_rehy
# manual correction irrigation
i ="C2M47-374"
plot(Weight_corr~decimalDay,df[df$idPot =="C2M47-374",],main=i)
df$Weight_corr[df$idPot == "C2M47-374" & df$ID > 38] <- df$Weight_corr[df$idPot == "C2M47-374" & df$ID > 38] - 0.28
plot(Weight_corr~decimalDay,df[df$idPot =="C2M47-374",],main=paste0(i," -corrected manually"))

  # corr 2:
i ="C2M47-306"
plot(Weight_corr~decimalDay,df[df$idPot ==i,],main=paste0(i))
df$Weight_corr[df$idPot == i & df$ID %in% c(290:296)] <- df$Weight_corr[df$idPot == i &  df$ID %in% c(290:296)] - 0.35
plot(Weight_corr~decimalDay,df[df$idPot ==i,],main=paste0(i," -corrected manually"))

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

# pots
#idpots <- unique(df$idPot)
# idpots <- unique(df$idPot)[1]
inpots = unique(df$idPot)
for (i in inpots){
  # i = unique(df$idPot)[1]
  # i = "C2M47-217"
  # i =inpots[1]
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
  
  for (j in hotspots){
    # check which points are real outliers at this timepoint
    # outliers detection based on input on the fit. Per hotspot only 1 point can be detected
    
    # j = hotspots[8]
    # j = 47
    
    new_outliers <- Outliers_v5(time = df$decimalDay[df$idPot == as.character(i)], # time input
                                weight = df$Weight_corr[df$idPot == as.character(i)], # weight input
                                hotspot = j, # ID of the hotspot
                                lightPeriod = df$lightPeriod[df$idPot == as.character(i)],
                                outlier = df$outlier[df$idPot == as.character(i)],
                                startdark = startdark,
                                period = darkperiod,
                                daynight = "yes",
                                min_around = 90)
    # df$outlier[df$idPotManip == i & df$ID %in% c(new_outliers) ] <- T
    outliers <- c(outliers,new_outliers)
  }
  
  # outliers <- c(outliers,out_sides)
  # outliers <- c(outliers)
  
  df$hotspot_outlier[df$idPot == i & df$ID %in% outliers ] <- T
  
}


df$outlier[df$side_outlier  | df$obvious_outlier | df$hotspot_outlier] <- T
df_lessstrict <- df

pdf(file.path(getwd(), spcs, idExp, "Figures", paste(idExp, "Step#01_correct_gravimetric_outliers.pdf", sep = "_")), width = 10, height = 10)
for (i in inpots){
  
  input <- df_lessstrict[df_lessstrict$idPot == i,]
  plot(Weight_corr~decimalDay,input,type= "n")
  points(Weight_corr~decimalDay,input[!input$outlier,])
  x = input$decimalDay
  y = input$Weight_corr
  rectangle(x,y)
  title(paste0(i," less_strict"))
  
  # outliers
  points(Weight_corr~decimalDay,input[input$side_outlier,],col= "red",pch = 16)
  points(Weight_corr~decimalDay,input[input$obvious_outlier,],col= "orange",pch = 16)
  points(Weight_corr~decimalDay,input[input$hotspot_outlier,],col= "violet",pch = 16)
  
  legend("topright",legend = c("side","obvious","hotspot"),fill = c("red","orange","violet"),title = "outliers")
  legend("bottomleft",legend = c("irrigation"),fill = c("black"))
  
  if(i %in% unique(corr1$idPot)){
    abline(v = corr1$decimalDay[corr1$idPot == i])
  }
}
dev.off()
#------------------------------------------------------------------------------#
#                           Output                                             #
#------------------------------------------------------------------------------#



