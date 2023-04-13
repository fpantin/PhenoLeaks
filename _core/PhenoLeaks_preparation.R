################################################################################
#                                                                              #
#                             PhenoLeaks - Data preparation                    #
#                                                                              #
#                 Functions to add leaf surface information,                   #
#             clean raw data for irrigations and gravimetric outliers          #
#                                                                              #
#                             Adriaan Westgeest, 2023                          #
#                                                                              #
################################################################################


#------------------------------------------------------------------------------#
#                       define specific parameters for data preparation        #
#------------------------------------------------------------------------------#

OFF <- as.numeric(strftime(as.POSIXct(strptime(lightsOFF,format= "%H:%M:%S"), tz = "UTC"),format="%H",tz = "UTC"))/(24) +
  as.numeric(strftime(as.POSIXct(lightsOFF,format= "%H:%M:%S", tz = "UTC"), format = "%M", tz='UTC'))/(24*60)

darkperiod <- (24-Pho_Per) / 24 
startdark <- as.numeric(round(OFF,10))
enddark <- startdark + darkperiod
enddark <- as.numeric(sub(".*\\.","0.",sprintf("%.10f", round(enddark,10)) ))
color <- c(rgb(0, 32, 96, alpha=75, maxColorValue=255),rgb(255, 192, 0, alpha=70, maxColorValue=255)) # night/day

# test whether "night" is during the normal night or inversed (during the day)
if (startdark + darkperiod < 1){
  enddark = startdark + darkperiod  # day and night inversed
  night = "inversed"
}else{
  enddark = startdark + darkperiod -1
  night = "normal"
}

ylab_tr = expression(paste("E"["rosette"], " (mmol m"^-2, " s"^-1, ")"))
ylab_tr_vpd = expression(paste("E"["rosette"], " (mmol m"^-2, " s"^-1, " KPa"^-1,")"))

#------------------------------------------------------------------------------#
# color vector to distinguish between individual pots
Colors <- c("brown3", "darkblue","darkorange", "cadetblue3", "red","darkolivegreen","chartreuse","darkgrey","chocolate","burlywood","darkkhaki", "darkorchid","deeppink","darkturquoise","darkred","orange")

#------------------------------------------------------------------------------#
#             Function to prepare gravimetric data                             #
#------------------------------------------------------------------------------#

prep_gravi <- function(){
  grv <- as.data.frame(splitstackshape::cSplit(grv, "idPotManip", sep="-",type.convert = F)) # split colum with idPotManip to get idPot, need to add as.data.frame to bring back from table to dataframe.
  colnames(grv)[match("idPotManip_2",colnames(grv))] <- "idPot" # rename column
  colnames(grv)[match("After.watering.weight",colnames(grv))] <- "weight" # rename column
  grv <- grv[,c("idPot","date","weight")]
  
  grv$date <- as.POSIXct(strptime(grv$date,format= "%Y-%m-%d %H:%M:%S", tz = "UTC")) 
  grv$decimalDay <- decimalDay(column=grv$date) # add decimal day
  grv$hourday <- as.numeric(sub(".*\\.","0.",sprintf("%.10f", round(grv$decimalDay,10)) )) # add decimal hour on day
  grv$dayofyear = as.integer(grv$decimalDay) # day of the year
  
  # filtering
  grv <- grv[!is.na(grv$weight),]# delete all NA's
  if(anyDuplicated(grv) > 0){ grv <- grv[!duplicated(grv),] }# remove identical lines
  grv <- grv[which(grv$date > from & grv$date < to),] #select the data from the start and end of gravimetric experiment
  
  # take out days with less than 4 measurements on a single day. Alternive is selecting on "genolist$idPot" (future)
  grv$idPotManipday <-  paste(grv$idPot, grv$dayofyear, sep="-") # new vector of idPotmanip 
  grv$idPotManipday <- as.factor(grv$idPotManipday)
  numreps <- as.data.frame(table(grv$idPotManipday))
  removepotsday <- unique(numreps[numreps$Freq < 4,]$Var1)
  grv <- grv[!grv$idPotManipday %in% removepotsday,] # take these pots out
  grv <- grv[grv$idPot %in% genolist$idPot,] # alternative for selecting only pots with more than 4 measurements per day. Better.
  
  # adding more information
  grv <- grv[order(as.numeric(grv$idPot),grv$decimalDay),] # order per pot and decimalDay
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
  grv <- grv[,!colnames(grv) %in% c("idPotManipday")]
  grv <- grv[order(as.numeric(grv$idPot),grv$decimalDay),]
  return(grv)
}

#------------------------------------------------------------------------------#
#             Function to calculate decimalDay from Posix data                 #
#------------------------------------------------------------------------------#
# can calculate the decimal day of date data but also datehour data

decimalDay <- function(column){
  # column needs to be a POSIX specified vector
  output <- as.numeric(strftime(column, format = "%j", 
                                          tz='UTC')) + (as.numeric(strftime(column, format = "%H", tz='UTC'))/(24)) +
    (as.numeric(strftime(column, format = "%M", tz='UTC'))/(24*60)) +
    (as.numeric(strftime(column, format = "%S", tz='UTC'))/(24*60*60)) 
  
  return(output)
}

#------------------------------------------------------------------------------#
#                       Add rectangle to figures                               #
#------------------------------------------------------------------------------#

rectangle <- function(x,y){
  # function to add rectangles to the plot
  days <- seq(min(unique(as.integer(x))),max(unique(as.integer(x))),1) 
  days = c(days,(min(days)-1),(max(days)+1))
  
  for (d in days){
    # i = 179
    rect(xleft= d + startdark, xright=d + startdark + darkperiod, ybottom= min(y,na.rm = T)-min(y,na.rm = T), ytop= max(y,na.rm = T)*2, density= NULL, col= color[1], border = NA)
    
  }
  
}

#------------------------------------------------------------------------------#
#             Function to test whether a vector with decimalDays               #  
#                              is at night or day                              #
#------------------------------------------------------------------------------#


dd_light <- function(dd){
  # works only with hours, not with day data
  
  hourday <- as.numeric(sub(".*\\.","0.",sprintf("%.10f", round(dd,10)) )) # extract only hour/min/sec data 
  
  lightp <- rep("light",length(hourday)) # create vector with information about the lightperiods
  
  if (startdark + darkperiod > 1){
    # night during night
    
    lightp[hourday > startdark | hourday <= enddark] <- "dark"
  }
  
  if (startdark + darkperiod < 1){
    # night during day
    
    lightp[hourday > startdark & hourday <= enddark] <- "dark"
    
  }
  
  return(lightp)
  
}

#------------------------------------------------------------------------------#
#             Function to test the best fit model to describe plant growth     #
#------------------------------------------------------------------------------#

# This function need a dataframe including "idPot", "outlier" (identified before manually), 
# "Areamm2", "decimalDay".

surf_fit <- function(surface){
  
  surf_coef <- data.frame(idPot = NULL, model = NULL,growth = NULL, from = NULL,
                          to = NULL, intercept = NULL, slope = NULL)
  
  mod = "test"
  growth = "continu"
  
  for (i in unique(surface$idPot)){
    # i = "C2M43-93"
    
    model = mod
    input = surface[surface$idPot == i & !surface$outlier,]
    # plot(Areamm2~decimalDay, data = input)
    
    
    if (nrow(input) > 1 & growth == "continu"){
      
      if(model == "lm"){
        lm <- lm(Areamm2~decimalDay, data = input)
      }
      if(model == "log"){
        lm <- lm(log(Areamm2)~ decimalDay, data = input)
      }
      if(model == "test"){
        
        lm_log <- lm(log(Areamm2)~ decimalDay, data = input)
        lm_lm <- lm(Areamm2~decimalDay, data = input)
        
        if (summary(lm_log)[[8]] > summary(lm_lm)[[8]]){
          lm <- lm_log
          model = "log"
          
        }else{
          lm <- lm_lm
          model = "lm"
        }
      }
      
      surf_coef <- rbind(surf_coef,data.frame(idPot = i, model = model,growth = growth, from = min(input$decimalDay), to = max(input$decimalDay), intercept = as.numeric(lm$coefficients[1]), slope = as.numeric(lm$coefficients[2])))
      
    }
  }
  
  return(surf_coef)

}


#------------------------------------------------------------------------------#
#             Function to test the best fit model to describe plant growth     #
#------------------------------------------------------------------------------#

surface_add <- function(input){
  input$surface <- 0
  no_surf <- c()
  for (i in unique(input$idPot)){
    
    if(i %in% unique(surf_coef$idPot)){
      i = as.character(i)
      
      mod <- surf_coef$model[surf_coef$idPot == i & surf_coef$growth == "continu"]
      
      slope = surf_coef$slope[surf_coef$idPot == i & surf_coef$growth == "continu"]
      
      intercept = surf_coef$intercept[surf_coef$idPot == i & surf_coef$growth == "continu"]
      
      if (mod == "lm"){
        
        input$surface[input$idPot == i] <- intercept +  input$decimalDay[input$idPot == i] * slope
      }
      if (mod == "log"){
        input$surface[input$idPot == i] <- exp(intercept +  input$decimalDay[input$idPot == i] * slope)
      }
    }else{
      input$surface[input$idPot == i] <- NA
      no_surf <- c(no_surf,i)
    }
    
  }
  #print(unique(no_surf))
  #print(length(unique(no_surf)))
  
  return(input)
}

#------------------------------------------------------------------------------#
#             Function to detect and correct for irrigations                   #
#------------------------------------------------------------------------------#

Rehy_Corr_v4 <- function(input,gap,jumps="positive",method = "automatic"){
  # if jumps is positive, correct only for irrigation (so not other jumps with negative values)
  # if jumps is both, takes into account negative and positive "irrigations"
  
  # method = "automatic" : automatic detection with the script based on gap. 
  # method = c() : insert the end dates of the cycle preceding the (manual) irrigation. 
  # works only when the irrigations happened at the same time for all pots!!!! Not the case for C2M47 for example
  
    #--------------------------------------------#
  # Detection
  selecpoint <- function(timer){
    
    middledater = input$date[j] # get the time
    # select values based on certain time frame
    beforedater <- middledater - 1*60* timer # add .. hour back in time
    afterdater <- middledater + 1*60* timer # add .. hour
    # select index values 
    dataindex <- allpoints[input$date[allpoints] >= beforedater & input$date[allpoints] <= afterdater]
    
    middledark = as.character(input$lightPeriod[j]) # get the lightperiod 
    # reselect the index base on light/dark (comparison only made with similar light/dark based on frequency)
    inputdark <- input$lightPeriod[dataindex] 
    
    if(middledark == "darklight"){
      dataindex1 <- dataindex[inputdark %in% c(input$lightPeriod[j-1],"darklight")] # take the lightperiod which was before! 
    }else{
      dataindex1 <- dataindex[inputdark %in% c(middledark,"darklight")]
    }
    
    dataindex2 <- dataindex1[dataindex1 > j] # and after the point of the difference
    info <- "after"
    
    if (length(dataindex2) < 2){ # when the irrigation is at the day night transition and middledark =! as points afterwards
      beforedater <- middledater - 1*60* timer*2 # add .. hour back in time
      afterdater <- middledater + 1*60* timer*2 # add .. hour
      dataindex <- allpoints[input$date[allpoints] >= beforedater & input$date[allpoints] <= afterdater]
      middledark = as.character(input$lightPeriod[j]) # get the lightperiod 
      # reselect the index base on light/dark (comparison only made with similar light/dark based on frequency)
      inputdark <- input$lightPeriod[dataindex] 
      if(middledark == "darklight"){
        dataindex1 <- dataindex[inputdark %in% c(input$lightPeriod[j-1],"darklight")] # take the lightperiod which was before! 
      }else{
        dataindex1 <- dataindex[inputdark %in% c(middledark,"darklight")]
      }
      dataindex2 <- dataindex1[dataindex1 > j] # and after the point of the difference
    }
    
    if (length(dataindex2) < 2){ # still smaller, means that in the same photoperiod no points
      beforedater <- middledater - 1*60* timer # add .. hour back in time
      afterdater <- middledater + 1*60* timer # add .. hour
      dataindex <- allpoints[input$date[allpoints] >= beforedater & input$date[allpoints] <= afterdater]
      dataindex2 <- dataindex[dataindex < j] # and before the point of the difference, NOT the point itself to avoid outlier sensibility
      info <- "before"
      }
    return(c(dataindex2,info))
  }
  
  if(!method[1] == "automatic"){
    # when we have irrigation dates
    method_dd <- decimalDay(method) # transform to decimaldays
  }
    
  rehydata <- data.frame(idPot = NULL, 
                         weight_before = NULL,
                         weight_after = NULL,
                         predicted_water_add = NULL,
                         from_index = NULL,
                         decimalDay = NULL,
                         note = NULL)
  
  for (id in unique(input$idPot)){
    # id = "217"
    # id = unique(input$idPot)[1]
    allpoints <- which(input$idPot == id) # select the index numbers (not ID col!). 
    
    if(method[1] == "automatic"){
      
      selec <- allpoints[-length(allpoints)] # take out last one if automatic
    }
    
    if(!method[1] == "automatic"){
      # when we have irrigations dates
      
      selec <- c()
      for (i in 1:length(method)){
        # for every irrigation, select the closest point before the irrigation.
        # i = 1
        selec <- c(selec,which(input$idPot == id & input$decimalDay < method_dd[i])[length(which(input$idPot == id & input$decimalDay < method_dd[i]))]) # select the index numbers (not ID col!). 
      }
    }
    
    #selec = selec[selec > 161]
    outlier <- c()
    for (j in selec){
      
      # for every measurement for this pot...
      # j are index numbers!! not ID numbers
      # j = 22
      if(!j %in% outlier){ # this is when there is one obvious outlier for which we do not calculate the difference for the next point
        beforeweight = input$weight[j]
        afterweight = input$weight[j + 1]
        diff1 = beforeweight -  afterweight
        # its a real rehydration if the diff between the beforeweight and multiple afterweights are bigger then the gap. 
        # diff2 = beforeweight -  input$weight[j - 1] # control to eliminate beforehand outliers. 
        
        if (abs(diff1) > gap){
          # if difference bigger then gap, can be an outlier or rehydration
          diff2 = abs(beforeweight -  input$weight[j + 2])
          diff3 = abs(beforeweight -  input$weight[j + 3])
          
          if(diff2 > gap & diff3 > gap){
            # if two points similar
            # if difference of before weight with other points afterwards is also bigger, then a rehydration
            
            # now need to know if j is the right one to start with. Most of the times there are bad points during the rehydration look to this.
            # How many points till normal again? Define new j
            # future: test also if the initial j is correct (may be already affected)
            
            # next point is good point 
            out <- selecpoint(timer = 160) #
            if(length(out)>1){
              # when points are found
              dataindex2 <- as.numeric(out[1:(length(out)-1)])
              info <- out[length(out)]
            }else{dataindex2 <- integer()}
            
            if (length(dataindex2) >= 2){
              
              inwei <- input$weight[dataindex2]
              indat <- input$decimalDay[dataindex2]
              
              ols1 <- lm(inwei ~ indat)
              
              if(info == "after"){ # points used after j (classic)
                predat <- input$decimalDay[j] # for the date j
                outpred <- predict(ols1,data.frame(indat = predat))
                differ <- outpred - beforeweight # compare prediction with real value for j
              }
              if(info == "before"){# points used before j (when not enough points in same photoperiod)
                predat <- input$decimalDay[j+1] # for the date j+1
                outpred <- predict(ols1,data.frame(indat = predat))
                differ <- afterweight - outpred # compare prediction with augmented value for j
              }
              if(differ > gap){ 
                rehydata <- rbind(rehydata, data.frame(idPot = id, 
                                                       weight_before = beforeweight,
                                                       weight_after = afterweight,
                                                       predicted_water_add = differ,
                                                       from_id = input$ID[j],
                                                       decimalDay = predat,
                                                       note = info))
              }
              # rm(differ,outpred,predat,ols1,indat,inwei)
            }
            
            
          }else{outlier <- c(outlier,j+1)}
          # rm(beforeweight,afterweight)
        }
      }
      
      
      
    }
    
  }
  
  
  #--------------------------------------------#
  # Correction
  input$Weight_corr <- input$weight
  if (jumps == "positive"){ # only real irrigations (increase in weight)
    rehydata <- rehydata[rehydata$predicted_water_add >0,]
  }
  
  for (i in 1:length(rehydata$idPot)){
    # i = 1
    idpotm <- as.character(rehydata$idPot[i])
    ind <- rehydata$from_id[i]
    wadd <- rehydata$predicted_water_add[i]
    
    input[which(input$idPot == idpotm & input$ID > ind),]$Weight_corr <-  input[which(input$idPot == idpotm & input$ID > ind),]$Weight_corr - wadd
    
    rm(idpotm,ind,wadd)
  }
  
  
  #--------------------------------------------#
  # Graphs
  mylist <- list("output"=input,"corr1" =rehydata)
  return(mylist)
  
}

#------------------------------------------------------------------------------#
#             Function to detect gravimetric outliers                          #
#------------------------------------------------------------------------------#

Outliers_v4_2 <- function(time,weight,startdark,darkperiod,min_dark,max_dark,min_light,max_light,surf_i,dd_surf_i = NULL){
  
  #----------------------------------------------------------------------#
  #----------------------explanation-------------------------------------#
  #----------------------------------------------------------------------#
  
  
  # This is an R script to analyse the outliers of a sequence of weight on ONE pot.
  # The script will return a vector with ID's for the type of outlier.
  # The method is based on biological limits of the difference between two points
  
  # inputs:
  #     - time: decimalDay vector 
  #     - weight: weight vector in g
  #     - startdark: start of the darkperiod in decimaldays
  #     - period: hours of dark period in decimaldays
  #     - min_dark: lower limit of transpiration during the dark period 
  #     - max_dark: upper limit of transpiration during the dark period 
  #     - min_light: lower limit of transpiration during the light period 
  #     - max_light: upper limit of transpiration during the light period 
  #     - surface: surface in mm2, if more surfaces (growth) then is a vector
  
  ###### NEW in version 4_2 in comparison with v4: adapted to Arabidopsis
  #     - more then 1 surface 
  #     - more then 3 lightperiods (first only 3, starting with day!)
  
  # to be improved:
  #     - first point of a period always taken out so no outlier detection on this one. 
  #----------------------------------------------------------------------#
  #----------------------environment-------------------------------------#
  #----------------------------------------------------------------------#
  
  # Fixed parameter:
  M <- 18.015 # molecular weight of water (g/mol)
  if (startdark + darkperiod < 1){
    enddark = startdark + darkperiod  # period is darkperiod!
    night = "inversed"
  }else{
    enddark = startdark + darkperiod - 1 # in hours
    night = "normal"
  }
  
  dayperiod <- 1 - darkperiod
  
  # Functions
  dd_light <- function(dd){
    # work only with hours!! not with day data
    
    hourday <- as.numeric(sub(".*\\.","0.",sprintf("%.10f", round(dd,10)) )) # decimal time on the day
    
    lightp <- rep("light",length(hourday))
    
    if (startdark + darkperiod > 1){
      # night during night
      
      lightp[hourday > startdark | hourday <= enddark] <- "dark"
    }
    
    if (startdark + darkperiod < 1){
      # night during day
      
      lightp[hourday > startdark & hourday <= enddark] <- "dark"
      
    }
    
    
    return(lightp)
    
  }
  
  
  #----------------------------------------------------------------------#
  #----------------------start-------------------------------------------#
  #----------------------------------------------------------------------#
  
  # Put the time vector in the good order
  input <- data.frame(time,weight)
  input <- input[order(input$time),]
  
  # add an ID colum 
  input$ids <- seq(1,nrow(input),by = 1)
  
  # slopes between points
  slopes <- diff(input$weight) / diff(input$time) 
  stopifnot(nrow(input)-1 == length(slopes))
  
  ids <- input$ids[-1]
  stopifnot(nrow(input)-1 == length(ids))
  
  transpi_times <- input$time[-nrow(input)] + (diff(input$time)/2)
  stopifnot(nrow(input)-1 == length(transpi_times))
  
  input$units <- (M * (10^-3) * 24* 60 * 60 * surf_i * 10^-6)
  
  units = input$units[-1]
  
  transpi <- slopes / units
  stopifnot(nrow(input)-1 == length(transpi))
  
  transpi_lightPeriod <- dd_light(dd= transpi_times)
  stopifnot(nrow(input)-1 == length(transpi_lightPeriod))
  
  cdata <- data.frame(ids,slopes,transpi,transpi_times,transpi_lightPeriod)
  
  # first point selection: what is the first good point?
  # limits depending if period starts in dark or in light!
  first_point <- c(1:3)
  # define limits for first points
  if(all(transpi_lightPeriod[first_point] == "light") | all(transpi_lightPeriod[first_point] == "dark")){
    # if all the first points are in 1 lightperiod either dark of light
    if(all(transpi_lightPeriod[first_point] == "light")){
      min_lim <- min_light
      max_lim <- max_light
    }
    if(all(transpi_lightPeriod[first_point] == "dark")){
      min_lim <- min_dark
      max_lim <- max_dark
    }
    
  }else{
    # if first points are a mix of 2 lightperiods, take maximum limits
    min_lim <- min_dark
    max_lim <- max_light
  }
  
  while( max(first_point) <= max(input$ids) - 3 &
         all(transpi[first_point] > -max_lim & transpi[first_point] < -min_lim)== F){
    # check when for the first time 3 points in a row show a good transpiration
    # diff(c(as.numeric(quantile(transpi[first_point])[2]),as.numeric(quantile(transpi[first_point])[4]))) > 0.5
    first_point <- first_point + 1 
  }
  
  # last point selection: what is the last good point?
  last_point <- max(ids) - c(3:1)
  if(all(transpi_lightPeriod[last_point] == "light") | all(transpi_lightPeriod[last_point] == "dark")){
    # if all the first points are in 1 lightperiod either dark of light
    if(all(transpi_lightPeriod[last_point] == "light")){
      min_lim <- min_light
      max_lim <- max_light
    }
    if(all(transpi_lightPeriod[last_point] == "dark")){
      min_lim <- min_dark
      max_lim <- max_dark
    }
    
  }else{
    # if first points are a mix of 2 lightperiods, take maximum limits
    min_lim <- min_dark
    max_lim <- max_light
  }
  
  while(min(last_point) >= min(input$ids)+3 & all(transpi[last_point] > -max_lim & transpi[last_point] < -min_lim +0.1)== F ){
    # check when for the last time 3 points in a row show a good transpiration
    # quantile accepts pots that are variable but stable at the beginning
    # diff(c(as.numeric(quantile(transpi[last_point])[2]),as.numeric(quantile(transpi[last_point])[4]))) > 0.5
    last_point <- last_point - 1 
  }
  
  first_point <- first_point[1] 
  last_point <- last_point[3]
  
  cdata <- cdata[cdata$ids > first_point & cdata$ids <= last_point,]
  # tail(cdata)
  ids_final_out <- c()
  
  
  ###################################################### 
  
  if (first_point + 3 > last_point){
    # in the case that no stable points are found from start to end (possibly only 2 are found)
    ids_final_out <- c()
    specific_outliers <- c()
    obvious_outliers_final <- c()
    
  }else{
    # based on stable first point
    # define periods
    days <- unique(as.integer(input$time))
    days <- days[order(days)]
    cdata$day <- as.integer(cdata$transpi_times)
    
    ### find out the start of the first period !!
    if(night == "inversed"){
      # start with the first day
      i = unique(cdata$day)[1] # not take days because we can only have points in second day (other outliers)
      if(!any(cdata$transpi_times[cdata$day ==  i] - i < startdark)){
        # test which period we find points to characterise the light period, if any go to else
        if(!any(cdata$transpi_times[cdata$day ==  i] - i > startdark & cdata$transpi_times[cdata$day ==  i] - i < enddark)){
          # if any found go to else
          if(!any(cdata$transpi_times[cdata$day ==  i] - i > enddark)){
            # if any found go to else
            warning(paste0("points are found but period not identified in pot "))
          }else{
            startperiod <- i + enddark
            endperiod <- i + enddark + dayperiod
          }
        }else{
          startperiod <- i + startdark
          endperiod <- i + enddark
        }
      }else{
        # if some points are found in the first period, this period will be period 1
        startperiod <- i
        endperiod <- i + startdark
      }
    }
    
    if(night == "normal"){
      # start with the first day
      i = unique(cdata$day)[1] # not take days because we can only have points in second day (other outliers)
      # determine in which period we start on the first day of the whole kinetic, output is the limits of this period
      if(!any(cdata$transpi_times[cdata$day ==  i] - i < enddark)){
        # Does it start in the morning dark period ?
        if(!any(cdata$transpi_times[cdata$day ==  i] - i > enddark & cdata$transpi_times[cdata$day ==  i] - i < startdark)){
          # does it start in the day light period ?
          if(!any(cdata$transpi_times[cdata$day ==  i] - i > startdark)){
            # Does it start in the evening night period ?
            warning(paste0("points are found but period not identified in pot "))
          }else{
            startperiod <- i + startdark
            endperiod <- i + startdark + darkperiod
          }
        }else{
          startperiod <- i + enddark
          endperiod <- i + startdark
        }
      }else{
        # if some points are found in the first period, this period will be period 1
        startperiod <- i
        endperiod <- i + enddark
      }
    }
    
    # now we have the start and end of the first period! From this starting point we can start a loop for every next lightperiod and put the ids of outliers in a list vector
    obvious_outliers <- c()
    specific_outliers <- c()
    transition_ids <- c()
    # ids_out <- c()
    
    
    while(max(cdata$transpi_times) > startperiod){
      # loop continues when there are still points outside a period.
      pointss <- cdata$ids[cdata$transpi_times > startperiod & cdata$transpi_times <= endperiod]
      transition_ids <- c(transition_ids,pointss[1]) # first one of a period, closure at start of night or increase start day, number of points dependent on the frequency of measurement and fast or slow closure/opening. These points will be later examined with less strict limits.
      select_ids <- pointss[-1] # take always out the first one because is acclimatization of the period! 
      
      
      if(length(select_ids) > 1){
        # more than 1 points in a specific period
        # use quantile to take into account the variability potentially observed within a pot. 
        select_period <- unique(cdata$transpi_lightPeriod[cdata$ids %in% select_ids])
        select_values_raw <- cdata$transpi[cdata$ids %in% select_ids]
        q1 = quantile(select_values_raw)[2]
        q3 = quantile(select_values_raw)[4]
        iqr = IQR(select_values_raw) # diff q1 and q3
        max_quantile =  (q1 - 1.5*iqr) # 1.5 factor to detect outliers https://towardsdatascience.com/why-1-5-in-iqr-method-of-outlier-detection-5d07fdc82097
        min_quantile =  (q3 + 1.5*iqr)
        
        if (select_period == "dark"){
          # obvious outliers
          obvious_outliers <- c(obvious_outliers,select_ids[(select_values_raw < -max_dark | select_values_raw > -min_dark)])
          specific_outliers <- c(specific_outliers,select_ids[(select_values_raw < max_quantile | select_values_raw > min_quantile)])
        }
        if (select_period == "light"){
          # obvious outliers
          obvious_outliers <- c(obvious_outliers,select_ids[(select_values_raw < -max_light | select_values_raw > -min_light)])
          specific_outliers <- c(specific_outliers,select_ids[(select_values_raw < max_quantile | select_values_raw > min_quantile)])
        }
      }
      
      # define new limits of the next period
      if (endperiod == as.integer(endperiod) + enddark){
        # if it is the end of the night
        startperiod <- endperiod
        endperiod <- endperiod + dayperiod
      }else{
        # if it is the end of day
        startperiod <- endperiod
        endperiod <- endperiod + darkperiod
      }
      
    }
    
    transition_ids_outliers <- cdata$ids[cdata$ids %in% transition_ids  & (cdata$transpi < -max_light | cdata$transpi > -min_dark)] # limit of both dark and night because transition
    
    
    obvious_outliers_final <- c()
    # obvious_outliers_final <- obvious_outliers
    
    # check here if it is really the second point that is the outlier or that it is the first point.
    # normally it is the second point because if not the point before would be already detected. 
    # this function also checks when two consecutive outliers are detected the second one is the result of the other first one.
    # only need to check how the ids are correctly identified. OK old ids are taken into account!
    if(length(obvious_outliers)> 0){
      # check obvious ones (always taken out directly!! So need to be sure)
      obvious_outliers <- obvious_outliers[order(obvious_outliers)]
      notsoobvious_outliers <- c()
      
      while(length(obvious_outliers) > 0 ){
        # to check if by taking out the point before the outlier, the outlier is still an outlier. Do this is a loop that if the point before is an outlier, this one is taken out and compared to the good one before.
        i = obvious_outliers[1]
        # print(i)
        realouttest_ids <- c(obvious_outliers_final, i - 1) # consisting out of outliers detected before in this loop and the point before!
        input1 <- input[!input$ids %in% c(realouttest_ids) & input$ids >= first_point & input$ids <= last_point + 1,] # redefine the df without the outlier
        
        slopes <- diff(input1$weight) / diff(input1$time)
        ids <- input1$ids[-1] # same ids as before
        units <- input1$units[-1]
        
        transpi <- slopes / units
        transpi_times <- input1$time[-nrow(input1)] + (diff(input1$time)/2)
        transpi_lightPeriod <- dd_light(dd= transpi_times)
        cdatab <- data.frame(ids,slopes,transpi,transpi_times,transpi_lightPeriod)
        cdatab <- cdatab[cdatab$ids >= first_point & cdatab$ids <= last_point,]
        
        dark_ids <- cdatab$ids[!cdatab$ids %in% transition_ids & cdatab$ids >= i & cdatab$transpi_lightPeriod == "dark" & (cdatab$transpi < -max_dark | cdatab$transpi > -min_dark )] # redefine other outliers. 
        
        light_ids <- cdatab$ids[!cdatab$ids %in% transition_ids & cdatab$ids >= i & cdatab$transpi_lightPeriod == "light" & (cdatab$transpi < -max_light | cdatab$transpi > -min_light)]# redefine other outliers. 
        
        obvious_outliers <- c(dark_ids,light_ids)
        obvious_outliers <- obvious_outliers[order(obvious_outliers)]
        # obvious_outliers <- obvious_outliers + 1 # take into account the shift in ids when one ID is taken out.
        
        if(i %in% obvious_outliers){
          # indicates that this one really is an outlier and not because of the first one
          obvious_outliers_final <- c(obvious_outliers_final,i)
          obvious_outliers <- obvious_outliers[!obvious_outliers %in% i] # take out the 'real' bad one so that the loop continues with the rest!
          obvious_outliers <- obvious_outliers[order(obvious_outliers)]
        }else{
          # less strong indication that it is really an outlier: how to test this?
          notsoobvious_outliers <- c(notsoobvious_outliers,i)
        }
        
        # now need actually a step that if the last one is detected as an outlier, by taking it out no new outliers are detected. Last step.
        if(length(obvious_outliers) == 0){ # the last one
          input1 <- input[!input$ids %in% c(obvious_outliers_final,i) & input$ids >= first_point & input$ids <= last_point + 1,] # redefine the df without the outlier
          
          slopes <- diff(input1$weight) / diff(input1$time)
          ids <- input1$ids[-1]
          units <- input1$units[-1]
          
          transpi <- slopes / units
          transpi_times <- input1$time[-nrow(input1)] + (diff(input1$time)/2)
          transpi_lightPeriod <- dd_light(dd= transpi_times)
          cdatab <- data.frame(ids,slopes,transpi,transpi_times,transpi_lightPeriod)
          cdatab <- cdatab[cdatab$ids >= first_point & cdatab$ids <= last_point,]
          
          dark_ids <- cdatab$ids[!cdatab$ids %in% transition_ids & cdatab$transpi_lightPeriod == "dark" & (cdatab$transpi < -max_dark | cdatab$transpi > -min_dark )]
          light_ids <- cdatab$ids[!cdatab$ids %in% transition_ids & cdatab$transpi_lightPeriod == "light" & (cdatab$transpi < -max_light | cdatab$transpi > -min_light)]
          obvious_outliers_final <- c(obvious_outliers_final,dark_ids,light_ids,transition_ids_outliers)
        }
        
      }
    }
  }
  
  
  # outliers on sides
  out_sides <- input$ids[input$ids < first_point | input$ids > last_point + 1]
  
  # specific outliers
  ids_final_out <- c(ids_final_out,specific_outliers)
  
  # obvious outliers
  # take out the specific ids outliers which are in the obvious to avoid hotspot calculation afterwards
  if(length(obvious_outliers_final) > 0){ ids_final_out <- ids_final_out[!ids_final_out %in% obvious_outliers_final]}
  
  # output with ID's
  return(list(hotspots = ids_final_out,out_sides = out_sides,out_obvious = obvious_outliers_final))
  
}

Outliers_v5 <- function(time,weight,lightPeriod,outlier, hotspot,startdark,period,min_around,daynight){
  #----------------------------------------------------------------------#
  #----------------------explanation-------------------------------------#
  #----------------------------------------------------------------------#
  
  # This is an R script to analyse the influence of points of previous selected hotspots of ONE pot.
  # The script will return a vector with ID's for the type of outlier, including information about Cooks, residuals and dfbetas.
  
  # inputs:
  #     - time: decimalDay vector 
  #     - weight: weight vector in g
  #     - startdark: start of the darkperiod in decimaldays
  #     - period: hours of dark period in decimaldays
  #     - min_around: minutes around the hotspot to select the interval  
  #     - daynight: taking into account the day and night when selecting the interval for points. 
  
  #----------------------------------------------------------------------#
  #----------------------environment-------------------------------------#
  #----------------------------------------------------------------------#
  
  # Functions
  dd_light <- function(dd){
    # work only with hours!! not with day data
    
    lightp_final <- c()
    
    for (i in dd){
      # i = dd[1]
      hourday <- as.numeric(sub(".*\\.","0.",sprintf("%.10f", round(i,10)) )) # decimal time on the day
      
      lightp <- rep("light",length(hourday))
      
      if (startdark + period > 1){
        # night during night
        
        lightp[hourday > startdark | hourday <= enddark] <- "dark"
        lightp_final <- c(lightp_final,lightp)
      }
      
      if (startdark + period < 1){
        # night during day
        
        lightp[hourday > startdark & hourday <= enddark] <- "dark"
        lightp_final <- c(lightp_final,lightp)
      }
      
      
    }
    
    return(lightp_final)
    
  }
  #----------------------------------------------------------------------#
  #----------------------start-------------------------------------------#
  #----------------------------------------------------------------------#
  
  # Put the time vector in the good order
  input <- data.frame(time,weight,lightPeriod)
  input <- input[order(input$time),]
  
  # add an ID colum 
  input$ids <- seq(1,nrow(input),by = 1)
  
  input <- input[!outlier,]
  
  selecpoint <- function(timer,daynight){
    
    middledater <- input$time[input$ids == hotspot] # get the time
    
    # select values based on certain time frame
    beforedater <- middledater - (timer / (24*60)) #  time in decimals hour
    afterdater <- middledater + (timer / (24*60)) # add 1 hour
    # select index values 1 hour before and 1 hour after
    dataindex <- input$ids[input$time >= beforedater & input$time <= afterdater]
    
    # Which one is the most close ?
    #closest <- which.min(as.numeric(Pesees$decimalDay - middledater)^2)
    
    
    if (length(dataindex) > 0){
      # middledark = Pesees$lightPeriod[closest] # get the lightperiod 
      middledark = dd_light(dd = middledater)
      
      inputdark = lightPeriod[dataindex]
      
      
      if (daynight == "yes"){
        
        # if(any(inputdark == "darklight")){
        #   
        #   if(match("darklight",inputdark) < match(middledark,inputdark)){
        #     # means that we are just after a transition and that we should not include the darklight
        #     dataindex1 <- dataindex[inputdark %in% c(middledark)] # darklight is special for points close to transition
        #   }else{
        #     # we have a transition before us...
        #     dataindex1 <- dataindex[inputdark %in% c(middledark,"darklight")] # darklight is special for points close to transition
        #     
        #   }
        # }else{
        # not close to transitions
        dataindex1 <- dataindex[inputdark %in% c(middledark,"darklight")]
        # }
        
        
      }
      if (daynight == "no"){
        dataindex1 <- dataindex
        
      }
      inputtime <- diff(range(input$time[input$ids %in% dataindex1]))
      
      
    }else{
      # if not enough points to run the model
      dataindex1 <- integer()
      inputtime <- integer()
      inputdark <- integer()
      middledark <- integer()
      
    }
    
    return(list(dataindex1,inputdark,middledark,inputtime))
    
  }
  
  timer = min_around  #  time around measurement taken into account (min) which is dependent on the frequency of the measurements (so to have at least 4, times 2)
  
  # selects points in the range based on interval and light/dark period
  output <- selecpoint(timer= timer,daynight = daynight)
  dataindex1 <- unlist(output[1])
  
  # when less or equal to 2 try to include points which are close to the border
  while (length(dataindex1) < 5 & timer < 180){
    timer = timer + 20
    output <- selecpoint(timer = timer,daynight = daynight)
    dataindex1 <- unlist(output[1])
  }
  
  if (length(dataindex1) > 4){
    
    inputweights = input$weight[input$ids %in% dataindex1]
    # decimal day
    inputdates = input$time[input$ids %in% dataindex1]
    #plot(inputweights~inputdates)
    
    # robust regression (weighted)
    cdata <- data.frame(dataindex1,inputdates,inputweights)
    ols <- lm(inputweights ~ inputdates, data = cdata)
    
    abc = influence.measures(ols, infl = influence(ols))
    cdata <- as.data.frame(abc$infmat)
    
    star <- rstandard(ols)# extract standardized residuals from a linear model
    stur <- rstudent(ols)# extract studentized residuals from a linear model
    
    cdata <- cbind(cdata,star,stur)
    
    cdata$outlier <- ifelse(abs(cdata$dffit) > (8*sqrt(2+1))/(nrow(cdata)-2-1),TRUE,FALSE) # limit dffit
    
    outputdark <- unlist(output[2])
    middledark <- unlist(output[3])
    
    # determine transition points: begin
    if(cdata$outlier[1] == T | is.na(cdata$outlier[1])){
      # so one at the side which is outlier
      
      
      if(min(dataindex1) %in% c(1,2)){
        cdata$outlier[1] <- ifelse(abs(cdata$dffit[1]) > (10*sqrt(2+1))/(nrow(cdata)-2-1),TRUE,FALSE)
        
      }else{
        if(!lightPeriod[dataindex1[1]-2] %in% c(middledark)){
          
          cdata$outlier[1] <- ifelse(abs(cdata$dffit[1]) > (50*sqrt(2+1))/(nrow(cdata)-2-1),TRUE,FALSE)
        }else{
          
          cdata$outlier[1] <- ifelse(abs(cdata$dffit[1]) > (15*sqrt(2+1))/(nrow(cdata)-2-1),TRUE,FALSE)
        }
      }
    }
    
    # end of period
    if(cdata$outlier[nrow(cdata)] == T | is.na(cdata$outlier[nrow(cdata)])){
      # so one at the side which is outlier
      
      if(!lightPeriod[max(dataindex1)+2] %in% c(middledark)){
        
        cdata$outlier[nrow(cdata)] <- ifelse(abs(cdata$dffit[nrow(cdata)]) > (50*sqrt(2+1))/(nrow(cdata)-2-1),TRUE,FALSE)
      }else{
        cdata$outlier[nrow(cdata)] <- ifelse(abs(cdata$dffit[nrow(cdata)]) > (15*sqrt(2+1))/(nrow(cdata)-2-1),TRUE,FALSE)
      }
      
      
    }
    
    
    outliers <- cdata$outlier
    
    
    # 2/sqrt(nrow(cdata)) # limit dfbetas
    
  }else{
    warning(paste("not enough hotspot points in ",i,hotspot,", only point transformed in outlier"))
    outliers <- rep(TRUE,length(dataindex1))
  }
  
  
  return(dataindex1[outliers])
  
  
}




#------------------------------------------------------------------------------#
#        Format the dataframe with a uniform time sequence for all pots        #
#------------------------------------------------------------------------------#

# This function written by FP was previously in "PhenoLeaks_outliers.R",
# but has moved here because it is now called at Step#01.
format_time <- function (all_dat,
                         Time_var = if ("decimalDay" %in% names(all_dat)) "decimalDay" else "Time", # the column name of the time axis
                         Trt_var, # the column name(s) of the statistical treatment(s), e.g. "idGenotype" or c("idGenotype", "idWatering")
                         time_step = 30) # time step between two observations (min) - strongly advised to be an integer >= 1
  {
  # Force showing warnings as they occur
  warn_opt <- options("warn")$warn
  options(warn = 1)
  
  # Set the number of digits for rounding time
  n.digits <- 7
  
  # Check that there is the same number of observations (even if missing) for all pots and at the right times.
  # This is normally done at Step#01. If not, the dataframe is extended and a warning message is shown for each line added.
  mini <- min(all_dat[, Time_var])
  maxi <- max(all_dat[, Time_var])
  time_vector <- seq(mini, maxi, by = time_step/60/24)
  nobs_per_pot <- aggregate(list(N = all_dat$idPot), by = list(idPot = all_dat$idPot), FUN = length)
  
  if (max(nobs_per_pot$N) != length(time_vector))
    {
    if (max(nobs_per_pot$N) > length(time_vector))
      {
      too_many_data <- nobs_per_pot$idPot[which(nobs_per_pot$N > length(time_vector))]
      stop (paste("More observations than expected based on the time step.\nCheck pots ", paste(too_many_data, collapse = ", "), ".", sep = ""))
      }
    else if (max(nobs_per_pot$N) < length(time_vector))
      {
      nobs_missing <- length(time_vector) - max(nobs_per_pot$N)
      warning (paste(nobs_missing, if (nobs_missing == 1) "observation is" else "observations are", "missing in all pots.", sep = " "))
      }
    }
  
  if (max(nobs_per_pot$N) < length(time_vector) | length(unique(nobs_per_pot$N)) != 1)
    {
    incomplete_pots <- nobs_per_pot$idPot[nobs_per_pot$N < length(time_vector)]
    for (pot in incomplete_pots)
      {
      time_pot <- round(all_dat[all_dat$idPot == pot, Time_var], n.digits)
      missing_times <- which(!round(time_vector, n.digits) %in% time_pot)
      n_miss <- length(missing_times)
      trt <- if (length(Trt_var) == 1) all_dat[all_dat$idPot == pot, Trt_var][1] else all_dat[all_dat$idPot == pot, Trt_var][1,]
      
      all_dat[nrow(all_dat)+1:n_miss, Time_var] <- time_vector[missing_times]
      all_dat[(nrow(all_dat)-n_miss+1):nrow(all_dat), "idPot"] <- pot
      all_dat[(nrow(all_dat)-n_miss+1):nrow(all_dat), Trt_var] <- trt
      warning (paste(n_miss, if (n_miss == 1) "line" else "lines", "added for pot", pot, sep = " "))
      }
    }
  
  # Order by pots and time, Rename lines and keep order in the 'idObs' column
  all_dat <- all_dat[order(all_dat$idPot, all_dat[, Time_var]), ]
  row.names(all_dat) <- 1:nrow(all_dat)
  all_dat <- cbind(data.frame(idObs = 1:nrow(all_dat)), all_dat)
  
  # # Round the time for comparison between pots (otherwise numeric errors may occur) 
  # all_dat[, Time_var] <- round(all_dat[, Time_var], n.digits)
  
  # Shift time so that it is 0 at Time_ON0
  all_dat$Time <- all_dat[, Time_var] - Time_ON0
  
  # Create a time column in minutes that corresponds exactly to the initial prediction time of transpiration
  # (avoid numeric issues that frequently arise when manipulating decimal numbers with infinite digits)
  all_dat$Exact_Time_min <- as.integer(round(all_dat$Time * 60 * 24))
  
  # Back to user warning options
  options(warn = warn_opt)
  
  # Return
  return(all_dat)
  }



