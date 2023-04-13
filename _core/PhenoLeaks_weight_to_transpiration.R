################################################################################
#                                                                              #
#                             PhenoLeaks - Transpiration                       #
#                                                                              #
#                 Function to estimate transpiration by statistical modelling  #
#                             Adriaan Westgeest, 2023                          #
#                                                                              #
################################################################################




#------------------------------------------------------------------------------#
#                                                       #
#------------------------------------------------------------------------------#
Transpi_calc_v7 <-  function(input, freq = 30, min_around = 120,nightperiod = 24,lightsOFF = "00:00:00",method = "lm", 
                             smooth.fac = 0.5, nop = 2,daynight = "yes",
                             max = 180,max_end = 120,gamma = 1, m  = "default",k=NA,variable = "E"){
  # This is the same transpiration script as v1 but is not corrected for the surface!
  # max = how much the total span can be maximum
  # v5 has an adaptive interval which increases till the max gradually. 
  
  # test:
  # input = df_outlier_2[!df_outlier_2$outlier,];lightsOFF = "18:00:00"; 
  # nightperiod = nightperiod; method = "lm";daynight = "yes";freq = 30;min_around = 90;
  # max = 180;max_end = 120;nop=2
  
  # input file:
  # 
  
  # input = df
  # freq = 30 (in minutes)
  # min_around = 90 (interval on the left/right, not the total interval)
  # lightsOFF = "00:00:00" (standard in the case there is no difference between day and night)
  # lightperiod = 24 (standard in the case there is no difference between day and night)
  # method = "lm" 
  # max = 180
  # nop = 4
  # daynight = "no"
  #--------------------------------------------#
  # eval(parse(text=paste0("input$Surface <- input$",surface)))
  # eval(parse(text=paste0("input$Weight_corr <- input$",weight)))
  
  #dayvec <- floor(dd) + OFF
  OFF <- as.numeric(strftime(as.POSIXct(strptime(lightsOFF,format= "%H:%M:%S"), tz = "UTC"),format="%H",tz = "UTC"))/(24) +
    as.numeric(strftime(as.POSIXct(lightsOFF,format= "%H:%M:%S", tz = "UTC"), format = "%M", tz='UTC'))/(24*60)
  period <- nightperiod / 24
  startdark <- as.numeric(round(OFF,10))
  enddark <- startdark + period
  enddark <- as.numeric(sub(".*\\.","0.",sprintf("%.10f", round(enddark,10)) ))
  
  # functions
  dd_light <- function(dd){
    # work only with hours!! not with day data
    
    hourday <- as.numeric(sub(".*\\.","0.",sprintf("%.10f", round(dd,10)) )) # decimal time on the day
    
    if (startdark + period > 1){
      # night during night
      if (hourday > startdark | hourday <= enddark){
        lightp <- "dark"
      }else{
        lightp <- "light"
      }
    }
    if (startdark + period < 1){
      # night during day
      if (hourday > startdark & hourday <= enddark){
        lightp <- "dark"
      }else{
        lightp <- "light"
      }
    }
    
    # but test here whether in VPD period on decimalDay
    if("VPD"%in% unique(input$lightPeriod)){
      # means that there are VPD periods
      
      # test for every vpd period
      for (i in seq(1,length(vpd_periods)-1,2)){
        if (dd > vpd_periods[i] & dd< vpd_periods[i+1]){
          # test if within a vpd period
          lightp <- "VPD"
        }
      }
      
      
    }
    
    
    return(lightp)
    
  }
  #--------------------------------------------#
  
  M <- 18.015 # molecular weight of water (g/mol)
  
  
  # new data frame
  transpi <- data.frame(idPot = NULL, 
                        decimalDay = NULL,
                        HourOnDay = NULL,
                        E = NULL,
                        nobs = NULL,
                        min_decimalDay= NULL,
                        max_decimalDay = NULL)
  
  #idpotlist <- c("C3M34B-636","C3M34B-500","C3M34B-499","C3M34B-570","C3M34B-637")
  idpotlist <- unique(input$idPot)
  
  selecpoint <- function(timer,daynight){
    
    middledater <- spectimevec[t] # get the time
    
    # select values based on certain time frame
    beforedater <- middledater - (timer / (24*60)) #  time in decimals hour
    afterdater <- middledater + (timer / (24*60)) # add 1 hour
    # select index values 1 hour before and 1 hour after
    dataindex <- which(Pesees$decimalDay >= beforedater & Pesees$decimalDay <= afterdater)
    
    # Which one is the most close ?
    #closest <- which.min(as.numeric(Pesees$decimalDay - middledater)^2)
    
    
    if (length(dataindex) > 0){
      # middledark = Pesees$lightPeriod[closest] # get the lightperiod 
      middledark = dd_light(dd = middledater)
      
      inputdark <- Pesees$lightPeriod[dataindex] 
      
      if (daynight == "yes"){
        
        dataindex1 <- dataindex[inputdark %in% c(middledark,"darklight")] # darklight is special for points close to transition
      }
      
      if (daynight == "no"){
        dataindex1 <- dataindex
      }
      
      if (daynight == "period"){
        
        if(middledark == "light"){
          # test in which period we are: light is easy
          beforedater <- as.integer(middledater) + enddark
          afterdater <- as.integer(middledater) + startdark
          
        }
        if (middledark == "dark"){
          # can be 2 periods, which one?
          if(middledater <= as.integer(middledater) + enddark){
            # in the morning
            beforedater <- as.integer(middledater) + enddark - period
            afterdater <- as.integer(middledater) + enddark
            
          }else{
            # in the evening
            beforedater <- as.integer(middledater) + startdark
            afterdater <- as.integer(middledater) + startdark + period
            
          }
        }
        dataindex1 <- which(Pesees$decimalDay >= beforedater & Pesees$decimalDay <= afterdater) #including the first point (imputed)
        inputdark <- Pesees$lightPeriod[dataindex1] 
        #stopifnot(all(inputdark %in% c(middledark,"darklight")))
      }
      
      if(length(dataindex1) > 0) {
        inputtime <- diff(range(Pesees$decimalDay[dataindex1]))
        
      }else{
        inputtime <- integer()
      }
      
    }else{
      # if not enough points to run the model
      dataindex1 <- integer()
      inputtime <- integer()
      inputdark <- integer()
      middledark <- integer()
      
    }
    
    return(list(dataindex1,inputdark,middledark,inputtime))
    
  }
  
  suppressPackageStartupMessages(library("svMisc")) #for progress function
  
  #suppressPackageStartupMessages(library("MASS"))
  
  ###########################################################################
  # modification: get the same number of transpiration values for every pot #
  ###########################################################################
  if(variable == "E"){
    
    a = seq(from= 0, to = 1440, by= freq) # vector by minutes for a day
    
    timevec <- a / (24* 60) # 24 hours
    timevec <- timevec[-length(timevec)] # get rid of the last one
    
    x = c(min(input$decimalDay,na.rm = T), max(input$decimalDay,na.rm = T))
    a = diff(x)
    
    # make a specific count up timevector for the whole range (decimalDay?)
    # select dec days
    
    unidays = unique(as.integer(input$decimalDay))
    spectimevec <- vector()
    
    # make a new vector containing every dec day in 30 min
    for (dec in 1:length(unidays)){
      spectimevec <- append(spectimevec, timevec + unidays[dec])
    }
    
    # select the index values which are in the range of i
    # at the start we do not want extrapolation, however at the end we would like to (last point of night)
    # select the last point of the night 
    interpol <- 4 / 24  #interpolation at the end of 2 hour
    trans_max <- as.integer(max(input$decimalDay)) + c(startdark,enddark)
    
    val_max <- max(input$decimalDay)
    if(max(input$decimalDay) <= trans_max[1] & max(input$decimalDay) >= trans_max[1] - interpol){
      # close to transition of day to night
      val_max = trans_max[1]
    } 
    if(max(input$decimalDay) <= trans_max[2] & max(input$decimalDay) >= trans_max[2] - interpol){
      # close to transition of day to night
      val_max = trans_max[2]
    }
    
    ndxvalues <- which(spectimevec >= (min(input$decimalDay)) & spectimevec <= val_max)
    
  }
  
  for (i in idpotlist){
    # i = idpotlist[1]
    ######## OUTLIERR
    #Pesees <- Workdata[which(Workdata$idPotManip == i & 
    #                         Workdata$Col == "blue"),]
    
    # numb <- match(i,unique(input$idPotManip))
    # progress(numb,max.value = max(order(idpotlist)))
    # if (numb == length(idpotlist) ) cat("Done!\n")
    # 
    # i = "C2M47-233"
    Pesees <- input[which(input$idPot == i),]
    Pesees <- Pesees[order(Pesees$decimalDay),]
    
    
    options(warn=2)
    
    for (t in ndxvalues){
      # piece-wise
      
      # select for every half an hour on the day the values 
      # t = 36
      #-----------------------------------------------#
      #             selection of measurements         #
      #-----------------------------------------------#
      # filter for time period
      #minlm <- min_around / nop
      #minlm <- min_around / nop  # minimum number of points as input for lm = 3 = preferred because have the end of night after for pgm, if not not!
      # disavantage =  more variation end of day/night.
      
      minlm = nop
      timer = min_around
      
      output <- selecpoint(timer = timer,daynight = daynight)
      dataindex1 <- unlist(output[1])
      #inputtime <- unlist(output[4])
      inputdark <- unlist(output[2])
      middledark <- unlist(output[3])
      
      # while (all(inputdark == middledark) & timer <= max){
      #   # increase the interval as long as all the points are in the same time period.
      #   # does not add something when darklight is added as increasing means inavoidable the selection
      #   # of points in the other period! 
      #   # avoid inequal intervals at both sides of the point!
      #   timer <- timer + 15
      #   output <- selecpoint(timer = timer,minlm = minlm,daynight = daynight)
      #   dataindex1 <- unlist(output[1])
      #   #inputtime <- unlist(output[4])
      #   inputdark <- unlist(output[2])
      #   middledark <- unlist(output[3])
      #   
      # }
      
      if (length(dataindex1) > 0 & daynight %in% c("yes")){
        # only if points are found and if daynight is yes
        time1 <- timer
        inputdarkness <- inputdark
        
        if(all(inputdarkness == middledark)){
          # condition when all input is dark:
          while ( all(inputdark == middledark) & time1 <= max){
            # increase the interval as long as all the points are in the same time period.
            # does not add something when darklight is added as increasing means inavoidable the selection
            # of points in the other period! 
            time1 <- time1 + 15
            output <- selecpoint(timer = time1,daynight = daynight)
            dataindex1 <- unlist(output[1])
            #inputtime <- unlist(output[4])
            inputdark <- unlist(output[2])
            middledark <- unlist(output[3])
            
          }
        }
        if(!inputdarkness[length(inputdarkness)] == middledark){ 
          # condition when at the end of a lightperiod:
          while (!inputdark[length(inputdark)] == middledark & timer <= max_end){
            # increase the interval as long as all the points are in the same time period.
            # does not add something when darklight is added as increasing means inavoidable the selection
            # of points in the other period! 
            timer <- timer + 15
            output <- selecpoint(timer = timer,daynight = daynight)
            dataindex1 <- unlist(output[1])
            #inputtime <- unlist(output[4])
            inputdark <- unlist(output[2])
            middledark <- unlist(output[3])
            
          }
        }
        
        
      }
      
      
      
      
      if (length(dataindex1) >= minlm){
        
        inputweights = Pesees$weight[dataindex1]
        inputdates = Pesees$decimalDay[dataindex1]
        cdata <- data.frame(inputdates,inputweights)
        
        if (method == "lm"){
          mod <- lm(inputweights ~ inputdates, data = cdata) 
          beta <- -1 * (as.numeric(mod$coefficients[2]))
          
        }
        if (variable == "E"){
          Ef <-  beta / (M * (10^-3) * 24* 60 * 60)
        }
        
        if (variable == "weight"){
          Ef <- predict(mod, data.frame(inputdates = spectimevec[t]))
          
        }
        min_decimalDay = min(cdata$inputdates)
        max_decimalDay = max(cdata$inputdates)
        
      }else{
        Ef <-  NA
        min_decimalDay <- NA
        max_decimalDay <- NA
      }
      
      tt <- spectimevec[t]
      hourday <- as.numeric(sub(".*\\.","0.",sprintf("%.10f", round(tt,10)) )) # decimal time on the day
      # hourconsec <- timevec[t] # consecutive time. 
      # idtreat <- Pesees$Treatment[1]
      # idgeno <- Pesees$idGenotype[1]
      transpi <- rbind(transpi,data.frame(idPot = i, decimalDay = tt, 
                                          HourOnDay = hourday, E = Ef,
                                          nobs = length(dataindex1), min_decimalDay = min_decimalDay, max_decimalDay = max_decimalDay))
      
      
    }
    
    
  }
  return(transpi) 
  
  
}