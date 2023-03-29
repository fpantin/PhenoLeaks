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

#------------------------------------------------------------------------------#
#             Function to calculate decimalDay from Posix data                 #
#------------------------------------------------------------------------------#
# can calculate the decimal day of date data but also datehour data

decimalDay <- function(column){
  input <- as.POSIXct(strptime(column, format= "%d/%m/%Y",tz = "UTC"))
  if(all(is.na(input))){
    input = as.POSIXct(strptime(column,format= "%Y-%m-%d %H:%M:%S", tz = "UTC")) # now the date/hour is not changing
  }
  output <- as.numeric(strftime(input, format = "%j", 
                                          tz='UTC')) + (as.numeric(strftime(input, format = "%H", tz='UTC'))/(24)) +
    (as.numeric(strftime(input, format = "%M", tz='UTC'))/(24*60)) +
    (as.numeric(strftime(input, format = "%S", tz='UTC'))/(24*60*60)) 
  
  return(output)
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

# This function need a dataframe including "idPotManip", "outlier" (identified before manually), 
# "Areamm2", "decimalDay".

surf_fit <- function(surface){
  
  surf_coef <- data.frame(idPotManip = NULL, model = NULL,growth = NULL, from = NULL,
                          to = NULL, intercept = NULL, slope = NULL)
  
  mod = "test"
  growth = "continu"
  
  for (i in unique(surface$idPotManip)){
    # i = "C2M43-93"
    
    model = mod
    input = surface[surface$idPotManip == i & !surface$outlier,]
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
      
      surf_coef <- rbind(surf_coef,data.frame(idPotManip = i, model = model,growth = growth, from = min(input$decimalDay), to = max(input$decimalDay), intercept = as.numeric(lm$coefficients[1]), slope = as.numeric(lm$coefficients[2])))
      
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
  for (i in unique(input$idPotManip)){
    
    if(i %in% unique(surf_coef$idPotManip)){
      i = as.character(i)
      
      mod <- surf_coef$model[surf_coef$idPotManip == i & surf_coef$growth == "continu"]
      
      slope = surf_coef$slope[surf_coef$idPotManip == i & surf_coef$growth == "continu"]
      
      intercept = surf_coef$intercept[surf_coef$idPotManip == i & surf_coef$growth == "continu"]
      
      if (mod == "lm"){
        
        input$surface[input$idPotManip == i] <- intercept +  input$decimalDay[input$idPotManip == i] * slope
      }
      if (mod == "log"){
        input$surface[input$idPotManip == i] <- exp(intercept +  input$decimalDay[input$idPotManip == i] * slope)
      }
    }else{
      input$surface[input$idPotManip == i] <- NA
      no_surf <- c(no_surf,i)
    }
    
  }
  print(unique(no_surf))
  print(length(unique(no_surf)))
  
  return(input)
}

