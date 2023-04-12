################################################################################
#                                                                              #
#         PhenoLeaks - Step #00 - Define the features of the experiment        #
#                                                                              #
#                   Script to set constants of the experiment                  #
#                           and graphical parameters                           #
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
#                                    C2M43B                                    #
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
for (pkg in c("igraph", "here"))
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


# Note that paths are built relative to the root directory of the PhenoLeaks project.
dir_PhenoLeaks <- here::here("_core")
source(file.path(dir_PhenoLeaks, "PhenoLeaks_generic.R"))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                              (3)  Set constants                              #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#------------------------------------------------------------------------------#
#                  Define and lock constants of the experiment                 #
#------------------------------------------------------------------------------#

# Set the constants within a user-defined function

set_constants_C2M43B <- function ()
  {
  # Plant species
  spcs <- "Arabidopsis"
  
  # Name of the experiment
  idExp <- "C2M43B"

  # Duration of the photoperiod (h)
  Pho_Per <- 12
  
  # Duration of the skotoperiod (h)
  # Note that the scripts will not work for diel cycles that do not match 24 h.
  Sko_Per <- 24 - Pho_Per
  
  # Chosen reference for initializing time at a night-to-day transition ("decimalDay" -> "Time" from 0) (d)
  Time_ON0 <- 103+8/24
  
  # Time when the experiment was set up and the balance started to be tested (d)
  #Time_setup <- -1 + 1.5/24 # useless here since Time_setup = Time_start_exp
  
  # Start of the experiment proper (d)
  Time_start_exp <- -1 + 1.5/24 - 1/60/24 # 1 min subtracted to avoid rounding issues (the 'real' value should be included in the interval Time_start_exp:Time_end_exp)
  
  # End of the proper experiment (d)
  Time_end_exp <- 5 + 2/24 + 1/60/24 # 1 min added to avoid rounding issues (the 'real' value should be included in the interval Time_start_exp:Time_end_exp)

  
  
  # Duration covering the rapid stomatal movements at the transitions (min)
  # Corresponds here to the first data point after a day/night transition.
  Rapid_dur <- 30
  
  # Duration by which the square wave (see PhenoLeaks_fit) is shifted compared to the dark/light cycle (min)
  # Should be close to, but lower than, 'Rapid_dur'.
  # Should match a time point where there is no data.
  # Corresponds here to the middle of the interval between the point at the transition and the one immediately after.
  SQW_shift <- 15 # (epsilon in Westgeest et al. 2022)
  
  

  # Lock the constants in the global environment and return them as a list
  #
  #       ---> --> ->   DO NOT CHANGE THIS LAST LINE OF CODE!   <- <-- <---
  #
  lock_constants(environment())
  } 

# Call the function
set_constants_C2M43B()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                         (4)  Set graphical parameters                        #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#------------------------------------------------------------------------------#
#            Set color and line type for treatments and environments           #
#------------------------------------------------------------------------------#

# Set the colors and line types within a user-defined function

set_colors_C2M43B <- function ()
  {
  ColorsTrt <- c(rgb(0, 0, 0, maxColorValue=255),  # Col-0 = black
                 rgb(237, 125, 49,maxColorValue=255), # pgm = orange
                 rgb(90, 181, 75, maxColorValue=255), # bam1 = dark green
                 rgb(65, 91, 39, maxColorValue=255), # bam1 bam3 = dark dark green
                 rgb(179, 218, 38, maxColorValue=255), # bam3 = lighter green yellow
                 rgb(21, 92, 171, maxColorValue=255)) # sex1 = bleuish
  
  ColorsTrt <- data.frame(idGenotype = c("Col-0", "pgm", "bam1", "bam1 bam3", "bam3", "sex1-3"),
                          idWatering = "WW",
                          col = ColorsTrt,
                          lty = 1)
  
  ColorsTrt$Trt <- paste(ColorsTrt$idGenotype, ColorsTrt$idWatering, sep = " - ")

  
  t1 <- c(-1, 3, 4, 5)
  t2 <- c(3, 4, 5, 6)
  ColorsPeriod <- data.frame(idPeriod = c("Control", "Low light 1", "Low light 2", "Low light 3"),
                             Time1 = t1,
                             Time2 = t2,
                             decimalDay1 = Time_ON0 + t1,
                             decimalDay2 = Time_ON0 + t2,
                             col = c("black", "burlywood4", "burlywood4", "burlywood4"))
  
  
  return (list(ColorsTrt = ColorsTrt, ColorsPeriod = ColorsPeriod))
  }


# Call the function
c(ColorsTrt, ColorsPeriod) := set_colors_C2M43B()



#------------------------------------------------------------------------------#
#                             Check irrigation data                            #
#------------------------------------------------------------------------------#

irrig_file <- file.path(here::here(), spcs, idExp, "Processed_data", "C2M43B_info_irrigation.csv")
if (file.exists(irrig_file))
  {
  # Import and manage irrigation data
  irrig <- read.csv(irrig_file)
  irrig <- irrig[order(irrig$idPot, irrig$decimalDay), ]
  
  # Check the number of irrigation events for each pot
  aggregate(list(N = irrig$decimalDay), by = list(idPot = irrig$idPot), FUN = length)
  # --> All pots have 1 irrigation event except pot 108 that has 2 (one additional event detected early in the kinetics)
  # Inspection of the transpiration data shows this pot will be removed (erratic data), so delete this early event
  irrig <- irrig[-which(irrig$idPot == "108" & irrig$decimalDay < 103), ]
  aggregate(list(N = irrig$decimalDay), by = list(idPot = irrig$idPot), FUN = length)
  # --> OK
  
  # Get the mean time of each irrigation event
  irrig$event <- NA
  for (pot in unique(irrig$idPot))
    {
    irrig$event[irrig$idPot == pot] <- 1:(nrow(irrig[irrig$idPot == pot,]))
    }
  mean_irrig_time <- aggregate(list(MeanTime = irrig$decimalDay), by = list(event = irrig$event), FUN = mean)
  }


#------------------------------------------------------------------------------#
#                        Depict environmental conditions                       #
#------------------------------------------------------------------------------#


# Define and draw upper rectangles representing environmental conditions
# It does not have to be run here, but it will be called by 'prepare_kin()' within 'PhenoLeaks_graphics.R'.
# It should be named as "rectangles_environment_" + the identifier of the experiment ('idExp')

rectangles_environment_C2M43B <- function (col.per = ColorsPeriod, # the dataframe containing the features of the different periods
                                          Time_var, # "Time" or "decimalDay"
                                          irrig_show_mode = "mean", # show an arrow at irrigation; use "mean" for the mean time(s) of irrigation across the experiment, "pot" for the irrigation specific to a given pot, or "none" for no arrow
                                          pot, # the pot id for showing irrigation (only if irrig_show_mode = "pot")
                                          cex.env = 0.55, # font size for labeling the environments
                                          inside = T, # should the rectangles be drawn inside or outside the box?
                                          y_prop = if (export_PPTX) 0.08 else 0.05, # height of the rectangles (in proportion to the total range of the y-axis)
                                          export_PPTX = F) # set to TRUE to export the graph as PPTX instead of PDF
  {
  usr <- par("usr")
  y1 <- usr[4]
  if (inside) y_prop <- -y_prop
  y2 <- y1 + y_prop*(y1 - usr[3])
  if (!export_PPTX) clip(usr[1], usr[2], y1, y2)
  rect(col.per[, paste(Time_var, 1, sep = "")],
       y1,
       col.per[, paste(Time_var, 2, sep = "")],
       y2,
       col = col.per$col, border = NA, density = NULL)
  if (!export_PPTX) clip(usr[1], usr[2], usr[3], usr[4])
  
  y.text <- (y1 + y2) / 2 - if (export_PPTX) 0.015 else 0
  t1.low.light <- col.per[2, paste(Time_var, 1, sep = "")]
  text(t1.low.light+Pho_Per/24/2+0:2, y.text, "Low light", cex = cex.env, xpd = T)
  segments(c(t1.low.light)+Pho_Per/24+0:1, y1, c(t1.low.light)+Pho_Per/24+0:1, y2, col = "gray", lty = "22", xpd = T)
  segments(c(t1.low.light)+1:2, y1, c(t1.low.light)+1:2, y2, col = "gray", lty = "22", xpd = T)
  
  if (irrig_show_mode != "none" & file.exists(irrig_file))
    {
    if (irrig_show_mode == "mean") t_irrig <- mean_irrig_time$MeanTime
    else if (irrig_show_mode == "pot") t_irrig <- irrig$decimalDay[irrig$idPot == pot]
    else stop ("'irrig_show_mode' should be \"none\", \"mean\" or \"pot\".")
    if (Time_var != "decimalDay") t_irrig <- t_irrig - Time_ON0
    text(t_irrig, usr[3] + 0.025*(usr[4] - usr[3]), "Rehydration", col = "black", cex = cex.env)
    #require(igraph)
    iArrows <- igraph:::igraph.Arrows
    iArrows(t_irrig, rep(usr[3] + 0.05*(usr[4] - usr[3]), length(t_irrig)),
            t_irrig, rep(usr[3] + 0.15*(usr[4] - usr[3]), length(t_irrig)),
            h.lwd = 1, sh.lwd = 1, sh.col = "black", curve = 0, width = 1, size = 0.5)
    }
  }

