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
#                                     C2M47                                    #
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

set_constants_C2M47 <- function ()
  {
  # Plant species
  spcs <- "Arabidopsis"
  
  # Name of the experiment
  idExp <- "C2M47"

  # Duration of the photoperiod (h)
  Pho_Per <- 12
  
  # Start of the gravimetric experiment (YYYY-MM-DD HH:MM:SS)
  from <- "2019-02-28 16:00:00"
  
  # End of the gravimetric experiment (YYYY-MM-DD HH:MM:SS)
  to <- "2019-03-08 23:00:00"
  
  # Duration of the skotoperiod (h)
  # Note that the scripts will not work for diel cycles that do not match 24 h.
  Sko_Per <- 24 - Pho_Per
  
  # Chosen reference for initializing time at a night-to-day transition ("decimalDay" -> "Time" from 0) (d)
  Time_ON0 <- 59.75 + 2
  
  # Time when the experiment was set up and the balance started to be tested (d)
  Time_setup <- -2
  
  # Start of the experiment proper (d)
  Time_start_exp <- Time_setup + 1 + Pho_Per/2/24
  
  # End of the proper experiment (d)
  Time_end_exp <- Time_setup + 8
  
  # time in UTC when the light was switched off and the darkperiod starts (HH:MM:SS)
  lightsOFF <-  "06:00:00" 
  
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
set_constants_C2M47()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                         (4)  Set graphical parameters                        #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#------------------------------------------------------------------------------#
#            Set color and line type for treatments and environments           #
#------------------------------------------------------------------------------#

# Set the colors and line types within a user-defined function

set_colors_C2M47 <- function ()
  {
  ColorsTrt <- c(rgb(0, 0, 0, maxColorValue=255),  # Col-0 = black
                 rgb(237, 125, 49,maxColorValue=255), # pgm = orange
                 rgb(56, 200, 176, maxColorValue=255), # abcb14-1 = pool
                 rgb(88, 168, 153, maxColorValue=255), # abcb14-2 = dark pool
                 rgb(25, 248, 8, maxColorValue=255), # amy3 = fel green
                 rgb(81, 211, 45, maxColorValue=255), # amy3 bam1 = green
                 rgb(90, 181, 75, maxColorValue=255), # bam1 = dark green
                 rgb(65, 91, 39, maxColorValue=255), # bam1 bam3 = dark dark green
                 rgb(179, 218, 38, maxColorValue=255), # bam3 = lighter green yellow
                 rgb(64, 64, 64, maxColorValue=255), # Col Prostarv = grey
                 rgb(231, 25, 187, maxColorValue=255), # dpe1 = fel rose
                 rgb(181, 49, 175, maxColorValue=255), # dpe2 = rose
                 rgb(177, 132, 244, maxColorValue=255), # isa1 = violet
                 rgb(32, 32, 244, maxColorValue=255), # mex = blue
                 rgb(245, 62, 31, maxColorValue=255), # pgi = red
                 rgb(21, 92, 171, maxColorValue=255), # sex1 = bleuish
                 rgb(92, 2, 190, maxColorValue=255)) # ss4 = violet
  
  ColorsTrt <- data.frame(idGenotype = c("Col-0", "pgm-1", "abcb14-1", "abcb14-2",
                                         "amy3-2", "amy3-2 bam1", "bam1", "bam1bam3", "bam3",
                                         "Col-0 ProStarv::LUC", "dpe1-2", "dpe2-5", "isa1-1 ProStarv::LUC",
                                         "mex1-1", "pgi1-1", "sex1-3", "ss4-3 ProStarv::LUC"),
                          idWatering = "WW",
                          col = ColorsTrt,
                          lty = 1)
  
  ColorsTrt <- rbind(ColorsTrt, data.frame(idGenotype = c("Col-0", "pgm-1"),
                                           idWatering = "WS",
                                           col = c(ColorsTrt$col[ColorsTrt$idGenotype == "Col-0"], ColorsTrt$col[ColorsTrt$idGenotype == "pgm-1"]),
                                           lty = 2))
  
  ColorsTrt$Trt <- paste(ColorsTrt$idGenotype, ColorsTrt$idWatering, sep = " - ")

  
  t1 <- c(-2, 2, 3, 4, 5)
  t2 <- c( 2, 3, 4, 5, 6)
  ColorsPeriod <- data.frame(idPeriod = c("Control", "High CO2", "Low light", "Recovery 1", "Recovery 2"),
                             Time1 = t1,
                             Time2 = t2,
                             decimalDay1 = Time_ON0 + t1,
                             decimalDay2 = Time_ON0 + t2,
                             #col = c("black", "coral3", "burlywood4", "gray50", "gray30"))
                             col = c("black", "coral3", "burlywood4", "gray30", "gray30"))
  
  
  return (list(ColorsTrt = ColorsTrt, ColorsPeriod = ColorsPeriod))
  }


# Call the function
c(ColorsTrt, ColorsPeriod) := set_colors_C2M47()



#------------------------------------------------------------------------------#
#                             Check irrigation data                            #
#------------------------------------------------------------------------------#

irrig_file <- file.path(here::here(), spcs, idExp, "Processed_data", "C2M47_starch_irrigations.csv")
if (file.exists(irrig_file))
  {
  # Import and manage irrigation data
  irrig <- read.csv(irrig_file)
  irrig$idPot <- sub("C2M47-", "", irrig$idPot)#################
  irrig <- irrig[order(irrig$idPot, irrig$decimalDay), ]
  
  # Check the number of irrigation events for each pot
  aggregate(list(N = irrig$decimalDay), by = list(idPot = irrig$idPot), FUN = length)
  # --> Variable, between 1 and 4 irrigation events depending on the pot (and not the treatment)
  #     It would therefore be confusing to plot an average time per experiment or per condition.
  
  # Nonetheless, many pots were irrigated on day 1 and all pots on day 3 in the morning:
  #hist(irrig$decimalDay - Time_ON0)
  # So for convenience we define a rough 'mean_irrig_time'
  mean_irrig_time <- c(1, 3) + mean(irrig$decimalDay - Time_ON0 - round(irrig$decimalDay - Time_ON0)) + Time_ON0
  # that will be avoided unless it can be specified in the figure legend that it is an approximation of the reality.
  }


#------------------------------------------------------------------------------#
#                        Depict environmental conditions                       #
#------------------------------------------------------------------------------#


# Define and draw upper rectangles representing environmental conditions
# It does not have to be run here, but it will be called by 'prepare_kin()' within 'PhenoLeaks_graphics.R'.
# It should be named as "rectangles_environment_" + the identifier of the experiment ('idExp')

rectangles_environment_C2M47 <- function (col.per = ColorsPeriod, # the dataframe containing the features of the different periods
                                          Time_var, # "Time" or "decimalDay"
                                          irrig_show_mode = "none", # show an arrow at irrigation; use "pot" for the irrigation specific to a given pot, or "none" for no arrow
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
  
  if (!export_PPTX)
    {
    clip(usr[1], usr[2], y1, y2)
    rect(col.per[, paste(Time_var, 1, sep = "")],
         y1,
         col.per[, paste(Time_var, 2, sep = "")],
         y2,
         col = col.per$col, border = NA, density = NULL)
    clip(usr[1], usr[2], usr[3], usr[4])
    }
  
  else
    {
    for (per in 1:nrow(col.per))
      {
      x1 <- col.per[per, paste(Time_var, 1, sep = "")]
      x2 <- col.per[per, paste(Time_var, 2, sep = "")]
      if ((x1+x2)/2 >= usr[1] & (x1+x2)/2 <= usr[2])
        {
        rect(x1, y1, x2, y2,
             col = col.per$col[per], border = NA, density = NULL)
        }
      }
    }
  
  y.text <- (y1 + y2) / 2 - if (export_PPTX) 0.015 else 0
  t.control <- 0.5 + if (Time_var == "decimalDay") Time_ON0 else 0
  if (t.control > usr[1] & t.control < usr[2]) text(t.control, y.text, "Control", cex = cex.env, xpd = T, col = "white")
  t1.high.co2 <- col.per[2, paste(Time_var, 1, sep = "")]
  t1.low.light <- col.per[3, paste(Time_var, 1, sep = "")]
  if (t1.high.co2 > usr[1] & t1.high.co2 < usr[2]) text(t1.high.co2+Pho_Per/24/2, y.text, expression("High CO"[2]), cex = cex.env, xpd = T)
  if (t1.low.light > usr[1] & t1.low.light < usr[2]) text(t1.low.light+Pho_Per/24/2, y.text, "Low light", cex = cex.env, xpd = T)
  segments(c(t1.high.co2, t1.low.light)+Pho_Per/24, y1, c(t1.high.co2, t1.low.light)+Pho_Per/24, y2, col = "gray", lty = "22", xpd = T)
  t.recovery <- 5 + if (Time_var == "decimalDay") Time_ON0 else 0
  if (t.recovery > usr[1] & t.recovery < usr[2]) text(t.recovery, y.text, "Recovery", cex = cex.env, xpd = T, col = "white")
  
  if (irrig_show_mode != "none" & file.exists(irrig_file))
    {
    if (irrig_show_mode == "mean") t_irrig <- mean_irrig_time
    else if (irrig_show_mode == "pot") t_irrig <- irrig$decimalDay[irrig$idPot == pot]
    else stop ("'irrig_show_mode' should be \"none\", \"mean\" or \"pot\".")
    if (Time_var != "decimalDay") t_irrig <- t_irrig - Time_ON0
    text(t_irrig, usr[3] + (if (!export_PPTX) 0.025 else 0.075)*(usr[4] - usr[3]),
         if (!export_PPTX) "ad hoc watering" else "ad hoc\nwatering",
         col = "black", cex = cex.env)
    #require(igraph)
    iArrows <- igraph:::igraph.Arrows
    iArrows(t_irrig, rep(usr[3] + (if (!export_PPTX) 0.05 else 0.15)*(usr[4] - usr[3]), length(t_irrig)),
            t_irrig, rep(usr[3] + (if (!export_PPTX) 0.15 else 0.25)*(usr[4] - usr[3]), length(t_irrig)),
            h.lwd = if (!export_PPTX) 1 else 0.75,
            sh.lwd = if (!export_PPTX) 1 else 0.75,
            sh.col = if (irrig_show_mode != "mean") "black" else "gray",
            curve = 0,
            width = if (!export_PPTX) 1 else 0.75,
            size = if (!export_PPTX) 0.5 else 0.3)
    }
  }

