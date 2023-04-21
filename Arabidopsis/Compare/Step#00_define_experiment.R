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
#                                    Compare                                   #
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
for (pkg in c("scales", "here"))
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

set_constants_Compare <- function ()
  {
  # Plant species
  spcs <- "Arabidopsis"
  
  # Name of the experiment
  idExp <- "Compare"
  
  # Duration of the photoperiod (h)
  Pho_Per <- 12
  
  # Duration of the skotoperiod (h)
  # Note that the scripts will not work for diel cycles that do not match 24 h.
  Sko_Per <- 24 - Pho_Per
  
  # Duration covering the rapid stomatal movements at the transitions (min)
  # Corresponds here to the first data point after a day/night transition.
  Rapid_dur <- 30
  
  # Fraction of days by which the square wave (see PhenoLeaks_fit) is shifted compared to the dark/light cycle (min)
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
set_constants_Compare()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                         (4)  Set graphical parameters                        #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#------------------------------------------------------------------------------#
#            Set color and line type for treatments and environments           #
#------------------------------------------------------------------------------#

# Set the colors and line types within a user-defined function

set_colors_Compare <- function ()
  {
  Colors_C3M31 <- rgb(0, 0, 0, maxColorValue=255)  # Col-0 = black
  Colors_C3M31 <- data.frame(idExperiment = "Experiment #1",
                             idGenotype = "Col-0",
                             idWatering = "WW",
                             col_Geno = Colors_C3M31,
                             lty = 1)
  
  Colors_C2M43A <- rgb(0, 0, 0, maxColorValue=255)  # Col-0 = black
  Colors_C2M43A <- data.frame(idExperiment = "Experiment #2",
                              idGenotype = "Col-0",
                              idWatering = "WW",
                              col_Geno = Colors_C2M43A,
                              lty = 1)
  
  Colors_C2M43B <- c(rgb(0, 0, 0, maxColorValue=255),  # Col-0 = black
                     rgb(237, 125, 49,maxColorValue=255), # pgm = orange
                     rgb(90, 181, 75, maxColorValue=255), # bam1 = dark green
                     rgb(65, 91, 39, maxColorValue=255), # bam1 bam3 = dark dark green
                     rgb(179, 218, 38, maxColorValue=255), # bam3 = lighter green yellow
                     rgb(21, 92, 171, maxColorValue=255)) # sex1 = bleuish
  Colors_C2M43B <- data.frame(idExperiment = "Experiment #3",
                              idGenotype = c("Col-0", "pgm-1", "bam1", "bam1bam3", "bam3", "sex1-3"),
                              idWatering = "WW",
                              col_Geno = Colors_C2M43B,
                              lty = 1)
  
  Colors_C2M47 <- c(rgb(0, 0, 0, maxColorValue=255),  # Col-0 = black
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
  Colors_C2M47 <- data.frame(idExperiment = "Experiment #4",
                             idGenotype = c("Col-0", "pgm-1", "abcb14-1", "abcb14-2",
                                            "amy3-2", "amy3-2 bam1", "bam1", "bam1bam3", "bam3",
                                            "Col-0 ProStarv::LUC", "dpe1-2", "dpe2-5", "isa1-1 ProStarv::LUC",
                                            "mex1-1", "pgi1-1", "sex1-3", "ss4-3 ProStarv::LUC"),
                             idWatering = "WW",
                             col_Geno = Colors_C2M47,
                             lty = 1)
  Colors_C2M47 <- rbind(Colors_C2M47, data.frame(idExperiment = "Experiment #4",
                                                 idGenotype = c("Col-0", "pgm-1"),
                                                 idWatering = "WS",
                                                 col_Geno = c(Colors_C2M47$col[Colors_C2M47$idGenotype == "Col-0"], Colors_C2M47$col[Colors_C2M47$idGenotype == "pgm-1"]),
                                                 lty = 2))
  
  ColorsTrt <- rbind(Colors_C3M31, Colors_C2M43A, Colors_C2M43B, Colors_C2M47)
  
  ColorsTrt$Trt <- paste(ColorsTrt$idExp, ColorsTrt$idGenotype, ColorsTrt$idWatering, sep = " - ")
  
  colors_expe <- hue_pal()(length(unique(ColorsTrt$idExperiment)))
  ColorsTrt$col_Expe <- colors_expe[match(ColorsTrt$idExperiment, unique(ColorsTrt$idExperiment))]
  
  # ONLY SET THE CONTROL PERIOD HERE (the only one that is shared between experiments)
  t1 <- c(-2)
  t2 <- c( 2)
  ColorsPeriod <- data.frame(idPeriod = c("Control"),
                             Time1 = t1,
                             Time2 = t2,
                             col = c("black"))
  
  
  return (list(ColorsTrt = ColorsTrt, ColorsPeriod = ColorsPeriod))
  }


# Call the function
c(ColorsTrt, ColorsPeriod) := set_colors_Compare()



#------------------------------------------------------------------------------#
#                        Depict environmental conditions                       #
#------------------------------------------------------------------------------#


# Define and draw upper rectangles representing environmental conditions
# It does not have to be run here, but it will be called by 'prepare_kin()' within 'PhenoLeaks_graphics.R'.
# It should be named as "rectangles_environment_" + the identifier of the experiment ('idExp')

rectangles_environment_Compare <- function (col.per = ColorsPeriod, # the dataframe containing the features of the different periods
                                            Time_var ="Time",
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
  t.control <- 0.5
  if (t.control > usr[1] & t.control < usr[2]) text(t.control, y.text, "Control", cex = cex.env, xpd = T, col = "white")
  }

