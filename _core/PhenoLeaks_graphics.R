################################################################################
#                                                                              #
#                             PhenoLeaks - GRAPHICS                            #
#                                                                              #
#              Graphical functions to plot the transpiration data              #
#                        and some associated statistics                        #
#                                                                              #
#                             Florent Pantin, 2022                             #
#                                                                              #
################################################################################



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                                (1)  Libraries                                #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#options(repos = "https://cran.rstudio.com/") # RStudio
for (pkg in c("plotrix", "ggplot2", "ggpmisc", "rstatix", "ggpubr", "scales", "viridis"))
  {
  if (!pkg %in% installed.packages()[, "Package"]) { install.packages(pkg) }
  #update.packages(pkg)
  library(pkg, character.only = T)
  }

if (!"export" %in% installed.packages()[, "Package"])
  { 
  if (!"devtools" %in% installed.packages()[, "Package"]) { install.packages("devtools") }
  library("devtools")
  devtools::install_github("tomwenseleers/export")
  }
library("export")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                                (2)  Functions                                #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



#------------------------------------------------------------------------------#
#                                 Plot kinetics                                #
#------------------------------------------------------------------------------#

# Prepare the kinetics plot
prepare_kin <- function (dat,
                         use_VPD = F, # set to TRUE when transpiration is normalized by VPD
                         Time_var = "Time", # the column name of the time axis
                         E_var = if (!use_VPD) "E_mmol_per_m2_s" else "E_mmol_per_m2_s_kPa", # the column name of the transpiration axis
                         xlab = if (x_axis_unit == "days") "Time (d)",
                         ylab = if (!use_VPD) expression(paste("E (mmol m"^-2, " s"^-1, ")", sep = "")) else expression(paste("E / VPD (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                         xlim = NULL,
                         ylim = NULL,
                         main = NULL,
                         add_VPD = use_VPD, # set to TRUE to plot VPD using the right y-axis
                         add_VPD_mode = "single", # if "single", only one VPD kinetics is drawn (averaged if several treatments); if "multiple", VPD kinetics are prepared but have to be drawn by another function (e.g. 'plot_avg_kin()' that draws one VPD kinetics per treatment)
                         ylim_VPD = c(0, 1.5), # the limits of the VPD axis (kPa)
                         add_SWC = F,  # set to TRUE to plot SWC using another right y-axis
                         ylim_SWC = c(0.5, 2), # the limits of the SWC axis (g water / g dry soil)
                         export_PPTX = F, # set to TRUE to export the graph as PPTX instead of PDF
                         draw_x_axis = T, # set to FALSE to draw only the ticks of the x-axis
                         x_axis_unit = "days", # "days or "hours"
                         write_y_axis_label = if (export_PPTX) F else T, # set to FALSE to omit the y-axis label. Useful when expressions are not supported by not supported by 'graph2ppt()', e.g. with superscripts, see https://github.com/tomwenseleers/export/issues/49
                         mar = if (export_PPTX) c(1.75, 2.5, 0.1, 0.1) else if (is.null(main)) c(2.5, 3.5, 0.5, if (add_VPD) 3.5 else 0.5) else c(2.5, 3.5, 2.5, if (add_VPD) 3.5 else 0.5), # plot margins
                         light_ON = if (Time_var == "decimalDay") Time_ON0 - floor(Time_ON0) else 0, # time of the diel cycle when the light switches on
                         night_rect_xleft = NULL, # vector of the x-positions for the left side of the grey rectangles (nighttime). Set to NULL for automatic positioning. Set to NA to ignore nighttime rectangles.
                         night_rect_xright = NULL, # vector of the x-positions for the right side of the grey rectangles (nighttime). Set to NULL for automatic positioning. Set to NA to ignore nighttime rectangles.
                         y.axis.break = if (par("usr")[3] < 0) NULL else par("usr")[3] + 0.4*(axTicks(2)[1] - par("usr")[3]), # position on the y-axis where a zigzag break should be drawn (optional),
                         draw_rect_env = T, # set to TRUE to draw upper rectangles representing environmental conditions
                         idExperiment = idExp, # identifier of the experiment, e.g."C2M47". Required if 'draw_rect_env' is TRUE.
                         ...) # other arguments to be passed to the 'rectangles_environment_[idExp]()' function
  {
  # x-limits
  if (is.null(xlim)) { xlim <- c(min(dat[, Time_var], na.rm = T), max(dat[, Time_var], na.rm = T)) }
 
  # y-limits - if only NAs
  if (length(which(!is.na(dat[, E_var]))) == 0) { ylim <- c(0, 1) }

  # Margins
  par(mar = mar)
  
  # Core plot
  plot(get(E_var) ~ get(Time_var), data = dat, type = "n", bty = "n",
       xlim = xlim, ylim = ylim,
       xlab = NA, ylab = NA, xaxt = "n", yaxt = "n")
  if (!is.null(main)) title(main, line = 1.5)
      
  # Axes
  if (draw_x_axis)
    {
    if (x_axis_unit == "hours")
      {
      axis(1, at = seq(0, 24, by = 2)/24, labels = seq(0, 24, by = 2), padj = if (export_PPTX) -2 else -1.25, tck = -0.025, cex.axis = 0.7)
      }
    else
      {
      axis(1, padj = if (export_PPTX) -2 else -1.25, tck = -0.025, cex.axis = 0.7)
      }
    mtext(xlab, side = 1, line = if (export_PPTX) (if (x_axis_unit == "hours") 1 else 0.7) else 1.25, cex = 0.9)
    }
  else
    {
    axis(1, at = axTicks(1), labels = rep("", length(axTicks(1))), tck = -0.025)
    }
  axis(2, hadj = if (export_PPTX) (if (x_axis_unit == "hours") 0.3 else 0.1) else NA, tck = -0.025, las = 1, cex.axis = 0.7)  
  if (write_y_axis_label)
    {
    mtext(ylab, side = 2, line = if (export_PPTX) 1.25 else 2, cex = 0.9)
    }
  
  # Draw nighttime rectangles
  usr <- par("usr")
  if (is.null(night_rect_xleft))
    {
    rect1 <- if (xlim[1] - light_ON == round(xlim[1] - light_ON)) xlim[1] else floor(xlim[1]) + light_ON
    if (usr[1] < rect1) rect1 <- rect1 - 1
    rectn <- if (xlim[2] - light_ON == round(xlim[2] - light_ON)) xlim[2] - 1 else ceiling(xlim[2]) + light_ON - 1
    if (usr[2] > rectn + 1 + Pho_Per/24) rectn <- rectn + 1
    night_rect_xleft <- rect1:rectn + Pho_Per/24
    }
  ybottom <- rep(usr[3], length(night_rect_xleft))
  if (is.null(night_rect_xright))
    {
    night_rect_xright <- night_rect_xleft + Sko_Per/24
    }
  ytop <- rep(usr[4], length(night_rect_xleft))

  if (!is.na(night_rect_xleft[1]) & !is.na(night_rect_xright[1]))
    {
    if (length(night_rect_xleft) != length(night_rect_xright))
      {
      stop ("'night_rect_xleft' and 'night_rect_xright' should have the same length.")
      }
    rect(night_rect_xleft, ybottom, night_rect_xright, ytop, col = rgb(230, 230, 230, maxColorValue = 255), border = NA, density = NULL)
    }
  
  # Draw upper rectangles to display environmental conditions
  # The function should be defined for each experiment on step#0 (e.g. 'Step#00_define_experiment.R') 
  if (draw_rect_env)
    {
    do.call(get(paste("rectangles_environment_", idExperiment, sep = "")), list(export_PPTX = export_PPTX, Time_var = Time_var, ...))
    }
  
  # Re-draw the box
  box()
  
  # Zigzag axis break
  #require(plotrix)
  if (!is.null(y.axis.break))
    {
    axis.break(2, y.axis.break, style = "zigzag")
    }

  # VPD plot (if required)
  if (add_VPD)
    {
    # SWC plot (if required).
    # For the moment only implemented with 'add_VPD' == TRUE, and for one pot
    if (add_SWC)
      {
      par(new = TRUE)
      plot(SWC ~ get(Time_var), data = dat, type = "l", bty = "n",
           xlim = xlim, ylim = ylim_SWC,
           xlab = NA, ylab = NA, xaxt = "n", yaxt = "n", col = "orange")
      axis(4, line = 4, hadj = if (export_PPTX) 0.1 else NA, tck = -0.025, las = 1, cex.axis = 0.7, col = "orange", col.axis = "orange")
      mtext(expression(paste("SWC (g"["water"], " g"["dry soil"]^"-1", ")", sep = "")), side = 4, line = 4.25 + if (export_PPTX) 1.25 else 2, cex = 0.9, col = "orange")
      }

    par(new = TRUE)    
    if (add_VPD_mode == "single")
      {
      dat_vpd <- aggregate(list(meanVPD = dat$meanVPD), by = list(Time = dat[, Time_var]), FUN = mean, na.rm = T) # no difference if only one treatment (except that NAs will be removed)
      plot(meanVPD ~ Time, data = dat_vpd, type = "l", col = "grey", bty = "n",
           xlim = xlim, ylim = ylim_VPD,
           xlab = NA, ylab = NA, xaxt = "n", yaxt = "n")
      axis(4, hadj = if (export_PPTX) 0.1 else NA, tck = -0.025, las = 1, cex.axis = 0.7, col = "grey", col.axis = "grey")
      mtext("VPD (kPa)", side = 4, line = if (export_PPTX) 1.25 else 2, cex = 0.9, col = "grey")
      
      # Re-plot empty transpiration so as to go back to the correct coordinates
      par(new = TRUE)
      plot(get(E_var) ~ get(Time_var), data = dat, type = "n", bty = "n",
           xlim = xlim, ylim = ylim,
           xlab = NA, ylab = NA, xaxt = "n", yaxt = "n")
      }
    else if (add_VPD_mode == "multiple")
      {
      plot(meanVPD ~ get(Time_var), data = dat, type = "n", bty = "n",
          xlim = xlim, ylim = ylim_VPD,
          xlab = NA, ylab = NA, xaxt = "n", yaxt = "n")
      axis(4, hadj = if (export_PPTX) 0.1 else NA, tck = -0.025, las = 1, cex.axis = 0.7, col = "grey", col.axis = "grey")
      mtext("VPD (kPa)", side = 4, line = if (export_PPTX) 1.25 else 2, cex = 0.9, col = "grey")
      }
    else
      {
      stop ("Argument 'add_VPD_mode' should be \"single\" or \"multiple\"")
      }
    }
  
  return (plot.status = if (add_VPD & add_VPD_mode == "multiple") "VPD" else "E")
  }


# Plot the average kinetics and error interval
plot_avg_kin <- function (avg_dat, # the average dataframe
                          treatments, # the (concatenated) name of the treatment(s) to be plotted, e.g. "Col-0 - WW"
                          Trt_var = "Trt", # the column name of the (concatenated) statistical treatment(s), e.g. "Trt" (or "idGenotype" if only one)
                          Time_var = "Time", # the column name of the time axis
                          E_var = if ("E_corr_per_kPa_offset_avg" %in% names(avg_dat)) "E_corr_per_kPa_offset_avg" else "E_corr_offset_avg", # the column name of the transpiration axis
                          error_bar = "SE", # either "SE", "SD" or "CI"
                          ignore_NA = F, # set to FALSE will interrupt the error interval where the error bar is NA
                          xlim = c(min(avg_dat[, Time_var], na.rm = T), max(avg_dat[, Time_var], na.rm = T)),
                          ylim = c(0.83*min(avg_dat[, E_var], na.rm = T), 1.06*max(avg_dat[, E_var], na.rm = T)),
                          add_VPD_mode = "multiple",
                          text.legend = trt, # the text for the legend
                          x.legend = par("usr")[1] + 0.521*(par("usr")[2]-par("usr")[1]),
                          y.legend = par("usr")[3] + 0.932*(par("usr")[4]-par("usr")[3]),
                          cex.legend = 0.9,
                          color_E_by_SWC = F, # set to TRUE to color transpiration by the soil water content
                          col_by_SWC = rep("black", length(treatments)), # color of the (semi-transparent) polygon of the error bar around transpiration when 'color_E_by_SWC' is TRUE
                          lty_by_SWC = rep(1, length(treatments)), # line type when 'color_E_by_SWC' is TRUE
                          title_by_SWC = "", # title for the legend when 'color_E_by_SWC' is TRUE
                          scale_bar_pos_by_SWC = c(x.legend + 0.5, y.legend - 0.17*(par("usr")[4]-par("usr")[3])), # position of the scale bar when 'color_E_by_SWC' is TRUE
                          ...) # other arguments to be passed to the 'prepare_kin()' function
  {
  
  # Select the data
  dat.trt <- avg_dat[avg_dat[, Trt_var] %in% treatments, ]

  # Prepare the plot
  plot.status <- prepare_kin(dat.trt, Time_var = Time_var, E_var = E_var, xlim = xlim, ylim = ylim,
                             #main = main, export_PPTX = export_PPTX,
                             add_VPD_mode = add_VPD_mode, ...)
 
  # VPD plot (if required)
  if (plot.status == "VPD")
    {
    for (trt in treatments)
      {
      dat <- dat.trt[dat.trt[, Trt_var] == trt, ]
      color =  ColorsTrt$col[match(trt, ColorsTrt[, Trt_var])]
      lty =  ColorsTrt$lty[match(trt, ColorsTrt[, Trt_var])]
      points(meanVPD ~ get(Time_var), data = dat, type = "l", lty = lty,
             col = rgb(t(col2rgb(color))/255, alpha = 0.7))
      rm(dat)
      }
    par(new = TRUE)
    plot(get(E_var) ~ get(Time_var), data = dat.trt, type = "n", bty = "n",
         xlim = xlim, ylim = ylim,
         xlab = NA, ylab = NA, xaxt = "n", yaxt = "n")
    }

  # Add the average and error interval of the treatments
  for (trt in treatments)
    {
    dat <- dat.trt[dat.trt[, Trt_var] == trt, ]
    
    if (!color_E_by_SWC)
      {
      color <- ColorsTrt$col[match(trt, ColorsTrt[, Trt_var])]
      lty <- ColorsTrt$lty[match(trt, ColorsTrt[, Trt_var])]
      }
    else
      {
      color <- col_by_SWC[which(treatments %in% trt)]
      lty <- lty_by_SWC[which(treatments %in% trt)]
      }
    if (ignore_NA)
      {
      xpol <- c(dat[, Time_var], rev(dat[, Time_var]))
      ypol <- c(dat[, E_var] + dat[, error_bar], rev(dat[, E_var] - dat[, error_bar]))
      xpol <- xpol[!is.na(ypol)]; ypol <- ypol[!is.na(ypol)]
      polygon(xpol, ypol, border = NA, col = rgb(t(col2rgb(color))/255, alpha = 0.3))
      }
    else # adapted from https://stackoverflow.com/questions/33372389/how-to-draw-a-polygon-around-na-values-in-r
      {
      enc <- rle(!is.na(dat[, error_bar]))
      endIdxs <- cumsum(enc$lengths)
      for(i in 1:length(enc$lengths))
        {
        if(enc$values[i])
          {
          endIdx <- endIdxs[i]
          startIdx <- endIdx - enc$lengths[i] + 1
          subtranspi <- dat[startIdx:endIdx, E_var]
          suberror <- dat[startIdx:endIdx, error_bar]
          subtime <- dat[startIdx:endIdx, Time_var]
          xpol <- c(subtime, rev(subtime))
          ypol <- c(subtranspi - suberror, rev(subtranspi + suberror))
          polygon(xpol, ypol, border = NA, col = rgb(t(col2rgb(color))/255, alpha = 0.3))
          }
        }  
      }
      
    if (!color_E_by_SWC)
      {
      points(get(E_var) ~ get(Time_var), data = dat, type = "l", lwd = 2, lty = lty,
             col = rgb(t(col2rgb(color))/255, alpha = 0.7))
      }
    else
      {
      #require("viridis")
      ncol <- 100
      col_palette <- rev(plasma(ncol))
      miniSWC <- 0.5 # g water / g dry soil
      maxiSWC <- 1.5 # g water / g dry soil
      x <- dat[, Time_var]
      y <- dat[, E_var]
      z <- dat[, "SWC"]
      xy <- approx(x, y, n = 1000, na.rm = F)
      xz <- approx(x, z, n = 1000, na.rm = F); names(xz)[2] <- "z"
      color_interpolate <- col_palette[pmax(pmin(1 + (xz$z - miniSWC)/maxiSWC * (ncol-1), ncol), 1)]
      segments(xy$x[1:length(xy$x)-1],
               xy$y[1:length(xy$x)-1],
               xy$x[2:length(xy$x)],
               xy$y[2:length(xy$x)],
               col = rgb(t(col2rgb(color_interpolate[1:length(xy$x)-1]))/255, alpha = 0.5),
               lwd = 4)
      points(get(E_var) ~ get(Time_var), data = dat, type = "l", lwd = 1, lty = lty,
             col = rgb(t(col2rgb(color))/255, alpha = 0.7))
      }
    rm(dat)
    }
  
  # Add the legend
  if (!color_E_by_SWC)
    {
    legend(x.legend, y.legend, legend = text.legend, cex = cex.legend,
           bty = "n", pch = NA,
           lty = ColorsTrt$lty[match(treatments, ColorsTrt[, Trt_var])], lwd = 2,
           col = ColorsTrt$col[match(treatments, ColorsTrt[, Trt_var])])
    }
  else
    {
    legend(x.legend, y.legend + 0.1, legend = text.legend, cex = cex.legend,
           bty = "n", pch = NA,
           lty = lty_by_SWC, lwd = 1,
           col = col_by_SWC,
           title = title_by_SWC, title.adj = 0)
    
    x1 <- scale_bar_pos_by_SWC[1]
    x2 <- x1 + 1
    y1 <- scale_bar_pos_by_SWC[2]
    y2 <- y1 + 0.05*(par("usr")[4]-par("usr")[3])
    slo <- (x2-x1)/ncol
    int <- x1
    #rect((0:(ncol-1))*slo + int, y1, (1:ncol)*slo + int, y2, col = rgb(t(col2rgb(col_palette[1:ncol]))/255, alpha = 0.5), border = NA, xpd = NA)
    rect((0:(ncol-1))*slo + int, y1, (1:ncol)*slo + int, y2, col = col_palette[1:ncol], border = NA, xpd = NA)
    rect(x1, y1, x2, y2, lwd = 0.5, border = "black")
    text((x1+x2)/2, y1, expression(paste("SWC (g"["water"], " g"["dry soil"]^"-1", ")", sep = "")), cex = cex.legend, col = "black", pos = 1)
    #text(x1,(y1+y2)/2, paste("\u2264", miniSWC, sep = " "), cex = cex.legend, pos = 2, col = "black") # to be solved: encoding problem, the symbol is not properly displayed (either a square or equal)
    text(x1,(y1+y2)/2, miniSWC, cex = cex.legend, pos = 2, col = "black")
    #text(x2,(y1+y2)/2, paste("\u2265", maxiSWC, sep = " "), cex = cex.legend, pos = 4, col = "black") # to be solved: encoding problem, the symbol is not properly displayed (either a square or equal)
    text(x2,(y1+y2)/2, maxiSWC, cex = cex.legend, pos = 4, col = "black")
    ticksScale <- c(0.75, 1, 1.25)
    ticksScale <- (ticksScale-miniSWC)/(maxiSWC-miniSWC)*ncol*slo + int
    segments(ticksScale, y1, ticksScale, y1+0.15*(y2-y1), lwd = 0.5, col = "black")
    segments(ticksScale, y2, ticksScale, y2-0.15*(y2-y1), lwd = 0.5, col = "black")
    }
  }



#------------------------------------------------------------------------------#
#                     Compute average kinetics by treatment                    #
#------------------------------------------------------------------------------#

# Compute the average (and error intervals)
compute_avg <- function (all_dat, # the whole dataframe
                         E_var = "E_corr",
                         Trt_var = "Trt", # the column name of the (concatenated) statistical treatment(s), e.g. "Trt" (or "idGenotype" if only one)
                         ALPHA = 0.05, # alpha risk for computing the confidence interval
                         avg_only = F) # set to TRUE to compute only the average, not the error bars
  {
  sortlist <- list(Exact_Time_min = all_dat$Exact_Time_min, Trt = all_dat[, Trt_var])

  if (avg_only)
    {
    avg_dat <- do.call(data.frame, aggregate(list(E_avg = all_dat[, E_var]),
                                             by = sortlist,
                                             FUN = function (x) AVG = mean(x, na.rm = T)))
    names(avg_dat) <- c("Exact_Time_min", Trt_var, paste(E_var, "_avg", sep = ""))
    }
  
  else
    {
    avg_dat <- do.call(data.frame, aggregate(list(E_avg = all_dat[, E_var]),
                                             by = sortlist,
                                             FUN = function (x) c(AVG = mean(x, na.rm = T),
                                                                  SD = sd(x, na.rm = T),
                                                                  N = length(na.omit(x)),
                                                                  SE = sd(x, na.rm = T) / sqrt(length(na.omit(x)) - 1),
                                                                  CI = qt((1 - ALPHA/2), length(na.omit(x)) - 1) * sd(x, na.rm = T)/sqrt(length(na.omit(x))))))
    
    names(avg_dat) <- c("Exact_Time_min", Trt_var, paste(E_var, "_avg", sep = ""), "SD", "N", "SE", "CI")
    }
  
  return (avg_dat)
  }


# Correct the average for missing values to avoid artifactual "jumps" in the kinetics due to punctual removal of outliers
offset_avg <- function (all_dat, # the whole dataframe
                        avg_dat, # the average dataframe
                        interval_est = 0.15) # time interval before and after the missing values to compute the correction (d)
  {
  Trt_var <- names(avg_dat)[2] # "Trt"
  E_var_avg <- names(avg_dat)[3] # "E_corr_avg"
  E_var <- sub("_avg", "", E_var_avg) # "E_corr"

  all_dat$DAY_or_NIGHT <- "DAY"
  all_dat$DAY_or_NIGHT[all_dat$Exact_Time_min/60/24 == round(all_dat$Exact_Time_min/60/24) | all_dat$Exact_Time_min/60/24 - floor(all_dat$Exact_Time_min/60/24) > Pho_Per/24] <- "NIGHT"
  
  all_dat[, paste(E_var, "_offset", sep = "")] <- all_dat[, E_var]
  
  for (pot in unique(all_dat$idPot))
    {
    dat <- all_dat[all_dat$idPot == pot, ]
    trt <- dat[1, Trt_var]
    missing.values <- which(is.na(dat[, E_var]))
    if (length(missing.values) > 0)
      {
      ifelse (length(missing.values) == 1,
              missing.groups <- 1,
              missing.groups <- c(1, which(diff(missing.values) != 1)+1))
      missing.groups <- c(missing.groups, length(missing.values)+1)
      for (miss.grp in 1:(length(missing.groups)-1))
        {
        miss.val <- missing.values[missing.groups[miss.grp]:(missing.groups[miss.grp+1]-1)]
        times <- dat$Exact_Time_min[miss.val]
        dn <- unique(dat$DAY_or_NIGHT[miss.val])
        
        dat_trt <- all_dat[all_dat[, Trt_var] == trt, ]
        local_dat <- dat_trt[dat_trt$Exact_Time_min > (min(times) - interval_est*24*60) &
                               dat_trt$Exact_Time_min < (max(times) + interval_est*24*60) &
                               dat_trt$DAY_or_NIGHT %in% dn, ]
        
        if (length(which(!is.na(local_dat[local_dat$idPot == pot, E_var]))) == 0)
          {
          # Release the day/night constraint.
          # Necessary if all the values of a photoperiod (resp. skotoperiod) have been discarded
          # while the previous and next points in the previous and next skotoperiods (resp. photoperiods)
          # have been kept (unadvised, but may be unavoidable).
          local_dat <- dat_trt[dat_trt$Exact_Time_min > (min(times) - interval_est*24*60) &
                                 dat_trt$Exact_Time_min < (max(times) + interval_est*24*60), ]
          if (length(which(!is.na(local_dat[local_dat$idPot == pot, E_var]))) == 0)
            {
            stop (paste("No data to estimate the missing values.\nCheck pot ", pot, " around time ", round(mean(times/60/24), 2), "d.", sep = ""))
            }
          # else
          #   {
          #   warning (paste("The day/night constraint has been released to estimate missing values.\nTHIS IS NOT ADVISED.\nCheck pot ", pot, " around time ", round(mean(times/60/24), 2), "d.", sep = ""))
          #   }
          }
        
        not_missing <- aggregate(list(N = local_dat[, E_var]), by = list(Exact_Time_min = local_dat$Exact_Time_min), FUN = function (x) length(na.omit(x)))
        times_with_at_least_one_data <- which(not_missing$Exact_Time_min %in% times & not_missing$N > 0)
        if (length(times_with_at_least_one_data) > 0)
          {
          local_dat$TimeFactor <- as.factor(local_dat$Exact_Time_min)
          local_dat$idPotFactor <- as.factor(local_dat$idPot)
          mod <- lm(get(E_var) ~ idPotFactor + TimeFactor, data = local_dat)
          pred.times <- as.factor(not_missing$Exact_Time_min[times_with_at_least_one_data])
          pred <- predict(mod, newdata = data.frame(idPotFactor = as.factor(pot), TimeFactor = pred.times))
          all_dat[all_dat$idPot == pot & all_dat$Exact_Time_min %in% pred.times, paste(E_var, "_offset", sep = "")] <- pred
          }
        }
      }
    }
  
  avg_offset <- compute_avg(all_dat, paste(E_var, "_offset", sep = ""), Trt_var, avg_only = T)
  n.Tr <- nrow(avg_dat)
  avg_dat <- merge(avg_dat, avg_offset, sort = F)
  if (nrow(avg_dat) != n.Tr) stop ("Inconsistent number of offset data ---> debug")
  
  # Retrieve time in days
  avg_dat$Time <- avg_dat$Exact_Time_min / 60 / 24
  
  # Return
  return (avg_dat)
  }


# Add the average VPD to an existing dataframe of average transpiration
append_VPD_avg <- function (all_dat, # the whole dataframe
                            avg_dat) # the average dataframe
  {
  Time_var <- names(avg_dat)[1] # "Time"
  Trt_var <- names(avg_dat)[2] # "Trt"

  VPD_avg <- compute_avg(all_dat, E_var = "meanVPD", Trt_var, avg_only = T)
  n.Tr <- nrow(avg_dat)
  avg_dat <- merge(avg_dat, VPD_avg, sort = F)
  if (nrow(avg_dat) != n.Tr) stop ("Inconsistent number of VPD data ---> debug")
  names(avg_dat)[names(avg_dat) == "meanVPD_avg"] <- "meanVPD"
  return(avg_dat)
  }


# Add the average SWC to an existing dataframe of average transpiration
append_SWC_avg <- function (all_dat, # the whole dataframe
                            avg_dat) # the average dataframe
  {
  Time_var <- names(avg_dat)[1] # "Time"
  Trt_var <- names(avg_dat)[2] # "Trt"
  
  SWC_avg <- compute_avg(all_dat, E_var = "SWC", Trt_var, avg_only = T)
  n.Tr <- nrow(avg_dat)
  avg_dat <- merge(avg_dat, SWC_avg, sort = F)
  if (nrow(avg_dat) != n.Tr) stop ("Inconsistent number of SWC data ---> debug")
  names(avg_dat)[names(avg_dat) == "SWC_avg"] <- "SWC"
  return(avg_dat)
  }


#------------------------------------------------------------------------------#
#                         Display statistical analyses                         #
#------------------------------------------------------------------------------#

# Transform a numeric p-value into an asterisk code
transform_pval <- function (pval)
  {
  ifelse(pval > 0.05, pval <- "ns", ifelse(pval > 0.01, pval <- "*", ifelse(pval > 0.001, pval <- "**", ifelse(pval > 0.0001, pval <- "***", pval <- "****"))))
  return (pval)
  }


# Get the string of the ANOVA effect sizes (generalized eta squared)
get_eta2g <- function (table.aov, # the ANOVA table returned by 'get_anova_table()'
                       export_PPTX = F)
  {

  # Main effect(s)
  mono_G <- mono_E <- mono_Expe <- mono_Irr <- ""
  pval_G <- pval_E <- pval_Expe <- pval_Irr <- ""
  ges_G <- ges_E <- ges_Expe <- ges_Irr <- ""
  if ("idGenotype" %in% table.aov$Effect)
    {
    mono_G <- " Geno "
    G <- which(table.aov$Effect == "idGenotype")
    pval_G <- transform_pval(table.aov$p[G])
    ges_G <- round(table.aov$ges[G], 2)
    }
  if ("idPeriod" %in% table.aov$Effect)
    {
    mono_E <- " Env "
    E <- which(table.aov$Effect == "idPeriod")
    pval_E <- transform_pval(table.aov$p[E])
    ges_E <- round(table.aov$ges[E], 2)
    }    
  if ("idExperiment" %in% table.aov$Effect)
    {
    mono_Expe <- " Expe "
    Expe <- which(table.aov$Effect == "idExperiment")
    pval_Expe <- transform_pval(table.aov$p[Expe])
    ges_Expe <- round(table.aov$ges[Expe], 2)
    }
  if ("idWatering" %in% table.aov$Effect)
    {
    mono_Irr <- " Irr "
    Irr <- which(table.aov$Effect == "idWatering")
    pval_Irr <- transform_pval(table.aov$p[Irr])
    ges_Irr <- round(table.aov$ges[Irr], 2)
    }
  
  # Double interaction(s)
  double_GE <- double_GExpe <- double_GIrr <- double_EExpe <- double_EIrr <- double_ExpeIrr <- ""
  pval_GE <- pval_GExpe <- pval_GIrr <- pval_EExpe <- pval_EIrr <- pval_ExpeIrr <- ""
  ges_GE <- ges_GExpe <- ges_GIrr <- ges_EExpe <- ges_EIrr <- ges_ExpeIrr <- ""
  if (any(c("idGenotype:idPeriod", "idPeriod:idGenotype") %in% table.aov$Effect))
    {
    double_GE <- " Geno×Env "
    GE <- which(table.aov$Effect %in% c("idGenotype:idPeriod", "idPeriod:idGenotype"))
    pval_GE <- transform_pval(table.aov$p[GE])
    ges_GE <- round(table.aov$ges[GE], 2)
    }
  if (any(c("idGenotype:idExperiment", "idExperiment:idGenotype") %in% table.aov$Effect))
    {
    double_GExpe <- " Geno×Expe "
    GExpe <- which(table.aov$Effect %in% c("idGenotype:idExperiment", "idExperiment:idGenotype"))
    pval_GExpe <- transform_pval(table.aov$p[GExpe])
    ges_GExpe <- round(table.aov$ges[GExpe], 2)
    }
  if (any(c("idGenotype:idWatering", "idWatering:idGenotype") %in% table.aov$Effect))
    {
    double_GIrr <- " Geno×Irr "
    GIrr <- which(table.aov$Effect %in% c("idGenotype:idWatering", "idWatering:idGenotype"))
    pval_GIrr <- transform_pval(table.aov$p[GIrr])
    ges_GIrr <- round(table.aov$ges[GIrr], 2)
    }
  if (any(c("idPeriod:idExperiment", "idExperiment:idPeriod") %in% table.aov$Effect))
    {
    double_EExpe <- " Env×Expe "
    EExpe <- which(table.aov$Effect %in% c("idPeriod:idExperiment", "idExperiment:idPeriod"))
    pval_EExpe <- transform_pval(table.aov$p[EExpe])
    ges_EExpe <- round(table.aov$ges[EExpe], 2)
    }  
  if (any(c("idPeriod:idWatering", "idWatering:idPeriod") %in% table.aov$Effect))
    {
    double_EIrr <- " Env×Irr "
    EIrr <- which(table.aov$Effect %in% c("idPeriod:idWatering", "idWatering:idPeriod"))
    pval_EIrr <- transform_pval(table.aov$p[EIrr])
    ges_EIrr <- round(table.aov$ges[EIrr], 2)
    }
  if (any(c("idWatering:idExperiment", "idExperiment:idWatering") %in% table.aov$Effect))
    {
    double_ExpeIrr <- " Expe×Irr "
    ExpeIrr <- which(table.aov$Effect %in% c("idWatering:idExperiment", "idExperiment:idWatering"))
    pval_ExpeIrr <- transform_pval(table.aov$p[ExpeIrr])
    ges_ExpeIrr <- round(table.aov$ges[ExpeIrr], 2)
    }    
  
  # Triple interaction
  triple_GEExpe <- triple_GEIrr <- triple_GExpeIrr <- triple_EExpeIrr <- ""
  pval_GEExpe <- pval_GEIrr <- pval_GExpeIrr <- pval_EExpeIrr <- ""
  ges_GEExpe <- ges_GEIrr <- ges_GExpeIrr <- ges_EExpeIrr <- ""
  case_GEExpe <- c("idGenotype:idPeriod:idExperiment", "idGenotype:idExperiment:idPeriod",
                   "idPeriod:idGenotype:idExperiment", "idPeriod:idExperiment:idGenotype",
                   "idExperiment:idGenotype:idPeriod", "idExperiment:idPeriod:idGenotype")
  if (any(case_GEExpe %in% table.aov$Effect))
    {
    triple_GEExpe <- " Geno×Env×Expe "
    GEExpe <- which(table.aov$Effect %in% case_GEExpe)
    pval_GEExpe <- transform_pval(table.aov$p[GEExpe])
    ges_GEExpe <- round(table.aov$ges[GEExpe], 2)
    }
  case_GEIrr <- c("idGenotype:idPeriod:idWatering", "idGenotype:idWatering:idPeriod",
                  "idPeriod:idGenotype:idWatering", "idPeriod:idWatering:idGenotype",
                  "idWatering:idGenotype:idPeriod", "idWatering:idPeriod:idGenotype")
  if (any(case_GEIrr %in% table.aov$Effect))
    {
    triple_GEIrr <- " Geno×Env×Irr "
    GEIrr <- which(table.aov$Effect %in% case_GEIrr)
    pval_GEIrr <- transform_pval(table.aov$p[GEIrr])
    ges_GEIrr <- round(table.aov$ges[GEIrr], 2)
    }
  case_GExpeIrr <- c("idGenotype:idExperiment:idWatering", "idGenotype:idWatering:idExperiment",
                     "idExperiment:idGenotype:idWatering", "idExperiment:idWatering:idGenotype",
                     "idWatering:idGenotype:idExperiment", "idWatering:idExperiment:idGenotype")
  if (any(case_GExpeIrr %in% table.aov$Effect))
    {
    triple_GExpeIrr <- " Geno×Expe×Irr "
    GExpeIrr <- which(table.aov$Effect %in% case_GExpeIrr)
    pval_GExpeIrr <- transform_pval(table.aov$p[GExpeIrr])
    ges_GExpeIrr <- round(table.aov$ges[GExpeIrr], 2)
    }
  case_EExpeIrr <- c("idWatering:idPeriod:idExperiment", "idWatering:idExperiment:idPeriod",
                     "idPeriod:idWatering:idExperiment", "idPeriod:idExperiment:idWatering",
                     "idExperiment:idWatering:idPeriod", "idExperiment:idPeriod:idWatering")
  if (any(case_EExpeIrr %in% table.aov$Effect))
    {
    triple_EExpeIrr <- " Env×Expe×Irr "
    EExpeIrr <- which(table.aov$Effect %in% case_EExpeIrr)
    pval_EExpeIrr <- transform_pval(table.aov$p[EExpeIrr])
    ges_EExpeIrr <- round(table.aov$ges[EExpeIrr], 2)
    }

  expr.aov <- bquote("h"["g"]^2 * .(mono_G) * .(ges_G)^.(pval_G) *
                       .(if (mono_G == "" | mono_E == "") "" else ",") * .(mono_E) * .(ges_E)^.(pval_E) *
                       .(if (mono_E == "" | mono_Expe == "") "" else ",") * .(mono_Expe) * .(ges_Expe)^.(pval_Expe) *
                       .(if (mono_Expe == "" | mono_Irr == "") "" else ",") * .(mono_Irr) * .(ges_Irr)^.(pval_Irr) *
                       .(if (double_GE == "") "" else ",") * .(double_GE) * .(ges_GE)^.(pval_GE) *
                       .(if (double_GExpe == "") "" else ",") * .(double_GExpe) * .(ges_GExpe)^.(pval_GExpe) *
                       .(if (double_GIrr == "") "" else ",") * .(double_GIrr) * .(ges_GIrr)^.(pval_GIrr) *
                       .(if (double_EExpe == "") "" else ",") * .(double_EExpe) * .(ges_EExpe)^.(pval_EExpe) *
                       .(if (double_EIrr == "") "" else ",") * .(double_EIrr) * .(ges_EIrr)^.(pval_EIrr) *
                       .(if (double_ExpeIrr == "") "" else ",") * .(double_ExpeIrr) * .(ges_ExpeIrr)^.(pval_ExpeIrr) *
                       .(if (triple_GEExpe == "") "" else ",") * .(triple_GEExpe)* .(ges_GEExpe)^.(pval_GEExpe) *
                       .(if (triple_GEIrr == "") "" else ",") * .(triple_GEIrr)* .(ges_GEIrr)^.(pval_GEIrr) *
                       .(if (triple_GExpeIrr == "") "" else ",") * .(triple_GExpeIrr)* .(ges_GExpeIrr)^.(pval_GExpeIrr) *
                       .(if (triple_EExpeIrr == "") "" else ",") * .(triple_EExpeIrr)* .(ges_EExpeIrr)^.(pval_EExpeIrr))

  # Special cases more suitable for PPTX export
  if (export_PPTX & nrow(table.aov) == 3 & "idGenotype" %in% table.aov$Effect & "idPeriod" %in% table.aov$Effect)
    {
    expr.aov <- bquote("h"["g"]^2 * " (G, E, G×E) = " * .(ges_G)^.(pval_G) * ", " * .(ges_E)^.(pval_E) * ", " * .(ges_GE)^.(pval_GE))
    }
  
  else if (export_PPTX & nrow(table.aov) == 7 & "idGenotype" %in% table.aov$Effect & "idPeriod" %in% table.aov$Effect & "idWatering" %in% table.aov$Effect)
    {
    # see https://stackoverflow.com/questions/50126447/reducing-spacing-between-lines-when-using-atop
    expr.aov <- bquote(atop(NA, atop(textstyle("h"["g"]^2 * .(mono_G) * .(ges_G)^.(pval_G) * "," * 
                                                 .(mono_E) * .(ges_E)^.(pval_E) * "," * 
                                                 .(mono_Irr) * .(ges_Irr)^.(pval_Irr) * ","),
                                     textstyle(.(double_GE) * .(ges_GE)^.(pval_GE) * "," * 
                                                 .(double_GIrr) * .(ges_GIrr)^.(pval_GIrr) * "," * 
                                                 .(double_EIrr) * .(ges_EIrr)^.(pval_EIrr) * "," * 
                                                 .(triple_GEIrr)* .(ges_GEIrr)^.(pval_GEIrr)))))
    } 

  else if (export_PPTX & nrow(table.aov) == 7 & "idGenotype" %in% table.aov$Effect & "idPeriod" %in% table.aov$Effect & "idExperiment" %in% table.aov$Effect)
    {
    # see https://stackoverflow.com/questions/50126447/reducing-spacing-between-lines-when-using-atop
    expr.aov <- bquote(atop(NA, atop(textstyle("h"["g"]^2 * .(mono_G) * .(ges_G)^.(pval_G) * "," * 
                                                 .(mono_E) * .(ges_E)^.(pval_E) * "," * 
                                                 .(mono_Expe) * .(ges_Expe)^.(pval_Expe) * ","),
                                     textstyle(.(double_GE) * .(ges_GE)^.(pval_GE) * "," * 
                                                 .(double_GExpe) * .(ges_GExpe)^.(pval_GExpe) * "," * 
                                                 .(double_EExpe) * .(ges_EExpe)^.(pval_EExpe) * "," * 
                                                 .(triple_GEExpe)* .(ges_GEExpe)^.(pval_GEExpe)))))
    } 


  # Return
  return (expr.aov)
  }


# Get the colors and labels for treatments
get_trt_info <- function (dat, col.per = ColorsPeriod, col.trt = ColorsTrt,
                          MAIN_FACTOR, GROUP, X_LABELS, LEGEND_LABELS,
                          export_PPTX = F)
  {
  EXPE_SELECT <- levels(dat$idExperiment)
  IRR_SELECT <- levels(dat$idWatering)
  GEN_SELECT <- levels(dat$idGenotype)
  PER_SELECT <- levels(dat$idPeriod)
  
  if (MAIN_FACTOR == "idPeriod")
    {
    if (missing(X_LABELS)) { X_LABELS <- if (!export_PPTX) PER_SELECT else element_blank() }
    }
  else if (MAIN_FACTOR == "idGenotype")
    {
    if (missing(X_LABELS)) { X_LABELS <- if (!export_PPTX) GEN_SELECT else element_blank() }
    }
  else if (MAIN_FACTOR == "idWatering")
    {
    if (missing(X_LABELS)) { X_LABELS <- if (!export_PPTX) IRR_SELECT else element_blank() }
    }
  else if (MAIN_FACTOR == "idExperiment")
    {
    if (missing(X_LABELS)) { X_LABELS <- if (!export_PPTX) EXPE_SELECT else element_blank() }
    }
  else
    {
    stop ("Argument 'MAIN_FACTOR' should be \"idPeriod\", \"idGenotype\", \"idWatering\" or \"idExperiment\"")
    }   
  
  trt.col <- if (missing(GROUP)) MAIN_FACTOR else GROUP
  if (trt.col == "idPeriod")
    {
    COLOR_VALUES <- col.per$col[match(PER_SELECT, col.per$idPeriod)]
    if (missing(LEGEND_LABELS)) { LEGEND_LABELS <- PER_SELECT }
    }
  else if (trt.col == "idGenotype")
    {
    COLOR_VALUES <- col.trt$col[match(GEN_SELECT, col.trt$idGenotype)]
    if (missing(LEGEND_LABELS)) { LEGEND_LABELS <- GEN_SELECT }
    }
  else if (trt.col == "idWatering")
    {
    COLOR_VALUES <- colorRampPalette(c("dodgerblue2", "firebrick1"))(length(IRR_SELECT))
    if (missing(LEGEND_LABELS)) { LEGEND_LABELS <- IRR_SELECT }
    }
  else if (trt.col == "idExperiment")
    {
    #require(scales)
    COLOR_VALUES <- hue_pal()(length(unique(dat$idExperiment)))
    if (missing(LEGEND_LABELS)) { LEGEND_LABELS <- EXPE_SELECT }
    }
  else
    {
    stop ("Argument 'GROUP' should be \"idPeriod\", \"idGenotype\", \"idWatering\" or \"idExperiment\"")
    }
  
  return (list(COLOR_VALUES = COLOR_VALUES, X_LABELS = X_LABELS, LEGEND_LABELS = LEGEND_LABELS))
  }


# Set the names and units of the parameters

set_names_units <- function (return_mode = "units") # can be "par", "names", "units", "units_kPa", or "check" (to print the labels in the console)
  {
  EXTRACTED_PAR <- c("diel_trend",
                     "E_diel_obs", "A_rapid_op_obs", "A_rapid_op_transition_obs", "A_rapid_op_stable_obs", "A_rapid_clo_obs", "E_EON_obs",
                     "E_day_mean_obs", "E_day_max_obs", "t_day_max_obs", "sigma_day_obs", "Delta_day_obs", "Sigma_preclo_obs",
                     "E_night_mean_obs", "E_night_min_obs", "t_night_min_obs", "sigma_night_obs", "Delta_night_obs", "Sigma_preop_obs",
                     "adjR2", "RMSE", "E_EON_mod", "E_mean", "A_SQW", "A1", "phi1", "A2", "phi2", "acclim_slope",
                     "E_day_mean_mod", "E_day_max_mod", "t_day_max_mod", "sigma_day_mod", "Delta_day_mod", "Sigma_preclo_mod",
                     "E_night_mean_mod", "E_night_min_mod", "t_night_min_mod", "sigma_night_mod", "Delta_night_mod", "Sigma_preop_mod",
                     "phi1_minus_phi2", "A1_plus_A2", "daily_amp_mod", "daily_amp_obs")
  
  PAR_NAME <- c(expression(sigma["ageing"]),
                expression("E"["diel"]),
                expression("A"["rapid op"]),
                expression("A"["rapid op transition"]),
                expression("A"["rapid op stable"]),
                expression("A"["rapid clo"]),
                expression("E"["end of night"]),
                expression("E"["day"]),
                expression("E"["day max"]),
                expression("t"["day max"]),
                expression(sigma["day"]),
                expression(Delta["day"]),
                expression(Sigma["preclo"]),
                expression("E"["night"]),
                expression("E"["night min"]),
                expression("t"["night min"]),
                expression(sigma["night"]),
                expression(Delta["night"]),
                expression(Sigma["preop"]),
                expression("R"["adj"]^"2"),
                expression("RMSE"),
                expression(hat("E")["end of night"]),
                expression("E"["mean"]),
                expression("A"["SQW"]),
                expression("A"["1"]),
                expression(varphi["1"]),
                expression("A"["2"]),
                expression(varphi["2"]),
                expression(alpha),
                expression(hat("E")["day"]),
                expression(hat("E")["day max"]),
                expression(hat("t")["day max"]),
                expression(hat(sigma)["day"]),
                expression(hat(Delta)["day"]),
                expression(hat(Sigma)["preclo"]),
                expression(hat("E")["night"]),
                expression(hat("E")["night min"]),
                expression(hat("t")["night min"]),
                expression(hat(sigma)["night"]),
                expression(hat(Delta)["night"]),
                expression(hat(Sigma)["preop"]),
                expression(varphi["1"] - varphi["2"]),
                expression("A"["1"]+"A"["2"]),
                expression(hat("A")["diel"]),
                expression("A"["diel"]))
  
  PAR_UNIT <- c(expression(paste(sigma["ageing"], " (mmol m"^-2, " s"^-1, " d"^-1, ")", sep = "")),
                expression(paste("E"["diel"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("A"["rapid op"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("A"["rapid op transition"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("A"["rapid op stable"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("A"["rapid clo"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("E"["end of night"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("E"["day"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("E"["day max"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("t"["day max"], " (h since start of day)", sep = "")),
                expression(paste(sigma["day"], " (mmol m"^-2, " s"^-1, " d"^-1, ")", sep = "")),
                expression(paste(Delta["day"], " (mol m"^-2, ")", sep = "")),
                expression(paste(Sigma["preclo"], " (mol m"^-2, ")", sep = "")),
                expression(paste("E"["night"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("E"["night min"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("t"["night min"], " (h since start of night)", sep = "")),
                expression(paste(sigma["night"], " (mmol m"^-2, " s"^-1, " d"^-1, ")", sep = "")),
                expression(paste(Delta["night"], " (mol m"^-2, ")", sep = "")),
                expression(paste(Sigma["preop"], " (mol m"^-2, ")", sep = "")),
                expression("R"["adj"]^"2"),
                expression(paste("RMSE (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste(hat("E")["end of night"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("E"["mean"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("A"["SQW"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("A"["1"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste(varphi["1"], " (h)", sep = "")),
                expression(paste("A"["2"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste(varphi["2"], " (h)", sep = "")),
                expression(paste(alpha, " (mmol m"^-2, " s"^-1, " d"^-1, ")", sep = "")),
                expression(paste(hat("E")["day"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste(hat("E")["day max"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste(hat("t")["day max"], " (h since start of day)", sep = "")),
                expression(paste(hat(sigma)["day"], " (mmol m"^-2, " s"^-1, " d"^-1, ")", sep = "")),
                expression(paste(hat(Delta)["day"], " (mol m"^-2, ")", sep = "")),
                expression(paste(hat(Sigma)["preclo"], " (mol m"^-2, ")", sep = "")),
                expression(paste(hat("E")["night"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste(hat("E")["night min"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste(hat("t")["night min"], " (h since start of night)", sep = "")),
                expression(paste(hat(sigma)["night"], " (mmol m"^-2, " s"^-1, " d"^-1, ")", sep = "")),
                expression(paste(hat(Delta)["night"], " (mol m"^-2, ")", sep = "")),
                expression(paste(hat(Sigma)["preop"], " (mol m"^-2, ")", sep = "")),
                expression(paste(varphi["1"] - varphi["2"], " (h)", sep = "")),
                expression(paste("A"["1"]+"A"["2"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste(hat("A")["diel"], " (mmol m"^-2, " s"^-1, ")", sep = "")),
                expression(paste("A"["diel"], " (mmol m"^-2, " s"^-1, ")", sep = "")))
  
  PAR_UNIT_kPa <- c(expression(paste(sigma["ageing"], " (mmol m"^-2, " s"^-1, " kPa"^-1, " d"^-1, ")", sep = "")),
                    expression(paste("E"["diel"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("A"["rapid op"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("A"["rapid op transition"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("A"["rapid op stable"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("A"["rapid clo"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("E"["end of night"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("E"["day"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("E"["day max"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("t"["day max"], " (h since start of day)", sep = "")),
                    expression(paste(sigma["day"], " (mmol m"^-2, " s"^-1, " kPa"^-1, " d"^-1, ")", sep = "")),
                    expression(paste(Delta["day"], " (mol m"^-2, " kPa"^-1, ")", sep = "")),
                    expression(paste(Sigma["preclo"], " (mol m"^-2, " kPa"^-1, ")", sep = "")),
                    expression(paste("E"["night"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("E"["night min"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("t"["night min"], " (h since start of night)", sep = "")),
                    expression(paste(sigma["night"], " (mmol m"^-2, " s"^-1, " kPa"^-1, " d"^-1, ")", sep = "")),
                    expression(paste(Delta["night"], " (mol m"^-2, " kPa"^-1, ")", sep = "")),
                    expression(paste(Sigma["preop"], " (mol m"^-2, " kPa"^-1, ")", sep = "")),
                    expression("R"["adj"]^"2"),
                    expression(paste("RMSE (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste(hat("E")["end of night"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("E"["mean"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("A"["SQW"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("A"["1"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste(varphi["1"], " (h)", sep = "")),
                    expression(paste("A"["2"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste(varphi["2"], " (h)", sep = "")),
                    expression(paste(alpha, " (mmol m"^-2, " s"^-1, " kPa"^-1, " d"^-1, ")", sep = "")),
                    expression(paste(hat("E")["day"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste(hat("E")["day max"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste(hat("t")["day max"], " (h since start of day)", sep = "")),
                    expression(paste(hat(sigma)["day"], " (mmol m"^-2, " s"^-1, " kPa"^-1, " d"^-1, ")", sep = "")),
                    expression(paste(hat(Delta)["day"], " (mol m"^-2, " kPa"^-1, ")", sep = "")),
                    expression(paste(hat(Sigma)["preclo"], " (mol m"^-2, " kPa"^-1, ")", sep = "")),
                    expression(paste(hat("E")["night"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste(hat("E")["night min"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste(hat("t")["night min"], " (h since start of night)", sep = "")),
                    expression(paste(hat(sigma)["night"], " (mmol m"^-2, " s"^-1, " kPa"^-1, " d"^-1, ")", sep = "")),
                    expression(paste(hat(Delta)["night"], " (mol m"^-2, " kPa"^-1, ")", sep = "")),
                    expression(paste(hat(Sigma)["preop"], " (mol m"^-2, " kPa"^-1, ")", sep = "")),
                    expression(paste(varphi["1"] - varphi["2"], " (h)", sep = "")),
                    expression(paste("A"["1"]+"A"["2"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste(hat("A")["diel"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")),
                    expression(paste("A"["diel"], " (mmol m"^-2, " s"^-1, " kPa"^-1, ")", sep = "")))
  
  if (length(EXTRACTED_PAR) != length(PAR_NAME) |
      length(EXTRACTED_PAR) != length(PAR_UNIT) | 
      length(EXTRACTED_PAR) != length(PAR_UNIT_kPa))
    {
    stop ("Inconsistent vector lengths, check the names and units in the set_names_units() function of the PhenoLeaks_graphics.R script")
    }
  
  if (return_mode == "check")
    {
    for (i in 1:length(EXTRACTED_PAR))
      {
      print(EXTRACTED_PAR[i])
      print(PAR_NAME[i])
      print(PAR_UNIT[i])
      print(PAR_UNIT_kPa[i])
      print("                  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                  ")
      }
    }
  else if (return_mode == "par")
    {
    return(EXTRACTED_PAR)
    }
  else if (return_mode == "names")
    {
    return(PAR_NAME)
    }
  else if (return_mode == "units")
    {
    return(PAR_UNIT)
    }
  else if (return_mode == "units_kPa")
    {
    return(PAR_UNIT_kPa)
    }
  else
    {
    stop ("The argument 'return_mode' should be \"par\", \"names\", \"units\", \"units_kPa\", or \"check\".")  
    }
  }

# One-way ANOVA and jitter plots
# for statistical comparisons of genotypes ("idGenotype"), aerial environments ("idPeriod"),
# irrigation regimes ("idWatering") or experiments ("idExperiment").
# One of these should be input as 'MAIN_FACTOR' (the variable for the main category on the x-axis).
# "idPeriod" can be treated as a within-subject factor (BUT NOT RECOMMENDED IF MISSING VALUES), i.e. repeated measures on the same pot,
# while the others are always between-subjects factors, i.e. one pot belongs to only one independent category.
anova1_jitter <- function (dat, # the dataframe
                           Y, # the name of the dependent variable (character)
                           use_within = F, # set to TRUE for treating "idPeriod" as a within-subject factor. This will increase the power of the statistical analysis, but remove ALL data (of a given variable) of one pot if only one period is missing
                           MAIN_FACTOR = "idPeriod", # "idPeriod", "idGenotype", "idWatering" or "idExperiment"
                           export_CSV = F, # set to TRUE to export the results of the ANOVA and the pairwise t-tests
                           file, # the path for the CSV file, if required
                           export_PPTX = F, # set to TRUE to export the graph as PPTX instead of PDF
                           use_VPD = F,  # set to TRUE when transpiration is normalized by VPD
                           Y_LABEL = if (export_PPTX) "" else (if (is.na(match(Y, set_names_units("par")))) Y else set_names_units(if (!use_VPD) "units" else "units_kPa")[match(Y, set_names_units("par"))]), # may be useful to omit the y-axis label (while keeping a space for it) when expressions are not supported by not supported by 'graph2ppt()', e.g. with superscripts, see https://github.com/tomwenseleers/export/issues/49
                           ...) # further arguments to be passed to 'get_trt_info()' i.e.:
                                  # X_LABELS,
                                  # LEGEND_LABELS,
                                  # col.per = ColorsPeriod,
                                  # col.trt = ColorsTrt
  {
  #require(ggplot2)
  #require(ggpubr)
  #require(rstatix)
  
  # Display progress
  print(paste("Current variable: '", Y, "' (", which(colnames(dat) == Y)-5, "/", ncol(dat)-5, ")", sep = ""))
  
  # One-way ANOVA
  if (MAIN_FACTOR == "idPeriod" & use_within)
    {
    res.aov <- anova_test(data = dat,
                          dv = all_of(Y),
                          within = idPeriod,
                          wid = idPot)
    }
  else
    {
    res.aov <- anova_test(data = dat,
                          dv = all_of(Y),
                          between = all_of(MAIN_FACTOR))
    }
  
  table.aov <- get_anova_table(res.aov)
  expr.aov <- get_eta2g(table.aov, export_PPTX)

  # CSV export
  if (export_CSV & !Y %in% c("adjR2", "RMSE"))
    {
    write.table(paste(rep("~", nchar(Y)), sep = "", collapse = ""), file, row.names = F, col.names = F, sep = ",", append = T)
    write.table(Y, file, row.names = F, col.names = F, sep = ",", append = T)
    write.table(paste(rep("~", nchar(Y)), sep = "", collapse = ""), file, row.names = F, col.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    
    write.table("ANOVA table", file, row.names = F, col.names = F, sep = ",", append = T)
    write.table(table.aov, file, row.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    
    # Pairwise t-tests
    write.table(paste(MAIN_FACTOR, "pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
    if (table.aov[1, "p"] >= 0.05) write.table("WARNING - The main effect is not significant", file, row.names = F, col.names = F, sep = ",", append = T)
    pwc <- dat %>%
      pairwise_t_test(., reformulate(MAIN_FACTOR, Y), p.adjust.method = "bonferroni")
    write.table(pwc, file, row.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    }
    
  # Pairwise t-tests
  pwc <- dat %>%
    pairwise_t_test(., reformulate(MAIN_FACTOR, Y), p.adjust.method = "bonferroni")
  pwc <- pwc %>% add_xy_position(x = MAIN_FACTOR, step.increase = 0.05)
  
  # Boxplot
  c(COLOR_VALUES, X_LABELS, LEGEND_LABELS) := get_trt_info(dat, MAIN_FACTOR = MAIN_FACTOR, export_PPTX = export_PPTX, ...)

  bxp <- ggplot(dat, aes(y = get(Y), x = get(MAIN_FACTOR), color = get(MAIN_FACTOR))) +
    geom_boxplot(outlier.shape = NA, width = 0.7, position = "dodge", lwd = 0.3) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.7), size = 0.25) +
    theme_bw() +
    labs(y = Y_LABEL, x = NULL) +
    scale_color_manual(values = COLOR_VALUES, labels = LEGEND_LABELS) +
    scale_x_discrete(labels = X_LABELS)
  
  bxp +
    stat_pvalue_manual(pwc, tip.length = 0.01, hide.ns = T, label.size = 2.5) +
    labs(title = expr.aov) +
    theme(plot.title.position = "plot",
          legend.title = element_blank(),
          legend.text.align = 0,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    if (!export_PPTX)
      theme(axis.title.y = element_text(size = 10),
            plot.title = element_text(hjust = 0, vjust = 0, size = 7))
    else
      theme(axis.title.y = element_text(size = 5),
            axis.text.y = element_text(size = 5.5, colour = "black"),
            plot.title = element_text(hjust = 0.7, vjust = -4, size = 5.5, margin = margin(-5, 0, 5, 0), colour = "black"),
            panel.border = element_rect(size = 0.35, colour = "black"),
            axis.ticks = element_line(size = 0.35, colour = "black"),
            plot.margin = unit(c(0,1,0,11), "points"),
            legend.margin = margin(t = 0, b = 0, unit = "line"),
            legend.text = element_text(margin = margin(r = 10, unit = "pt")))
  }


# Two-way (mixed) ANOVA and jitter plots
# for statistical comparisons of genotypes ("idGenotype"), aerial environments ("idPeriod"),
# irrigation regimes ("idWatering") or experiments ("idExperiment").
# One of these should be input as 'MAIN_FACTOR' (the variable for the main category on the x-axis),
# and another as 'GROUP' (the variable for the pairwise t-tests).
# "idPeriod" can be treated as a within-subject factor (BUT NOT RECOMMENDED IF MISSING VALUES), i.e. repeated measures on the same pot,
# while the others are always between-subjects factors, i.e. one pot belongs to only one independent category.
anova2_jitter <- function (dat, # the dataframe
                           Y, # the name of the dependent variable (character)
                           use_within = F, # set to TRUE for treating "idPeriod" as a within-subject factor. This will increase the power of the statistical analysis, but remove ALL data (of a given variable) of one pot if only one period is missing
                           MAIN_FACTOR = "idPeriod",
                           GROUP = "idGenotype",
                           export_CSV = F, # set to TRUE to export the results of the ANOVA and the pairwise t-tests
                           file, # the path for the CSV file, if required
                           export_PPTX = F, # set to TRUE to export the graph as PPTX instead of PDF
                           use_VPD = F,  # set to TRUE when transpiration is normalized by VPD
                           Y_LABEL = if (export_PPTX) "" else (if (is.na(match(Y, set_names_units("par")))) Y else set_names_units(if (!use_VPD) "units" else "units_kPa")[match(Y, set_names_units("par"))]), # may be useful to omit the y-axis label (while keeping a space for it) when expressions are not supported by not supported by 'graph2ppt()', e.g. with superscripts, see https://github.com/tomwenseleers/export/issues/49
                           ...) # further arguments to be passed to 'get_trt_info()' i.e.:
                                  # X_LABELS,                                
                                  # LEGEND_LABELS,
                                  # col.per = ColorsPeriod,
                                  # col.trt = ColorsTrt
  {
  #require(ggplot2)
  #require(ggpubr)
  #require(rstatix)
  
  # Display progress
  print(paste("Current variable: '", Y, "' (", which(colnames(dat) == Y)-5, "/", ncol(dat)-5, ")", sep = ""))
  
  # Two-way (mixed) ANOVA
  if ("idPeriod" %in% c(MAIN_FACTOR, GROUP) & use_within)
    {
    # Type II if balanced design, type III if unbalanced (if one missing value at a given period for one pot, all values of the other periods are removed for this pot)
    between_var <- c(MAIN_FACTOR, GROUP)[!c(MAIN_FACTOR, GROUP) %in% "idPeriod"]
    res.aov <- anova_test(data = dat,
                          dv = all_of(Y),
                          between = all_of(between_var),
                          within = idPeriod,
                          wid = idPot)
    }
  else
    {
    # Type II if balanced design, type III if unbalanced (if one missing value at a given period for one pot, the values of the other periods are kept)
    res.aov <- anova_test(data = dat,
                          dv = all_of(Y),
                          between = c(all_of(MAIN_FACTOR), all_of(GROUP)))
    # Always type II
    #res.aov <- dat %>% anova_test(., reformulate(paste(GROUP, MAIN_FACTOR, sep = "*"), Y))
    }
  
  table.aov <- get_anova_table(res.aov)
  expr.aov <- get_eta2g(table.aov, export_PPTX)
  
  # CSV export
  if (export_CSV & !Y %in% c("adjR2", "RMSE"))
    {
    write.table(paste(rep("~", nchar(Y)), sep = "", collapse = ""), file, row.names = F, col.names = F, sep = ",", append = T)
    write.table(Y, file, row.names = F, col.names = F, sep = ",", append = T)
    write.table(paste(rep("~", nchar(Y)), sep = "", collapse = ""), file, row.names = F, col.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)

    write.table("ANOVA table", file, row.names = F, col.names = F, sep = ",", append = T)
    write.table(table.aov, file, row.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    
    # Pairwise t-tests
    for (effect in 1:2)
      {
      write.table(paste(table.aov$Effect[effect], "pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
      if (table.aov[table.aov == table.aov$Effect[effect], "p"] >= 0.05) write.table("WARNING - The main effect is not significant", file, row.names = F, col.names = F, sep = ",", append = T)
      if (table.aov[3, "p"] < 0.05) write.table("WARNING - There is a significant interaction term", file, row.names = F, col.names = F, sep = ",", append = T)
      pwc <- dat %>%
        pairwise_t_test(., reformulate(table.aov$Effect[effect], Y), p.adjust.method = "bonferroni")
      write.table(pwc, file, row.names = F, sep = ",", append = T)
      write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
      }
 
    write.table(paste(GROUP, "conditional pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
    if (table.aov[3, "p"] >= 0.05) write.table("WARNING - The interaction term is not significant", file, row.names = F, col.names = F, sep = ",", append = T)
    pwc <- dat %>%
      group_by(get(MAIN_FACTOR)) %>%
      pairwise_t_test(., reformulate(GROUP, Y), p.adjust.method = "bonferroni")
    colnames(pwc)[1] <- MAIN_FACTOR
    write.table(pwc, file, row.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    
    write.table(paste(MAIN_FACTOR, "conditional pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
    if (table.aov[3, "p"] >= 0.05) write.table("WARNING - The interaction term is not significant", file, row.names = F, col.names = F, sep = ",", append = T)
    pwc <- dat %>%
      group_by(get(GROUP)) %>%
      pairwise_t_test(., reformulate(MAIN_FACTOR, Y), p.adjust.method = "bonferroni")
    colnames(pwc)[1] <- GROUP
    write.table(pwc, file, row.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    
    write.table(paste(rep("_", 150), sep = "", collapse = ""), file, row.names = F, col.names = F, sep = ",", append = T)  
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    }
  
  # Pairwise t-tests
  pwc <- dat %>%
    group_by(get(MAIN_FACTOR)) %>%
    pairwise_t_test(., reformulate(GROUP, Y), p.adjust.method = "bonferroni")
  colnames(pwc)[1] <- MAIN_FACTOR
  pwc <- pwc %>% add_xy_position(x = MAIN_FACTOR, step.increase = 0.05)
  
  # Boxplot
  c(COLOR_VALUES, X_LABELS, LEGEND_LABELS) := get_trt_info(dat, MAIN_FACTOR = MAIN_FACTOR, GROUP = GROUP, export_PPTX = export_PPTX, ...)

  bxp <- ggplot(dat, aes(y = get(Y), x = get(MAIN_FACTOR), color = get(GROUP))) +
    geom_boxplot(outlier.shape = NA, width = 0.7, position = "dodge", lwd = 0.3) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.7), size = 0.25) +
    theme_bw() +
    labs(y = Y_LABEL, x = NULL) +
    scale_color_manual(values = COLOR_VALUES, labels = LEGEND_LABELS) +
    scale_x_discrete(labels = X_LABELS)
  
  bxp +
    stat_pvalue_manual(pwc, tip.length = 0.01, hide.ns = T, label.size = 2.5) +
    labs(title = expr.aov) +
    theme(plot.title.position = "plot",
          legend.title = element_blank(),
          legend.text.align = 0,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    if (!export_PPTX)
      theme(axis.title.y = element_text(size = 10),
            plot.title = element_text(hjust = 0, vjust = 0, size = 7))
    else
      theme(axis.title.y = element_text(size = 5),
            axis.text.y = element_text(size = 5.5, colour = "black"),
            plot.title = element_text(hjust = 1, vjust = -4, size = 5.5, margin = margin(-5, 0, 5, 0), colour = "black"),
            panel.border = element_rect(size = 0.35, colour = "black"),
            axis.ticks = element_line(size = 0.35, colour = "black"),
            plot.margin = unit(c(0,1,0,11), "points"),
            legend.margin = margin(t = 0, b = 0, unit = "line"),
            legend.text = element_text(margin = margin(r = 10, unit = "pt")))
  
  }


# Three-way mixed ANOVA and jitter plots
# for statistical comparisons of genotypes ("idGenotype"), aerial environments ("idPeriod"),
# irrigation regimes ("idWatering") or experiments ("idExperiment").
# One of these should be input as 'MAIN_FACTOR' (the variable for the main category on the x-axis),
# another as 'GROUP' (the variable for the pairwise t-tests),
# and another as 'FACET' (one series of boxplots per facet level).
# "idPeriod" can be treated as a within-subject factor (BUT NOT RECOMMENDED IF MISSING VALUES), i.e. repeated measures on the same pot,
# while the others are always between-subjects factors, i.e. one pot belongs to only one independent category.
anova3_jitter <- function (dat, # the dataframe
                           Y, # the name of the dependent variable (character)
                           use_within = F, # set to TRUE for treating "idPeriod" as a within-subject factor. This will increase the power of the statistical analysis, but remove ALL data (of a given variable) of one pot if only one period is missing
                           MAIN_FACTOR = "idPeriod",
                           GROUP = "idGenotype",
                           FACET = "idExperiment", # or "idWatering"
                           export_CSV = F, # set to TRUE to export the results of the ANOVA and the pairwise t-tests
                           file, # the path for the CSV file, if required
                           export_PPTX = F, # set to TRUE to export the graph as PPTX instead of PDF
                           use_VPD = F,  # set to TRUE when transpiration is normalized by VPD
                           Y_LABEL = if (export_PPTX) "" else (if (is.na(match(Y, set_names_units("par")))) Y else set_names_units(if (!use_VPD) "units" else "units_kPa")[match(Y, set_names_units("par"))]), # may be useful to omit the y-axis label (while keeping a space for it) when expressions are not supported by not supported by 'graph2ppt()', e.g. with superscripts, see https://github.com/tomwenseleers/export/issues/49
                           ...) # further arguments to be passed to 'get_trt_info()' i.e.:
                                  # X_LABELS,                                
                                  # LEGEND_LABELS,
                                  # col.per = ColorsPeriod,
                                  # col.trt = ColorsTrt
  {
  #require(ggplot2)
  #require(ggpubr)
  #require(rstatix)

  # Display progress
  print(paste("Current variable: '", Y, "' (", which(colnames(dat) == Y)-5, "/", ncol(dat)-5, ")", sep = ""))
  
  # Three-way (mixed) ANOVA
  if ("idPeriod" %in% c(MAIN_FACTOR, GROUP, FACET) & use_within)
    {
    between_var <- c(MAIN_FACTOR, GROUP, FACET)[!c(MAIN_FACTOR, GROUP, FACET) %in% "idPeriod"]
    res.aov <- anova_test(data = dat,
                          dv = all_of(Y),
                          between = c(all_of(between_var[1]), all_of(between_var[2])),
                          within = idPeriod,
                          wid = idPot)
    }
  else
    {
    res.aov <- anova_test(data = dat,
                          dv = all_of(Y),
                          between = c(all_of(MAIN_FACTOR), all_of(GROUP), all_of(FACET)))
    }

  table.aov <- get_anova_table(res.aov)
  expr.aov <- get_eta2g(table.aov, export_PPTX)
  
  # CSV export
  if (export_CSV & !Y %in% c("adjR2", "RMSE"))
    {
    write.table(paste(rep("~", nchar(Y)), sep = "", collapse = ""), file, row.names = F, col.names = F, sep = ",", append = T)
    write.table(Y, file, row.names = F, col.names = F, sep = ",", append = T)
    write.table(paste(rep("~", nchar(Y)), sep = "", collapse = ""), file, row.names = F, col.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    
    write.table("ANOVA table", file, row.names = F, col.names = F, sep = ",", append = T)
    write.table(table.aov, file, row.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    
    # Pairwise t-tests
    for (effect in 1:3)
      {
      write.table(paste(table.aov$Effect[effect], "pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
      if (table.aov[table.aov == table.aov$Effect[effect], "p"] >= 0.05) write.table("WARNING - The main effect is not significant", file, row.names = F, col.names = F, sep = ",", append = T)
      if (any(table.aov[4:7, "p"] < 0.05)) write.table("WARNING - There is a significant interaction term", file, row.names = F, col.names = F, sep = ",", append = T)
      pwc <- dat %>%
        pairwise_t_test(., reformulate(table.aov$Effect[effect], Y), p.adjust.method = "bonferroni")
      write.table(pwc, file, row.names = F, sep = ",", append = T)
      write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
      }
    
    write.table(paste(GROUP, "conditional pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
    if (table.aov[7, "p"] >= 0.05) write.table("WARNING - The triple interaction term is not significant", file, row.names = F, col.names = F, sep = ",", append = T)
    pwc <- dat %>%
      group_by(get(FACET), get(MAIN_FACTOR)) %>%
      pairwise_t_test(., reformulate(GROUP, Y), p.adjust.method = "bonferroni")
    colnames(pwc)[1] <- FACET
    colnames(pwc)[2] <- MAIN_FACTOR
    write.table(pwc, file, row.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    
    write.table(paste(MAIN_FACTOR, "conditional pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
    if (table.aov[7, "p"] >= 0.05) write.table("WARNING - The triple interaction term is not significant", file, row.names = F, col.names = F, sep = ",", append = T)
    pwc <- dat %>%
      group_by(get(FACET), get(GROUP)) %>%
      pairwise_t_test(., reformulate(MAIN_FACTOR, Y), p.adjust.method = "bonferroni")
    colnames(pwc)[1] <- FACET
    colnames(pwc)[2] <- GROUP
    write.table(pwc, file, row.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    
    write.table(paste(FACET, "conditional pairwise comparisons"), file, row.names = F, col.names = F, sep = ",", append = T)
    if (table.aov[7, "p"] >= 0.05) write.table("WARNING - The triple interaction term is not significant", file, row.names = F, col.names = F, sep = ",", append = T)
    pwc <- dat %>%
      group_by(get(GROUP), get(MAIN_FACTOR)) %>%
      pairwise_t_test(., reformulate(FACET, Y), p.adjust.method = "bonferroni")
    colnames(pwc)[1] <- GROUP
    colnames(pwc)[2] <- MAIN_FACTOR
    write.table(pwc, file, row.names = F, sep = ",", append = T)
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    
    write.table(paste(rep("_", 150), sep = "", collapse = ""), file, row.names = F, col.names = F, sep = ",", append = T)  
    write.table("", file, row.names = F, col.names = F, sep = ",", append = T)
    }

  # Pairwise t-tests
  pwc <- dat %>%
    group_by(get(FACET), get(MAIN_FACTOR)) %>%
    pairwise_t_test(., reformulate(GROUP, Y), p.adjust.method = "bonferroni")
  colnames(pwc)[1] <- FACET
  colnames(pwc)[2] <- MAIN_FACTOR
  pwc <- pwc %>% add_xy_position(x = MAIN_FACTOR, step.increase = 0.05)
  
  # Boxplot
  c(COLOR_VALUES, X_LABELS, LEGEND_LABELS) := get_trt_info(dat, MAIN_FACTOR = MAIN_FACTOR, GROUP = GROUP, export_PPTX = export_PPTX, ...)

  bxp <- ggplot(dat, aes(y = get(Y), x = get(MAIN_FACTOR), color = get(GROUP))) +
    facet_wrap(~ get(FACET)) +
    geom_boxplot(outlier.shape = NA, width = 0.7, position = "dodge", lwd = 0.3) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.7), size = 0.25) +
    theme_bw() +
    labs(y = Y_LABEL, x = NULL) +
    scale_color_manual(values = COLOR_VALUES, labels = LEGEND_LABELS) +
    scale_x_discrete(labels = X_LABELS)

  bxp +
    stat_pvalue_manual(pwc, tip.length = 0.01, hide.ns = T, label.size = 2.5) +
    labs(title = expr.aov) +
    theme(plot.title.position = "plot",
          legend.title = element_blank(),
          legend.text.align = 0,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    if (!export_PPTX)
      theme(axis.title.y = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 7))
    else
      theme(axis.title.y = element_text(size = 5),
            axis.text.y = element_text(size = 5.5, colour = "black"),
            plot.title = element_text(hjust = 1, vjust = -2, size = 5.5, margin = margin(-5, 0, 5, 0), colour = "black"),
            panel.border = element_rect(size = 0.35, colour = "black"),
            axis.ticks = element_line(size = 0.35, colour = "black"),
            plot.margin = unit(c(0,1,0,11), "points"),
            legend.margin = margin(t = 0, b = 0, unit = "line"),
            legend.text = element_text(margin = margin(r = 10, unit = "pt")))

  }


# Perform ANOVAs and jitter plots as batch
# by calling either 'anova2_jitter()' or 'anova3_jitter()'
anova_batch <- function (results_fit, # the dataframe containing the results of the fit
                         anova_jitter_function = anova2_jitter, # either anova2_jitter or anova3_jitter, i.e. the name of the function for ANOVAs and jitter plots
                         VAR_SELECT = NULL, # the name of the variable(s) to be considered (if NULL, all dependent variables will be considered)
                         EXPE_SELECT = unique(results_fit$idExperiment), # the experiment(s) to be considered (use several only with anova3_jitter and FACET = "idExperiment", otherwise data of several experiments will be pooled)
                         PER_SELECT = unique(results_fit$idPeriod), # the period(s) to be considered
                         GEN_SELECT = unique(results_fit$idGenotype), # the genotype(s) to be considered
                         IRR_SELECT = "WW", # the watering regime(s) to be considered (use several only with anova3_jitter and FACET = "idWatering", otherwise data of several watering regimes will be pooled)
                         filter_phi_for_low_A = T, # assign NA to phase values if the magnitude is very small (making it hard to estimate the phase)
                         ...) # further arguments to be passed to 'anova_jitter_function' i.e.:
                                # use_within = T, # to treat "idPeriod" as a within-subject factor
                                # MAIN_FACTOR = "idPeriod",
                                # GROUP = "idGenotype",
                                # FACET = "idExperiment" or "idWatering" - only with anova3_jitter
                                # export_CSV = F, # set to TRUE to export the results of the ANOVA and the pairwise t-tests
                                # file, # the path for the CSV file, if required
                                # export_PPTX = F # set to TRUE to export the graph as PPTX instead of PDF
                                # use_VPD = F # set to TRUE when transpiration is normalized by VPD
                                # Y_LABEL = # label of the y-axis (default depends on 'export_PPTX' and 'use_VPD')
                                # X_LABELS, # the legend for the x-labels (MAIN_FACTOR)
                                # LEGEND_LABELS, # the legend for the GROUP (or MAIN_FACTOR again if anova1_jitter)
                                # col.per = ColorsPeriod,
                                # col.trt = ColorsTrt
  {
  dat_select <- results_fit[results_fit$idExperiment %in% EXPE_SELECT &
                              results_fit$idWatering %in% IRR_SELECT &
                              results_fit$idGenotype %in% GEN_SELECT &
                              results_fit$idPeriod %in% PER_SELECT, ]
  if (filter_phi_for_low_A)
    {
    # filter_A1overA2 <- which(dat_select$A1/dat_select$A2 < 1 & !is.na(dat_select$A1/dat_select$A2))
    # if (length(filter_A1overA2) > 0)
    #   {
    #   dat_select[filter_A1overA2, c("A1", "A2", "A1_plus_A2", "phi1", "phi2", "phi1_minus_phi2")] <- NA
    #   print(paste("Doubtful fit: ", length(filter_A1overA2), " 'A1', 'A2', phi1' and 'phi2' removed because 'A1' lower than 'A2':", sep = ""))
    #   for (i in 1:length(filter_A1overA2)) print(dat_select[filter_A1overA2[i], which(names(dat_select) %in% c("idExperiment", "idGenotype", "idWatering", "idPot", "idPeriod"))])
    #   }
    filter_A1 <- which(dat_select$A1 < 0.05 & !is.na(dat_select$A1))
    if (length(filter_A1) > 0)
      {
      dat_select[filter_A1, c("phi1", "phi2", "phi1_minus_phi2")] <- NA
      print(paste("Amplitude too low for reliable phase estimate: ", length(filter_A1), " 'phi1' (and 'phi2') removed because of low 'A1':", sep = ""))
      for (i in 1:length(filter_A1)) print(dat_select[filter_A1[i], which(names(dat_select) %in% c("idExperiment", "idGenotype", "idWatering", "idPot", "idPeriod"))])
      }
    filter_A2 <- which(dat_select$A2 < 0.02 & !is.na(dat_select$A2))
    if (length(filter_A2) > 0)
      {
      dat_select[filter_A2, c("phi2", "phi1_minus_phi2")] <- NA
      print(paste("Amplitude too low for reliable phase estimate: ", length(filter_A2), " 'phi2' removed because of low 'A2':", sep = ""))
      for (i in 1:length(filter_A2)) print(dat_select[filter_A2[i], which(names(dat_select) %in% c("idExperiment", "idGenotype", "idWatering", "idPot", "idPeriod"))])
      }
    }
  dat_select$idGenotype <- factor(as.character(dat_select$idGenotype), levels = GEN_SELECT, ordered = T)
  dat_select$idPeriod <- factor(as.character(dat_select$idPeriod), levels = PER_SELECT, ordered = T)
  dat_select$idExperiment <- factor(as.character(dat_select$idExperiment), levels = EXPE_SELECT, ordered = T)
  dat_select$idWatering <- factor(as.character(dat_select$idWatering), levels = IRR_SELECT, ordered = T)
  
  if (is.null(VAR_SELECT)) VAR_SELECT <- colnames(dat_select)[!colnames(dat_select) %in% c("idExperiment", "idGenotype", "idWatering", "SWC_drop", "idPot", "idPeriod", "diel_trend", "diel_trend_baseline", "t1_fit", "t2_fit", "E_EON_before_obs", "E_EON_before_mod", "A_rapid_op_transition_obs", "A_rapid_op_stable_obs")]
  dat_select <- dat_select[, colnames(dat_select) %in% c("idExperiment", "idGenotype", "idWatering", "idPot", "idPeriod", VAR_SELECT)]
  
  myplots <- lapply(VAR_SELECT, anova_jitter_function, dat = dat_select, ...)
  
  return (list(myplots = myplots, VAR_SELECT = VAR_SELECT))
  }

