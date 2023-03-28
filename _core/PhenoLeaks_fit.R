################################################################################
#                                                                              #
#                               PhenoLeaks - FIT                               #
#                                                                              #
#           Core functions to fit a sine model to transpiration data           #
#                        and extract dynamic parameters                        #
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
for (pkg in c("pracma", "fBasics", "gsignal"))
  {
  if (!pkg %in% installed.packages()[, "Package"]) { install.packages(pkg) }
  #update.packages(pkg)
  library(pkg, character.only = T)
  }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                                (2)  Functions                                #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#------------------------------------------------------------------------------#
#                                     Model                                    #
#------------------------------------------------------------------------------#


# Function to model the diel transpiration
# using a a square + sine (first and second harmonics) wave
E_diel <- function (TIME, # time (days)
                    E_mean, # mean transpiration over a 24-h period (mmol/m2/s)
                    A_SQW, # semi-amplitude of the square wave (mmol/m2/s)
                    A1, phi1, A2, phi2, # semi-amplitude (mmol/m2/s) and phase (days) of the first and second harmonics
                    use_shifted_time = T, # set to TRUE if time has been shifted to match 0 when light switches on at the reference cycle (typically done at Step#02)
                    light_ON = if (use_shifted_time) 0 else Time_ON0 - floor(Time_ON0)) # time of the diel cycle when the light switches on
  {
  E_mod <- E_mean + A_SQW*sign(sin((TIME-light_ON-SQW_shift/60/24)*2*pi)) + A1*sin((TIME-light_ON+phi1)*2*pi) + A2*sin((TIME-light_ON+phi2)*4*pi)
  return (E_mod)
  }


# Function to convert the parameters of a sine wave
# from a*sin(x) + b*cos(x) to A*sin(x+phi)
transform_sine <- function (a, b) # coefficients from the linear regression
  {
  A <- sqrt(a^2+b^2)
  phi <- ifelse(a > 0, atan(b/a), atan(b/a) + pi)
  return (list(A = as.numeric(A), phi = as.numeric(phi)))
  }


# Function to perform a linear regression fitting the transpiration data
# to a square + sine (first and second harmonics) wave
lm_SQW_SIN2 <- function (xy, # the dataframe with the time (x) and transpiration (y) values
                         use_shifted_time = T, # set to TRUE if time has been shifted to match 0 when light switches on at the reference cycle (typically done at Step#02)
                         light_ON = if (use_shifted_time) 0 else Time_ON0 - floor(Time_ON0)) # time of the diel cycle when the light switches on
  {
  mod_sq_sin2 <- lm(y ~ sign(sin((x-light_ON-SQW_shift/60/24)*2*pi)) + sin((x-light_ON)*2*pi) + cos((x-light_ON)*2*pi) + sin((x-light_ON)*4*pi) + cos((x-light_ON)*4*pi),
                    data = xy, na.action = na.exclude)
  return(mod_sq_sin2)
  }


# Function to produce a linear pulse (between -1 to 0)
# during a defined acclimation period ('acclim_range')
acclim <- function (x,
                    acclim_range = c(0, 1), # lower and upper time limits of the acclimation period (automatically shifted by 'SQW_shift')
                    use_shifted_time = T, # set to TRUE if time has been shifted to match 0 when light switches on at the reference cycle (typically done at Step#02)
                    light_ON = if (use_shifted_time) 0 else Time_ON0 - floor(Time_ON0)) # time of the diel cycle when the light switches on
  {
  #require("fBasics") # for 'Boxcar()'
  #require("gsignal") # for 'sawtooth()'
  if (acclim_range[2] <= acclim_range[1]) { stop ("Invalid 'acclim_range': upper limit should be higher than lower limit") }
  range_pulse <- acclim_range[2] - acclim_range[1]
  boxcar_pulse <- Boxcar(x - (acclim_range[1] + range_pulse/2) - light_ON - SQW_shift/60/24, a = -range_pulse/2)
  sawtooth_periodic <- -(0.5+0.5*sawtooth(-2*pi*(x - light_ON - SQW_shift/60/24)/range_pulse))
  linear_pulse <- sawtooth_periodic * boxcar_pulse * range_pulse
  return (linear_pulse)
  }


# Function to model the diel transpiration
# using a a square + sine (first and second harmonics) wave
# and a linear pulse during a defined acclimation period
E_diel_acclim <- function (TIME, # time (days)
                           E_mean, # mean transpiration over a 24-h period (mmol/m2/s)
                           A_SQW, # semi-amplitude of the square wave (mmol/m2/s)
                           A1, phi1, A2, phi2, # semi-amplitude (mmol/m2/s) and phase (days) of the first and second harmonics
                           acclim_slope, # slope of the linear pulse during acclimation (mmol/m2/s/d)
                           acclim_range = c(0, 1), # lower and upper time limits of the acclimation period (automatically shifted by 'SQW_shift')
                           use_shifted_time = T, # set to TRUE if time has been shifted to match 0 when light switches on at the reference cycle (typically done at Step#02)
                           light_ON = if (use_shifted_time) 0 else Time_ON0 - floor(Time_ON0)) # time of the diel cycle when the light switches on
  {
  E_mod_acclim <- E_diel(TIME, E_mean, A_SQW, A1, phi1, A2, phi2, use_shifted_time, light_ON) + acclim_slope * acclim(TIME, acclim_range, use_shifted_time, light_ON)
  return (E_mod_acclim)
  }


# Function to perform a linear regression fitting the transpiration data
# to a square + sine (first and second harmonics) wave
# and a linear pulse during a defined acclimation period
lm_SQW_SIN2_acclim <- function (xy, # the dataframe with the time (x) and transpiration (y) values
                                acclim_range = c(0, 1), # lower and upper time limits of the acclimation period (automatically shifted by 'SQW_shift')
                                use_shifted_time = T, # set to TRUE if time has been shifted to match 0 when light switches on at the reference cycle (typically done at Step#02)
                                light_ON = if (use_shifted_time) 0 else Time_ON0 - floor(Time_ON0)) # time of the diel cycle when the light switches on
  {
  mod_sq_sin2_acclim <- lm(y ~ sign(sin((x-light_ON-SQW_shift/60/24)*2*pi)) + sin((x-light_ON)*2*pi) + cos((x-light_ON)*2*pi) + sin((x-light_ON)*4*pi) + cos((x-light_ON)*4*pi) + acclim(x, acclim_range, use_shifted_time, light_ON),
                           data = xy, na.action = na.exclude)
  
  return(mod_sq_sin2_acclim)
  }


# Function to retrieve the fitted biological parameters
# from the linear regression (standard fit or fit with linear acclimation)
retrieve_pars <- function (mod) # the fitted linear model
  {
  fit <- coef(mod)
  E_mean <- fit[1] # mean transpiration over a 24-h period (y-intercept of the linear model)
  A_SQW <- fit[2] # semi-amplitude of the square wave accounting for the rapid (< 30 min) opening (resp. closing) in the light (resp. dark)
  c(A1, phi1) := transform_sine(fit[3], fit[4]) # semi-amplitude and phase of the first harmonics
  phi1 <- phi1/2/pi # switch from radians to days
  if (phi1 <= -0.5) { phi1 <- phi1 + 1 } else if (phi1 > 0.5) { phi1 <- phi1 - 1 } # wrap the phase
  c(A2, phi2) := transform_sine(fit[5], fit[6]) # semi-amplitude and phase of the second harmonics
  phi2 <- phi2/4/pi # switch from radians to days
  if (phi2 <= -0.25) { phi2 <- phi2 + 0.5 } else if (phi2 > 0.25) { phi2 <- phi2 - 0.5 } # wrap the phase
  
  acclim_slope <- if (length(fit) == 6) NA else as.numeric(fit[7])
  
  return (list(E_mean = as.numeric(E_mean), A_SQW = as.numeric(A_SQW),
                 A1 = A1, phi1 = phi1, A2 = A2, phi2 = phi2,
                 acclim_slope = acclim_slope))
  }




# # Function to retrieve the fitted biological parameters
# # from the linear regression
# retrieve_pars <- function (mod) # the fitted linear model
#   {
#   fit <- coef(mod)
#   E_mean <- fit[1] # mean transpiration over a 24-h period (y-intercept of the linear model)
#   A_SQW <- fit[2] # semi-amplitude of the square wave accounting for the rapid (< 30 min) opening (resp. closing) in the light (resp. dark)
#   c(A1, phi1) := transform_sine(fit[3], fit[4]) # semi-amplitude and phase of the first harmonics
#   phi1 <- phi1/2/pi # switch from radians to days
#   if (phi1 <= -0.5) { phi1 <- phi1 + 1 } else if (phi1 > 0.5) { phi1 <- phi1 - 1 } # wrap the phase
#   c(A2, phi2) := transform_sine(fit[5], fit[6]) # semi-amplitude and phase of the second harmonics
#   phi2 <- phi2/4/pi # switch from radians to days
#   if (phi2 <= -0.25) { phi2 <- phi2 + 0.5 } else if (phi2 > 0.25) { phi2 <- phi2 - 0.5 } # wrap the phase
#   return (list(E_mean = as.numeric(E_mean), A_SQW = as.numeric(A_SQW),
#                A1 = A1, phi1 = phi1, A2 = A2, phi2 = phi2))
#   }



#------------------------------------------------------------------------------#
#                             Parameter extraction                             #
#------------------------------------------------------------------------------#


# Function to extract some numeric values of interest
# for transpiration based on the modelled data
E_val <- function (E_mean, A_SQW, A1, phi1, A2, phi2, # see E_diel()
                   acclim_slope = 0, acclim_range = c(0, 1),  # see E_diel_acclim() - note that with 'acclim_slope = 0', E_diel(...) == E_diel_acclim(...)
                   use_shifted_time = T, # set to TRUE if time has been shifted to match 0 when light switches on at the reference cycle (typically done at Step#02)
                   light_ON = if (use_shifted_time) 0 else Time_ON0 - floor(Time_ON0)) # time of the diel cycle when the light switches on
  {
  
  # Initial check on the acclimation slope
  if (is.na(acclim_slope)) 
    {
    # This may occur if 'lm_SQW_SIN2_acclim()' has been fitted using only "acclimated" data
    # but no "acclimating" data within 'acclim_range' (e.g. idPot "248" in idExp "C2M47").
    # In this case the fit parameters are the same as with those obtained by the standard 'lm_SQW_SIN2()',
    # and it is not possible to estimate the values of interest over the acclimation range.
    return (list(E_EON_mod = NA,
                 E_day_mean_mod = NA, E_day_max_mod = NA, t_day_max_mod = NA,
                 sigma_day_mod = NA, Delta_day_mod = NA, Sigma_preclo_mod = NA,
                 E_night_mean_mod = NA, E_night_min_mod = NA, t_night_min_mod = NA,
                 sigma_night_mod = NA, Delta_night_mod = NA, Sigma_preop_mod = NA))
    }
  
  # Initial check on the acclimation range
  range_pulse <- acclim_range[2] - acclim_range[1]
  if (acclim_range[1] - light_ON != round(acclim_range[1] - light_ON) | range_pulse != 1) { stop ("'acclim_range' should start when the light switches on, and ends 24 h later (other cases not implemented yet)") }
  acclim_range <- acclim_range - acclim_range[1] # translate back to 0 for compatibility with the rest of the function
  
  
  ## DAY ##
  
  # Numeric vector of day times, 30 min after the start of the day
  t_SOD <- Rapid_dur/60/24
  t_EOD <- Pho_Per/24
  time_day_vector <- seq(t_SOD, t_EOD, length.out = 100*2*(Pho_Per/24)/(Pho_Per/24+Sko_Per/24))
  
  # Numeric vector of daytime transpiration at time_day_vector
  E_day_vector <- E_diel_acclim(TIME = time_day_vector, E_mean, A_SQW, A1, phi1, A2, phi2, acclim_slope, acclim_range, use_shifted_time = T, light_ON = 0)
  
  # Transpiration at the start and the end of the day
  E_SOD <- E_day_vector[1] # == E_diel_acclim(TIME = t_SOD, E_mean, A_SQW, A1, phi1, A2, phi2, acclim_slope, acclim_range, use_shifted_time = T, light_ON = 0)
  E_EOD <- E_day_vector[length(E_day_vector)] # == E_diel_acclim(TIME = t_EOD, E_mean, A_SQW, A1, phi1, A2, phi2, acclim_slope, acclim_range, use_shifted_time = T, light_ON = 0)
  
  # Maximal daytime transpiration
  day_max <- which.max(E_day_vector) # index
  E_day_max <- E_day_vector[day_max] # mmol/m2/s
  
  # Time when maximal daytime transpiration occurs
  t_day_max <- time_day_vector[day_max] # fraction of day so that 0 is the start of the day and Pho_Per/24 [e.g. 0.5] the end of the day
  
  # Average daytime transpiration
  E_day_mean <- mean(E_day_vector) # mmol/m2/s
  
  # Average slope of change in transpiration throughout the daytime after 30-min illumination (sigma_day),
  # equals to the integral of the time-derivative of transpiration (so back to the E function) divided by time interval,
  # which also equals to the linear slope between the start (after the square wave) and the end of the day
  sigma_day <- (E_EOD - E_SOD) / (t_EOD - t_SOD) # mmol/m2/s/d
  # and can also be estimated by averaging the time-derivative over the time interval:
  #day_deriv <- A1*2*pi*cos((time_day_vector+phi1)*2*pi) + A2*4*pi*cos((time_day_vector+phi2)*4*pi) + acclim_slope
  #sigma_day <- mean(day_deriv)
  
  # Cumulative transpiration dynamics above the daytime trend (Delta_day, "absolute daytime dynamics"),
  # equals to the area above the trapezoid formed by E_SOD and E_EOD (i.e. the diurnal trend),
  # which is usually positive (yet some parts may be negative if the curve is convex)
  primitive_t_SOD <- -(A1/(2*pi))*cos((t_SOD+phi1)*2*pi) - (A2/(4*pi))*cos((t_SOD+phi2)*4*pi) + (E_mean+A_SQW)*t_SOD + 0.5*acclim_slope*t_SOD^2 - acclim_slope*(range_pulse+SQW_shift/60/24)*t_SOD
  primitive_t_EOD <- -(A1/(2*pi))*cos((t_EOD+phi1)*2*pi) - (A2/(4*pi))*cos((t_EOD+phi2)*4*pi) + (E_mean+A_SQW)*t_EOD + 0.5*acclim_slope*t_EOD^2 - acclim_slope*(range_pulse+SQW_shift/60/24)*t_EOD
  trapezoid_area <- (E_SOD+E_EOD)*(t_EOD-t_SOD)/2
  Delta_day <- primitive_t_EOD - primitive_t_SOD - trapezoid_area # (mmol/m2/s) * d
  # that can also be computed by numerical integration like in 'E_val_obs()':
  #Delta_day <- trapz(time_day_vector, E_day_vector) - trapz(c(Rapid_dur/60/24, Pho_Per/24), c(E_SOD, E_EOD)) # (mmol/m2/s) * d
  Delta_day <- Delta_day*24*60*60/1000 # mol/m2
  
  # # (not used anymore)
  # # "Relative daytime dynamics":
  # # algebraic sum (may be positive or negative) of
  # # the (positive) area above E_SOD (i.e. the initial baseline) and
  # # the (negative) area below E_SOD
  # rectangle_area <- E_SOD*(t_EOD-t_SOD)
  # rel_dyn <- primitive_t_EOD - primitive_t_SOD - rectangle_area # (mmol/m2/s) * d
  # rel_dyn <- rel_dyn*24*60*60/1000 # mol/m2
  
  # Cumulative preclosure once maximal transpiration has been reached (Sigma_preclo),
  # equals to the opposite of the integral of transpiration between t_day_max and t_EOD
  # (y-translated so that transpiration is 0 at t_day_max)
  primitive_t_day_max_transl <- -(A1/(2*pi))*cos((t_day_max+phi1)*2*pi) - (A2/(4*pi))*cos((t_day_max+phi2)*4*pi) + (E_mean+A_SQW-E_day_max)*t_day_max + 0.5*acclim_slope*t_day_max^2 - acclim_slope*(range_pulse+SQW_shift/60/24)*t_day_max
  primitive_t_EOD_transl <- -(A1/(2*pi))*cos((t_EOD+phi1)*2*pi) - (A2/(4*pi))*cos((t_EOD+phi2)*4*pi) + (E_mean+A_SQW-E_day_max)*t_EOD + 0.5*acclim_slope*t_EOD^2 - acclim_slope*(range_pulse+SQW_shift/60/24)*t_EOD
  Sigma_preclo <- -(primitive_t_EOD_transl - primitive_t_day_max_transl) # (mmol/m2/s) * d
  # that can also be computed by numerical integration like in 'E_val_obs()':
  #Sigma_preclo <- -trapz(time_day_vector[time_day_vector >= t_day_max], E_day_vector[time_day_vector >= t_day_max] - E_day_max) # (mmol/m2/s) * d
  Sigma_preclo <- Sigma_preclo*24*60*60/1000 # mol/m2  

  
  ## NIGHT ##
  
  # Numeric vector of night times, 30 min after the start of the night
  t_SON <- Pho_Per/24+Rapid_dur/60/24
  t_EON <- (Pho_Per+Sko_Per)/24
  time_night_vector <- seq(t_SON, t_EON, length.out = 100*2*(Sko_Per/24)/(Pho_Per/24+Sko_Per/24))
  
  # Numeric vector of nighttime transpiration at time_night_vector
  E_night_vector <- E_diel_acclim(TIME = time_night_vector, E_mean, A_SQW, A1, phi1, A2, phi2, acclim_slope, acclim_range, use_shifted_time = T, light_ON = 0)
  
  # Transpiration at the start and the end of the night
  E_SON <- E_night_vector[1] # == E_diel_acclim(TIME = t_SON, E_mean, A_SQW, A1, phi1, A2, phi2, acclim_slope, acclim_range, use_shifted_time = T, light_ON = 0)
  E_EON <- E_night_vector[length(E_night_vector)] # == E_diel_acclim(TIME = t_EON, E_mean, A_SQW, A1, phi1, A2, phi2, acclim_slope, acclim_range, use_shifted_time = T, light_ON = 0)
  
  # Minimal nighttime transpiration
  night_min <- which.min(E_night_vector) # index
  E_night_min <- E_night_vector[night_min] # mmol/m2/s
  
  # Time when minimal nighttime transpiration occurs
  t_night_min <- time_night_vector[night_min] # fraction of day so that Pho_Per/24 [e.g. 0.5] is the start of the night and 1 [i.e. Pho_Per/24+Sko_Per/24] the end of the night
  
  # Average nighttime transpiration
  E_night_mean <- mean(E_night_vector) # mmol/m2/s
  
  # Average slope of change in transpiration throughout the nighttime after 30-min darkness (sigma_night),
  # equals to the integral of the time-derivative of transpiration (so back to the E function) divided by time interval,
  # which also equals to the linear slope between the start (after the square wave) and the end of the night
  sigma_night <- (E_EON - E_SON) / (t_EON - t_SON) # mmol/m2/s/d
  # and can also be estimated by averaging the time-derivative over the time interval:
  #night_deriv <- A1*2*pi*cos((time_night_vector+phi1)*2*pi) + A2*4*pi*cos((time_night_vector+phi2)*4*pi) + acclim_slope
  #sigma_night <- mean(night_deriv)
  
  # Cumulative transpiration dynamics above the nighttime trend (Delta_night, "absolute nighttime dynamics"),
  # equals to the area above (usually, below) the trapezoid formed by E_SON and E_EON (i.e. the nocturnal trend),
  # which is usually negative (yet some parts may be positive if the curve is concave)
  primitive_t_SON <- -(A1/(2*pi))*cos((t_SON+phi1)*2*pi) - (A2/(4*pi))*cos((t_SON+phi2)*4*pi) + (E_mean-A_SQW)*t_SON + 0.5*acclim_slope*t_SON^2 - acclim_slope*(range_pulse+SQW_shift/60/24)*t_SON
  primitive_t_EON <- -(A1/(2*pi))*cos((t_EON+phi1)*2*pi) - (A2/(4*pi))*cos((t_EON+phi2)*4*pi) + (E_mean-A_SQW)*t_EON + 0.5*acclim_slope*t_EON^2 - acclim_slope*(range_pulse+SQW_shift/60/24)*t_EON
  trapezoid_area <- (E_SON+E_EON)*(t_EON-t_SON)/2
  Delta_night <- primitive_t_EON - primitive_t_SON - trapezoid_area # (mmol/m2/s) * d
  # that can also be computed by numerical integration like in 'E_val_obs()':
  #Delta_night <- trapz(time_night_vector, E_night_vector) - trapz(c(Pho_Per/24+Rapid_dur/60/24, Pho_Per/24+Sko_Per/24), c(E_SON, E_EON)) # (mmol/m2/s) * d
  Delta_night <- Delta_night*24*60*60/1000 # mol/m2
  
  # Cumulative preopening once minimal transpiration has been reached (Sigma_preop),
  # equals to the integral of transpiration between tmin and tend
  # (y-translated so that transpiration is 0 at tmin)
  primitive_t_night_min <- -(A1/(2*pi))*cos((t_night_min+phi1)*2*pi) - (A2/(4*pi))*cos((t_night_min+phi2)*4*pi) + (E_mean-A_SQW)*t_night_min + 0.5*acclim_slope*t_night_min^2 - acclim_slope*(range_pulse+SQW_shift/60/24)*t_night_min
  rectangle_area <- E_night_min*(t_EON-t_night_min)
  Sigma_preop <- primitive_t_EON - primitive_t_night_min - rectangle_area # (mmol/m2/s) * d
  # which is equal to a similar calculus as for the day:
  #primitive_t_night_min_transl <- -(A1/(2*pi))*cos((t_night_min+phi1)*2*pi) - (A2/(4*pi))*cos((t_night_min+phi2)*4*pi) + (E_mean-A_SQW-E_night_min)*t_night_min + 0.5*acclim_slope*t_night_min^2 - acclim_slope*(range_pulse+SQW_shift/60/24)*t_night_min
  #primitive_t_EON_transl <- -(A1/(2*pi))*cos((t_EON+phi1)*2*pi) - (A2/(4*pi))*cos((t_EON+phi2)*4*pi) + (E_mean-A_SQW-E_night_min)*t_EON + 0.5*acclim_slope*t_EON^2 - acclim_slope*(range_pulse+SQW_shift/60/24)*t_EON
  #Sigma_preop <- primitive_t_EON_transl - primitive_t_night_min_transl # (mmol/m2/s) * d
  # and can also be computed by numerical integration like in 'E_val_obs()':
  #Sigma_preop <- trapz(time_night_vector[time_night_vector >= t_night_min], E_night_vector[time_night_vector >= t_night_min] - E_night_min) # (mmol/m2/s) * d
  Sigma_preop <- Sigma_preop*24*60*60/1000 # mol/m2
  

  ## OUTPUT ##
  
  return (list(E_EON_mod = E_EON,
               E_day_mean_mod = E_day_mean, E_day_max_mod = E_day_max, t_day_max_mod = t_day_max,
               sigma_day_mod = sigma_day, Delta_day_mod = Delta_day, Sigma_preclo_mod = Sigma_preclo, #rel_dyn_mod = rel_dyn,
               E_night_mean_mod = E_night_mean, E_night_min_mod = E_night_min, t_night_min_mod = t_night_min,
               sigma_night_mod = sigma_night, Delta_night_mod = Delta_night, Sigma_preop_mod = Sigma_preop))
  
  }


# Function to extract some numeric values of interest
# for transpiration based on the observed, experimental data
E_val_obs <- function (xy, # the dataframe with the time (x) and transpiration (y) values, appended with information from 'get_cycles()'
                       E_EON_before) # value of transpiration at the end of the previous night
  {
  #require("pracma") # for 'trapz()'
  
  ## DAY ##
  
  xy_day <- xy[xy$DN == "DAY", ] # xy[xy$time_fraction_min <= Pho_Per*60, ]

  E_SOD <- E_EOD <- NA
  E_day_mean <- E_day_max <- t_day_max <- sigma_day <- Delta_day <- Sigma_preclo <- NA #rel_dyn <-
  if (length(which(!is.na(xy_day$y))) > 0)
    {
    # Average by 'time_fraction_min' (there could be several periods)
    day_mean <- aggregate(list(y = xy_day$y), by = list(x = xy_day$time_fraction_min), FUN = mean, na.rm = T)

    # Get the start and end values
    E_SOD <- if (Rapid_dur %in% day_mean$x) day_mean$y[day_mean$x == Rapid_dur] else NA # transpiration at the start of the day (30 min after the night-to-day transition)
    E_EOD <- if ((Pho_Per*60) %in% day_mean$x) day_mean$y[day_mean$x == Pho_Per*60] else NA # transpiration at the end of the day
    
    if (length(which(!is.na(xy_day$y))) > 10)
      {
      # Back to days
      day_mean$x <- day_mean$x/60/24
      
      # Maximal daytime transpiration
      day_max <- which.max(day_mean$y) # index
      E_day_max <- day_mean$y[day_max] # mmol/m2/s
      
      # Time when maximal daytime transpiration occurs
      t_day_max <- day_mean$x[day_max] # fraction of day so that 0 is the start of the day and Pho_Per/24 [i.e. 0.5] the end of the day
      
      # Average daytime transpiration
      E_day_mean <- mean(day_mean$y, na.rm = T) # mmol/m2/s
      
      # Average slope of change in transpiration throughout the daytime after 30-min illumination (sigma_day)
      sigma_day <- (E_EOD - E_SOD) / (Pho_Per/24 - Rapid_dur/60/24) # mmol/m2/s/d
      
      # Cumulative transpiration dynamics above the daytime trend (Delta_day, "absolute daytime dynamics"),
      # equals to the area above the trapezoid formed by E_SOD and E_EOD (i.e. the diurnal trend),
      # which is usually positive (yet some parts may be negative if the curve is convex)
      y_dyn <- day_mean$y
      x_dyn <- day_mean$x
      x_dyn <- x_dyn[!is.na(y_dyn)]
      y_dyn <- y_dyn[!is.na(y_dyn)]
      Delta_day <- trapz(x_dyn, y_dyn) - trapz(c(Rapid_dur/60/24, Pho_Per/24), c(E_SOD, E_EOD)) # (mmol/m2/s) * d
      Delta_day <- Delta_day*24*60*60/1000 # mol/m2
        
      # # (not used anymore)
      # # "Relative daytime dynamics":
      # # algebraic sum (may be positive or negative) of
      # # the (positive) area above E_SOD (i.e. the initial baseline) and
      # # the (negative) area below E_SOD
      # # Note that trapz(c(Rapid_dur/60/24, Pho_Per/24), c(E_SOD, E_SOD)) == E_SOD*(Pho_Per/24-Rapid_dur/60/24)
      # rel_dyn <- trapz(x_dyn, y_dyn) - trapz(c(Rapid_dur/60/24, Pho_Per/24), c(E_SOD, E_SOD)) # (mmol/m2/s) * d
      # rel_dyn <- rel_dyn*24*60*60/1000 # mol/m2
      
      # Cumulative preclosure once maximal transpiration has been reached (Sigma_preclo),
      y_clo <- day_mean$y[day_max:nrow(day_mean)] - day_mean$y[day_max]
      x_clo <- day_mean$x[day_max:nrow(day_mean)]
      x_clo <- x_clo[!is.na(y_clo)]
      y_clo <- y_clo[!is.na(y_clo)]
      Sigma_preclo <- -trapz(x_clo, y_clo) # (mmol/m2/s) * d
      Sigma_preclo <- Sigma_preclo*24*60*60/1000 # mol/m2
      }
    }
  
  # If extrema have been removed by 'remove_incomplete_periods()', get them back using the backup
  if (is.na(E_SOD))
    if (Rapid_dur %in% xy$time_fraction_min)
      if (length(which(!is.na(xy$y_backup[xy$time_fraction_min %in% Rapid_dur]))) > 0)
        E_SOD <- mean(xy$y_backup[xy$time_fraction_min %in% Rapid_dur], na.rm = T)
  
  if (is.na(E_EOD))
    if ((Pho_Per*60) %in% xy$time_fraction_min)
      if (length(which(!is.na(xy$y_backup[xy$time_fraction_min %in% (Pho_Per*60)]))) > 0)
        E_SOD <- mean(xy$y_backup[xy$time_fraction_min %in% (Pho_Per*60)], na.rm = T)
  
  
  ## NIGHT ##
  
  xy_night <- xy[xy$DN == "NIGHT", ] # xy[xy$time_fraction_min > Pho_Per*60, ]

  E_SON <- E_EON <- NA
  E_night_mean <- E_night_min <- t_night_min <- sigma_night <- Delta_night <- Sigma_preop <- NA
  if (length(which(!is.na(xy_night$y))) > 0)
    {
    # Average by 'time_fraction_min' (there could be several periods)
    night_mean <- aggregate(list(y = xy_night$y), by = list(x = xy_night$time_fraction_min), FUN = mean, na.rm = T)

    # Get the start and end values
    E_SON <- if ((Pho_Per*60 + Rapid_dur) %in% night_mean$x) night_mean$y[night_mean$x == Pho_Per*60 + Rapid_dur] else NA # transpiration at the start of the night (30 min after the day-to-night transition)
    E_EON <- if ((Pho_Per*60 + Sko_Per*60) %in% night_mean$x) night_mean$y[night_mean$x == Pho_Per*60 + Sko_Per*60] else NA # transpiration at the end of the night
    
    if (length(which(!is.na(xy_night$y))) > 10)
      {
      # Back to days
      night_mean$x <- night_mean$x/60/24
      
      # Minimal nighttime transpiration
      night_min <- which.min(night_mean$y) # index
      E_night_min <- night_mean$y[night_min] # mmol/m2/s
      
      # Time when minimal nighttime transpiration occurs
      t_night_min <- night_mean$x[night_min] # fraction of day so that Pho_Per/24 [e.g. 0.5] is the start of the night and 1 [i.e. Pho_Per/24+Sko_Per/24] the end of the night
      
      # Average nighttime transpiration
      E_night_mean <- mean(night_mean$y, na.rm = T) # mmol/m2/s
      
      # Average slope of change in transpiration throughout the nighttime after 30-min darkness (sigma_night)
      sigma_night <- (E_EON - E_SON) / ((Pho_Per/24+Sko_Per/24) - (Pho_Per/24+Rapid_dur/60/24)) # mmol/m2/s/d
      
      # Cumulative transpiration dynamics above the nighttime trend (Delta_night, "absolute nighttime dynamics"),
      # equals to the area above (usually, below) the trapezoid formed by E_SON and E_EON (i.e. the nocturnal trend),
      # which is usually negative (yet some parts may be positive if the curve is concave)
      y_dyn <- night_mean$y
      x_dyn <- night_mean$x
      x_dyn <- x_dyn[!is.na(y_dyn)]
      y_dyn <- y_dyn[!is.na(y_dyn)]
      Delta_night <- trapz(x_dyn, y_dyn) - trapz(c(Pho_Per/24+Rapid_dur/60/24, Pho_Per/24+Sko_Per/24), c(E_SON, E_EON)) # (mmol/m2/s) * d
      Delta_night <- Delta_night*24*60*60/1000 # mol/m2
      
      # Cumulative preopening once minimal transpiration has been reached (Sigma_preop)
      y_reop <- night_mean$y[night_min:nrow(night_mean)] - night_mean$y[night_min]
      x_reop <- night_mean$x[night_min:nrow(night_mean)]
      x_reop <- x_reop[!is.na(y_reop)]
      y_reop <- y_reop[!is.na(y_reop)]
      Sigma_preop <- trapz(x_reop, y_reop) # (mmol/m2/s) * d
      Sigma_preop <- Sigma_preop*24*60*60/1000 # mol/m2
      }
    }
  
  # If extrema have been removed by 'remove_incomplete_periods()', get them back using the backup
  if (is.na(E_SON))
    if ((Pho_Per*60 + Rapid_dur) %in% xy$time_fraction_min)
      if (length(which(!is.na(xy$y_backup[xy$time_fraction_min %in% (Pho_Per*60 + Rapid_dur)]))) > 0)
        E_SON <- mean(xy$y_backup[xy$time_fraction_min %in% (Pho_Per*60 + Rapid_dur)], na.rm = T)
  
  if (is.na(E_EON))
    if ((Pho_Per*60 + Sko_Per*60) %in% xy$time_fraction_min)
      if (length(which(!is.na(xy$y_backup[xy$time_fraction_min %in% (Pho_Per*60 + Sko_Per*60)]))) > 0)
        E_SOD <- mean(xy$y_backup[xy$time_fraction_min %in% (Pho_Per*60 + Sko_Per*60)], na.rm = T)  
  
  
  ## DIEL ##
  
  E_mean <- mean(c(E_day_mean, E_night_mean))
  A_rapid_op_transition <- if (missing(E_EON_before)) NA else (E_SOD - E_EON_before) # rapid opening at the night-to-day transition when the environmental condition was different on the previous day (WARNING, E_SOD may have been averaged over several days, depending on the duration of the current period!)
  A_rapid_op_stable <- E_SOD - E_EON # rapid opening at the night-to-day transition in "stable" conditions (considering transpiration at the start of the day is the same from one day to the other, which may be wrong if the environmental conditions have just changed)
  # For compatibility with earlier versions, also compute 'A_rapid_op',
  # which is 'A_rapid_op_stable' for the initial condition and 'A_rapid_op_transition' otherwise:
  if (missing(E_EON_before))
    {
    E_EON_before <- E_EON # when E at the end of the current night is the same as at the end of the previous night (e.g. for the control period, when environmental conditions are stable)
    }
  A_rapid_op <- E_SOD - E_EON_before
  A_rapid_clo <- E_EOD - E_SON
  

  ## OUTPUT ##
  return (list(E_diel_obs = E_mean, A_rapid_op_obs = A_rapid_op,
               A_rapid_op_transition_obs = A_rapid_op_transition, A_rapid_op_stable_obs = A_rapid_op_stable,
               A_rapid_clo_obs = A_rapid_clo, E_EON_obs = E_EON,
               E_day_mean_obs = E_day_mean, E_day_max_obs = E_day_max, t_day_max_obs = t_day_max,
               sigma_day_obs = sigma_day, Delta_day_obs = Delta_day, Sigma_preclo_obs = Sigma_preclo, #rel_dyn_obs = rel_dyn,
               E_night_mean_obs = E_night_mean, E_night_min_obs = E_night_min, t_night_min_obs = t_night_min,
               sigma_night_obs = sigma_night, Delta_night_obs = Delta_night, Sigma_preop_obs = Sigma_preop))
  }


# Identify individual diel cycles within a larger period of stable environmental conditions,
# as well as and days and nights
get_cycles <- function (dat) # the dataframe of one pot
  {
  dat$cycle <- floor(dat$Exact_Time_min/60/24) # diel cycle as an integer
  dat$time_fraction <- dat$Exact_Time_min/60/24 - floor(dat$Exact_Time_min/60/24) # fraction of day so that 0 is the start of the day and 1 the end of the night
  if (length(which(dat$time_fraction == 0)) > 0) # the points at the night-to-day transitions...
    {
    dat$time_fraction[dat$time_fraction == 0] <- 1 # ...are part of the night (EON, end of the night)
    dat$cycle[dat$time_fraction == 1] <- dat$cycle[dat$time_fraction == 1] - 1
    }
  dat$time_fraction_min <- as.integer(round(dat$time_fraction * 60 * 24)) # Same process as in 'format_time()' in 'PhenoLeaks_outliers.R': back to integers, otherwise numeric issues are observed when comparing numbers with infinite digits (using 'aggregate()' in 'E_val_obs()')
  dat$DN <- "DAY"
  dat$DN[dat$time_fraction_min > Pho_Per*60] <- "NIGHT"
  return (dat)
  }


# Remove data of incomplete periods, namely periods that have:
#   - less than 'min_obs' observations per day or night (defaults to 10 observations),
#   - more than 'max_skew' min difference between the observed and expected mean time in the day or night (defaults to 1 h)
# If 'remove_mode' is:
#   - "global" (less stringent), the number and time-distribution of observations will be evaluated by pooling days and nights of different cycles (if several available)
#   - "day_or_night" (intermediate), only the day or night period that violates one of the rules will be removed
#   - "day_and_night" (most stringent), both the day and the consecutive night of the same cycle will be removed if one of them violates one of the rules
remove_incomplete_periods <- function (xy,
                                       min_obs = 10, # discard periods with less than 10 observations per day or night
                                       max_skew = 1*60, # discard periods with more than 1-h difference between the observed and expected mean time in the day or night
                                       remove_mode = "day_and_night") # should be "global", "day_or_night" or "day_and_night"
  {
  # Number of available observations
  if (remove_mode == "global") # less stringent - does not distinguish between individual cycles (which are pooled if several are available). Equivalent to "day_or_night" if only one cycle
    {
    N_per_DN <- aggregate(list(N = xy$y), by = list(DN = xy$DN), FUN = function (x) length(na.omit(x)))
    not_enough_data <- N_per_DN$DN[N_per_DN$N < min_obs]
    xy$y[xy$DN %in% not_enough_data] <- NA
    }
  else
    {
    N_per_period_DN <- aggregate(list(N = xy$y), by = list(cycle = xy$cycle, DN = xy$DN), FUN = function (x) length(na.omit(x)))
    if (remove_mode == "day_or_night") # intermediate, look at each cycle separately and remove only the day or night period with insufficient data
      {
      not_enough_data <- N_per_period_DN$N < min_obs
      for (i in which(not_enough_data))
        {
        xy$y[xy$cycle == N_per_period_DN$cycle[i] & xy$cycle == N_per_period_DN$DN[i]] <- NA
        }
      }
    else if (remove_mode == "day_and_night") # most stringent, look at each cycle separately and remove both the day or night period if either has insufficient data
      {
      for (cyc in N_per_period_DN$cycle) # if either day or night is absent (i.e. truly absent, not NA) from the input, then it is absent in the aggregated dataframe, so need to add it in order to identify the other as an incomplete period in 'not_enough_data'
        {
        dn <- N_per_period_DN$DN[N_per_period_DN$cycle == cyc]
        if (length(dn) == 1)
          {
          N_per_period_DN <- rbind(N_per_period_DN, data.frame(cycle = cyc, DN = c("DAY", "NIGHT")[!c("DAY", "NIGHT") %in% dn], N = 0))
          }
        }
      not_enough_data <- N_per_period_DN$cycle[N_per_period_DN$N < min_obs]
      xy$y[xy$cycle %in% not_enough_data] <- NA
      }
    else
      {
      stop("The argument 'remove_mode' should be \"global\", \"day_or_night\" or \"day_and_night\"")
      }
    }
      
  # Time-distribution of the remaining observations
  if (length(which(!is.na(xy$y))) > 0)
    {
    if (remove_mode == "global")
      {
      MeanTime_per_DN <- aggregate(list(MeanTime = xy$time_fraction_min[!is.na(xy$y)]), by = list(DN = xy$DN[!is.na(xy$y)]), FUN = mean, na.rm = T)      
      MeanTime_per_DN$MeanDiff[MeanTime_per_DN$DN == "DAY"] <- abs(MeanTime_per_DN$MeanTime[MeanTime_per_DN$DN == "DAY"] - (Rapid_dur + (Pho_Per*60-Rapid_dur)/2))
      MeanTime_per_DN$MeanDiff[MeanTime_per_DN$DN == "NIGHT"] <- abs(MeanTime_per_DN$MeanTime[MeanTime_per_DN$DN == "NIGHT"] - (Rapid_dur + Pho_Per*60 + (Sko_Per*60-Rapid_dur)/2))
      biased_time <- MeanTime_per_DN$DN[MeanTime_per_DN$MeanDiff > max_skew]
      xy$y[xy$DN %in% biased_time] <- NA
      }
    else
      {
      MeanTime_per_period_DN <- aggregate(list(MeanTime = xy$time_fraction_min[!is.na(xy$y)]), by = list(cycle = xy$cycle[!is.na(xy$y)], DN = xy$DN[!is.na(xy$y)]), FUN = mean, na.rm = T)      
      MeanTime_per_period_DN$MeanDiff[MeanTime_per_period_DN$DN == "DAY"] <- abs(MeanTime_per_period_DN$MeanTime[MeanTime_per_period_DN$DN == "DAY"] - (Rapid_dur + (Pho_Per*60-Rapid_dur)/2))
      MeanTime_per_period_DN$MeanDiff[MeanTime_per_period_DN$DN == "NIGHT"] <- abs(MeanTime_per_period_DN$MeanTime[MeanTime_per_period_DN$DN == "NIGHT"] - (Rapid_dur + Pho_Per*60 + (Sko_Per*60-Rapid_dur)/2))
      if (remove_mode == "day_or_night")
        {
        biased_time <- MeanTime_per_period_DN$MeanDiff > max_skew
        for (i in which(biased_time))
          {
          xy$y[xy$cycle == MeanTime_per_period_DN$cycle[i] & xy$DN == MeanTime_per_period_DN$DN[i]] <- NA
          }
        }
      else if (remove_mode == "day_and_night")
        {
        biased_time <- MeanTime_per_period_DN$cycle[MeanTime_per_period_DN$MeanDiff > max_skew]
        xy$y[xy$cycle %in% biased_time] <- NA
        }
      }
    }
    
  return (xy)
  }

    
#------------------------------------------------------------------------------#
#                      Run the fit and extract parameters                      #
#------------------------------------------------------------------------------#


# This function runs the fit and extract parameters
run_fit <- function (all_dat, # the whole dataframe
                     use_VPD = F, # set to TRUE when transpiration is normalized by VPD
                     Time_var = "Time", # the column name of the time axis
                     E_var = if (!use_VPD) "E_corr" else "E_corr_per_kPa", # the column name of the transpiration axis
                     Trt_var, # the column name(s) of the statistical treatment(s), e.g. "idGenotype" or c("idGenotype", "idWatering")
                     use_acclim = F, # set to TRUE to integrate an acclimation period when the environmental conditions change (and get an additional fit parameter, 'acclim_slope')
                     col.per = if (!use_acclim) ColorsPeriod else ColorsPeriod_acclim, # the dataframe containing the features of the different periods
                     linear_detrend = T, # set to TRUE to remove the linear trend before the fit - requires shifted time
                     period_for_detrend = col.per$idPeriod[1], # the reference period(s) that will be used to compute the linear trend to be subtracted to all periods
                     ...) # other arguments to be passed to the 'remove_incomplete_periods()' function ('min_obs', 'max_skew', 'remove_mode')
                          # or to the 'E_val()' and 'lm_SQW_SIN2()' functions ('use_shifted_time', 'light_ON')
  {
  # Initialize the results dataframe.
  # Note that the units related to transpiration are per kPa if the transpiration input was normalized by the VPD.
  results <- data.frame(idPot = NULL, # pot id
                        idPeriod = NULL, # diel period id (as defined in 'ColorsPeriod' at Step#00)
                        diel_trend = NULL, # the mean slope of the daytime and nighttime linear trend over the selected period(s) (mmol/m2/s/d)
                        diel_trend_baseline = NULL, # the mean intercept of the daytime and nighttime linear trend over the selected period(s) (mmol/m2/s)
                        E_diel_obs = NULL, # observed average diel transpiration (E, mmol/m2/s)
                        A_rapid_op_obs = NULL, # observed amplitude of stomatal opening at the night-to-day transition (mmol/m2/s), which is 'A_rapid_op_stable' for the initial condition and 'A_rapid_op_transition' otherwise
                        A_rapid_op_transition_obs = NULL, # observed amplitude of stomatal opening at the night-to-day transition when the environmental condition was different on the previous day, i.e. E_SOD - E_EON_before (mmol/m2/s)
                        A_rapid_op_stable_obs = NULL, # observed amplitude of stomatal opening at the night-to-day transition in "stable" conditions, i.e. E_SOD - E_EON (mmol/m2/s)
                        A_rapid_clo_obs = NULL, # observed amplitude of stomatal closure at the day-to-night transition (mmol/m2/s)
                        E_EON_before_obs = NULL, # observed transpiration at the end of the previous night and at the start of the current day (mmol/m2/s)
                        E_EON_obs = NULL, # observed transpiration at the end of the current night and at the start of the next day (mmol/m2/s)
                        E_day_mean_obs = NULL, # observed average daytime transpiration (mmol/m2/s)
                        E_day_max_obs = NULL, # observed maximal daytime transpiration (mmol/m2/s)
                        t_day_max_obs = NULL, # observed time when maximal daytime transpiration occurs (fraction of day so that 0 is the start of the day and Pho_Per/24 [e.g. 0.5] the end of the day)
                        sigma_day_obs = NULL, # observed average temporal change in daytime transpiration (mmol/m2/s/d)
                        Delta_day_obs = NULL, # observed absolute daytime dynamics, i.e. the area above the diurnal trend (mol/m2)
                        Sigma_preclo_obs = NULL, # observed cumulative daytime (afternoon) preclosure (mol/m2)
                        E_night_mean_obs = NULL, # observed average nighttime transpiration (mmol/m2/s)
                        E_night_min_obs = NULL, # observed minimal nighttime transpiration (mmol/m2/s)
                        t_night_min_obs = NULL, # observed time when minimal nighttime transpiration occurs (fraction of day so that Pho_Per/24 [e.g. 0.5] is the start of the night and 1 [i.e. Pho_Per/24+Sko_Per/24] the end of the night)
                        sigma_night_obs = NULL, # observed average temporal change in nighttime transpiration (mmol/m2/s/d)
                        Delta_night_obs = NULL, # observed absolute nighttime dynamics, i.e. the area above the nocturnal trend (mol/m2)
                        Sigma_preop_obs = NULL, # observed cumulative nighttime preopening (mol/m2)
                        t1_fit = NULL, # start of the time interval used for the fit
                        t2_fit = NULL, # end of the time interval used for the fit
                        adjR2 = NULL, # adjusted R squared of the linear model
                        RMSE = NULL, # root mean squared error
                        E_mean = NULL, # fitted mean transpiration over a 24-h period (mmol/m2/s)
                        A_SQW = NULL, # fitted semi-amplitude of the square wave (mmol/m2/s)
                        A1 = NULL, # fitted semi-amplitude of the first harmonics (mmol/m2/s)
                        phi1 = NULL, # fitted phase of the first harmonics (days)
                        A2 = NULL, # fitted semi-amplitude of the second harmonics (mmol/m2/s)
                        phi2 = NULL, # fitted phase of the second harmonics (days)
                        acclim_slope = NULL, # slope of the linear pulse during acclimation (mmol/m2/s/d)
                        E_EON_before_mod = NULL, # modelled transpiration at the end of the previous night and at the start of the current day (mmol/m2/s)
                        E_EON_mod = NULL, # modelled transpiration at the end of the current night and at the start of the next day (mmol/m2/s)
                        E_day_mean_mod = NULL, # modelled average daytime transpiration (mmol/m2/s)
                        E_day_max_mod = NULL, # modelled maximal daytime transpiration (mmol/m2/s)
                        t_day_max_mod = NULL, # modelled time when maximal daytime transpiration occurs (fraction of day so that 0 is the start of the day and Pho_Per/24 [e.g. 0.5] the end of the day)
                        sigma_day_mod = NULL, # modelled average temporal change in daytime transpiration (mmol/m2/s/d)
                        Delta_day_mod = NULL, # modelled absolute daytime dynamics, i.e. the area above the diurnal trend (mol/m2 - usually positive since most day curves are concave)
                        Sigma_preclo_mod = NULL, # modelled cumulative daytime (afternoon) preclosure (mol/m2)
                        E_night_mean_mod = NULL, # modelled average nighttime transpiration (mmol/m2/s)
                        E_night_min_mod = NULL, # modelled minimal nighttime transpiration (mmol/m2/s)
                        t_night_min_mod = NULL, # modelled time when minimal nighttime transpiration occurs (fraction of day so that Pho_Per/24 [e.g. 0.5] is the start of the night and 1 [i.e. Pho_Per/24+Sko_Per/24] the end of the night)
                        sigma_night_mod = NULL, # modelled average temporal change in nighttime transpiration (mmol/m2/s/d)
                        Delta_night_mod = NULL, # modelled absolute nighttime dynamics, i.e. the area above the nocturnal trend (mol/m2 - usually negative since most night curves are convex)
                        Sigma_preop_mod = NULL) # modelled cumulative nighttime preopening (mol/m2)
  
  for (pot in unique(all_dat$idPot))
    {
    # Select the pot data
    dat <- all_dat[all_dat$idPot == pot, ]
    
    # Get period information (diel cycles, time fraction, day/night)
    dat <- get_cycles(dat)
      
    # Detrend based on the mean of the daytime and nighttime linear trend 
    if (linear_detrend)
      {
      detrend_period <- vector()
      for (period in period_for_detrend)
        {
        t1 <- col.per[col.per$idPeriod == period, paste(Time_var, 1, sep = "")]
        t2 <- col.per[col.per$idPeriod == period, paste(Time_var, 2, sep = "")]
        t1_min <- 24 * 60 * (t1 - if (is.null(list(...)$light_ON)) 0 else list(...)$light_ON)
        t2_min <- 24 * 60 * (t2 - if (is.null(list(...)$light_ON)) 0 else list(...)$light_ON)
        
        if (period == col.per$idPeriod[1])
          {
          detrend_period <- c(detrend_period, which(dat$Exact_Time_min >= t1_min & dat$Exact_Time_min <= t2_min))
          #detrend_period <- c(detrend_period, which(dat[, Time_var] >= t1 & dat[, Time_var] <= t2))
          }
        else
          {
          detrend_period <- c(detrend_period, which(dat$Exact_Time_min > t1_min & dat$Exact_Time_min <= t2_min))
          #detrend_period <- c(detrend_period, which(dat[, Time_var] > t1 & dat[, Time_var] <= t2))
          }
        }
    
      dat_for_detrend <- dat[detrend_period, ]
      day_trend <- coef(lm(get(E_var) ~ get(Time_var), data = dat_for_detrend[dat_for_detrend$DN == "DAY", ]))
      night_trend <- coef(lm(get(E_var) ~ get(Time_var), data = dat_for_detrend[dat_for_detrend$DN == "NIGHT", ]))
      diel_trend <- mean(c(day_trend[2], night_trend[2]))
      diel_trend_baseline <- mean(c(day_trend[1], night_trend[1])) # for graphical purpose only
      }
      
    else
      {
      diel_trend_baseline <- diel_trend <- 0 
      }
    
    dat$y <- dat[, E_var] - diel_trend*dat[, Time_var]

    # Initialize the parameters
    E_EON_before_mod <- E_EON_before_obs <- NA
    E_mean <- A_SQW <- A1 <- phi1 <- A2 <- phi2 <- acclim_slope <- NA 
    RMSE <- adjR2 <- modsqsin2 <- NA 
    E_EON_mod <- E_day_mean_mod <- E_night_mean_mod <- E_day_max_mod <- t_day_max_mod <- sigma_day_mod <- Delta_day_mod <- Sigma_preclo_mod <- E_night_min_mod <- t_night_min_mod <- sigma_night_mod <- Delta_night_mod <- Sigma_preop_mod <- NA
    E_diel_obs <- A_rapid_op_obs <- A_rapid_op_transition_obs <- A_rapid_op_stable_obs <- A_rapid_clo_obs <- E_EON_obs <- E_day_mean_obs <- E_day_max_obs <- t_day_max_obs <- sigma_day_obs <- Delta_day_obs <- Sigma_preclo_obs <- E_night_mean_obs <- E_night_min_obs <- t_night_min_obs <- sigma_night_obs <- Delta_night_obs <- Sigma_preop_obs <- NA
    
    for (idperiod in 1:nrow(col.per))
      {
      # Select data of the relevant period
      # Be careful not to input a point from a night-to-day condition belonging to the previous cycle if the latter was another condition,
      # e.g. if the current diel cycle is high CO2 while the previous cycle was control,
      # do not input the point at the night-to-day transition between the control night and the high CO2 day.
      t1 <- col.per[idperiod, paste(Time_var, 1, sep = "")]
      t2 <- col.per[idperiod, paste(Time_var, 2, sep = "")]
      t1_min <- 24 * 60 * (t1 - if (is.null(list(...)$light_ON)) 0 else list(...)$light_ON)
      t2_min <- 24 * 60 * (t2 - if (is.null(list(...)$light_ON)) 0 else list(...)$light_ON)
      
      if (idperiod == 1) # initial, stable conditions (so t1 can be kept)
        {
        period <- dat$Exact_Time_min >= t1_min & dat$Exact_Time_min <= t2_min
        #period <- dat[, Time_var] >= t1 & dat[, Time_var] <= t2
        }
      else # then make sure the end of the previous night (t1) does not constrain the fit
        {
        period <- dat$Exact_Time_min > t1_min & dat$Exact_Time_min <= t2_min
        #period <- dat[, Time_var] > t1 & dat[, Time_var] <= t2
        }

      # create the dataframe
      xy <- dat[period, c(Time_var, "y", "cycle", "time_fraction_min", "DN", "Exact_Time_min")]
      colnames(xy)[1] <- "x"
      
      # Remove data of incomplete periods that could generate bias
      xy$y_backup <- xy$y # backup that 'E_val_obs()' may use to recover specific values that do not require complete dynamics
      args <- list(xy = xy, min_obs = list(...)$min_obs, max_skew = list(...)$max_skew, remove_mode = list(...)$remove_mode)
      args <- Filter(Negate(is.null), args)
      xy <- do.call(remove_incomplete_periods, args)

      # Set the fit type and the acclimation range, if required
      fit.type <- "standard"
      if (use_acclim)
        {
        t1_acclim <- col.per[idperiod, paste(Time_var, "1_acclim", sep = "")]
        if (!is.na(t1_acclim)) # note that if NA is specified for the acclimation range of the current 'idperiod' (initial, steady conditions), then the standard fit is run
          {
          t1_acclim_min <- 24 * 60 * (t1_acclim - if (is.null(list(...)$light_ON)) 0 else list(...)$light_ON)
          t2_acclim <- col.per[idperiod, paste(Time_var, "2_acclim", sep = "")]
          t2_acclim_min <- 24 * 60 * (t2_acclim - if (is.null(list(...)$light_ON)) 0 else list(...)$light_ON)
          acclim_range <- c(t1_acclim, t2_acclim)
          fit.type <- if (t2_acclim == t2) "acclimating" else c("acclimating", "acclimated") #  need to distinguish between "acclimating" and "acclimated" because the extracted parameters by 'E_val()' and 'E_val_obs()' won't be the same in both periods (yet there is only one fit for both periods, so the fitted parameters are the same)
          }
        }
       
      # Run the fit if enough data
      check.fit <- F
      if (length(which(!is.na(xy$y))) > 0)
        {
        # Check if data are available in the day and in the night (useful only if 'remove_mode' was "global" or "day_or_night")
        N_per_DN <- aggregate(list(N = xy$y), by = list(DN = xy$DN), FUN = function (x) length(na.omit(x)))
        if (nrow(N_per_DN) == 2 & min(N_per_DN$N) > 0) # at least 10 observations in the day and in the night
          {
          # Potentially get additional arguments to be passed to 'E_val()' and 'lm_SQW_SIN2()' (or their '_acclim' version)
          args <- list(use_shifted_time = list(...)$use_shifted_time, light_ON = list(...)$light_ON)
          args <- Filter(Negate(is.null), args)
          
          # Square + sine (first and second harmonics) wave 
          if (fit.type[1] == "standard")
            {
            modsqsin2 <- do.call(lm_SQW_SIN2, append(list(xy = xy), args))
            #modsqsin2 <- lm_SQW_SIN2(xy, ...)
            }
          else # add an 'acclimation pulse' over the specified range
            {
            modsqsin2 <- do.call(lm_SQW_SIN2_acclim, append(list(xy = xy, acclim_range = acclim_range), args))
            }
          adjR2 <- summary(modsqsin2)$adj.r.squared
          RMSE <- sqrt(mean((xy$y - predict(modsqsin2))^2, na.rm = T))
          #c(E_mean, A_SQW, A1, phi1, A2, phi2) := retrieve_pars(modsqsin2)
          c(E_mean, A_SQW, A1, phi1, A2, phi2, acclim_slope) := retrieve_pars(modsqsin2)
          check.fit <- T
          }
        }
      
      for (f in fit.type)
        {
        # Extract parameters - modelled
        if (check.fit)
          {
          if (f != "acclimating")
            {
            c(E_EON_mod, E_day_mean_mod, E_day_max_mod, t_day_max_mod, sigma_day_mod, Delta_day_mod, Sigma_preclo_mod, E_night_mean_mod, E_night_min_mod, t_night_min_mod, sigma_night_mod, Delta_night_mod, Sigma_preop_mod) := do.call(E_val, append(list(E_mean, A_SQW, A1, phi1, A2, phi2), args))
            #c(E_EON_mod, E_day_mean_mod, E_day_max_mod, t_day_max_mod, sigma_day_mod, Delta_day_mod, Sigma_preclo_mod, E_night_mean_mod, E_night_min_mod, t_night_min_mod, sigma_night_mod, Delta_night_mod, Sigma_preop_mod) := E_val(E_mean, A_SQW, A1, phi1, A2, phi2, ...)
            }
          else
            {
            c(E_EON_mod, E_day_mean_mod, E_day_max_mod, t_day_max_mod, sigma_day_mod, Delta_day_mod, Sigma_preclo_mod, E_night_mean_mod, E_night_min_mod, t_night_min_mod, sigma_night_mod, Delta_night_mod, Sigma_preop_mod) := do.call(E_val, append(list(E_mean, A_SQW, A1, phi1, A2, phi2, acclim_slope, acclim_range), args))  
            }
          }
      
        # Extract parameters - observed
        if (f == "standard")
          {
          c(E_diel_obs, A_rapid_op_obs, A_rapid_op_transition_obs, A_rapid_op_stable_obs, A_rapid_clo_obs, E_EON_obs, E_day_mean_obs, E_day_max_obs, t_day_max_obs, sigma_day_obs, Delta_day_obs, Sigma_preclo_obs, E_night_mean_obs, E_night_min_obs, t_night_min_obs, sigma_night_obs, Delta_night_obs, Sigma_preop_obs) := if (idperiod == 1) E_val_obs(xy) else E_val_obs(xy, E_EON_before_obs)
          }
        else
          {
          if (f == "acclimating")
            {
            subperiod <- xy$Exact_Time_min > t1_acclim_min & xy$Exact_Time_min <= t2_acclim_min
            }
          else # f == "acclimated"
            {
            subperiod <- xy$Exact_Time_min > t2_acclim_min
            }
          c(E_diel_obs, A_rapid_op_obs, A_rapid_op_transition_obs, A_rapid_op_stable_obs, A_rapid_clo_obs, E_EON_obs, E_day_mean_obs, E_day_max_obs, t_day_max_obs, sigma_day_obs, Delta_day_obs, Sigma_preclo_obs, E_night_mean_obs, E_night_min_obs, t_night_min_obs, sigma_night_obs, Delta_night_obs, Sigma_preop_obs) := if (idperiod == 1) E_val_obs(xy[subperiod, ]) else E_val_obs(xy[subperiod, ], E_EON_before_obs)
          }

        # Save the results
        results <- rbind(results, data.frame(idPot = pot,
                                             idPeriod = paste(col.per$idPeriod[idperiod], if (f == "standard") "" else paste(" ", f, sep = ""), sep = ""),
                                             diel_trend = diel_trend, diel_trend_baseline = diel_trend_baseline,
                                             E_diel_obs = E_diel_obs,
                                             A_rapid_op_obs = A_rapid_op_obs, A_rapid_op_transition_obs = A_rapid_op_transition_obs, A_rapid_op_stable_obs = A_rapid_op_stable_obs,
                                             A_rapid_clo_obs = A_rapid_clo_obs,
                                             E_EON_before_obs = E_EON_before_obs, E_EON_obs = E_EON_obs,
                                             E_day_mean_obs = E_day_mean_obs, E_day_max_obs = E_day_max_obs, t_day_max_obs = t_day_max_obs,
                                             sigma_day_obs = sigma_day_obs, Delta_day_obs = Delta_day_obs, Sigma_preclo_obs = Sigma_preclo_obs,
                                             E_night_mean_obs = E_night_mean_obs, E_night_min_obs = E_night_min_obs, t_night_min_obs = t_night_min_obs,
                                             sigma_night_obs = sigma_night_obs, Delta_night_obs = Delta_night_obs, Sigma_preop_obs = Sigma_preop_obs,
                                             t1_fit = t1, t2_fit = t2, adjR2 = adjR2, RMSE = RMSE,
                                             E_EON_before_mod = E_EON_before_mod, E_EON_mod = E_EON_mod,
                                             E_mean = E_mean, A_SQW = A_SQW,
                                             A1 = A1, phi1 = phi1, A2 = A2, phi2 = phi2,
                                             acclim_slope = if (f != "acclimated") acclim_slope else 0, 
                                             E_day_mean_mod = E_day_mean_mod, E_day_max_mod = E_day_max_mod, t_day_max_mod = t_day_max_mod,
                                             sigma_day_mod = sigma_day_mod, Delta_day_mod = Delta_day_mod, Sigma_preclo_mod = Sigma_preclo_mod,
                                             E_night_mean_mod = E_night_mean_mod, E_night_min_mod = E_night_min_mod, t_night_min_mod = t_night_min_mod,
                                             sigma_night_mod = sigma_night_mod, Delta_night_mod = Delta_night_mod, Sigma_preop_mod = Sigma_preop_mod)) 
        
        # Re-initialize the parameters calculated at each subperiod
        E_EON_before_mod <- E_EON_mod
        E_EON_before_obs <- E_EON_obs
        E_EON_mod <- E_day_mean_mod <- E_night_mean_mod <- E_day_max_mod <- t_day_max_mod <- sigma_day_mod <- Delta_day_mod <- Sigma_preclo_mod <- E_night_min_mod <- t_night_min_mod <- sigma_night_mod <- Delta_night_mod <- Sigma_preop_mod <- NA
        E_diel_obs <- A_rapid_op_obs <- A_rapid_op_transition_obs <- A_rapid_op_stable_obs <- A_rapid_clo_obs <- E_EON_obs <- E_day_mean_obs <- E_day_max_obs <- t_day_max_obs <- sigma_day_obs <- Delta_day_obs <- Sigma_preclo_obs <- E_night_mean_obs <- E_night_min_obs <- t_night_min_obs <- sigma_night_obs <- Delta_night_obs <- Sigma_preop_obs <- NA
        }
      
      # Re-initialize the fit parameters
      E_mean <- A_SQW <- A1 <- phi1 <- A2 <- phi2 <- acclim_slope <- NA
      RMSE <- adjR2 <- modsqsin2 <- NA
      }
    }
  
  # Add treatment information, move to the first column(s) and return
  results <- merge(results, all_dat[!duplicated(all_dat$idPot), c("idPot", Trt_var)], by = "idPot")
  results <- results[c((ncol(results)-length(Trt_var)+1):ncol(results), 1:(ncol(results)-length(Trt_var)))]
  return (results)
  
  }


# This function adds the fit over an existing plot for one given pot
curve_fit <- function (dat, # the pot dataframe
                       use_VPD = F, # set to TRUE when transpiration is normalized by VPD
                       Time_var = "Time", # the column name of the time axis
                       E_var = if (!use_VPD) "E_corr" else "E_corr_per_kPa", # the column name of the transpiration axis
                       results, # the dataframe containing the results of the fit
                       use_acclim = F, # set to TRUE to integrate an acclimation period when the environmental conditions change (and get an additional fit parameter, 'acclim_slope')
                       col.per = if (!use_acclim) ColorsPeriod else ColorsPeriod_acclim, # the dataframe containing the features of the different periods
                       export_PPTX = F,
                       add_points = T,
                       acclim_only = F,
                       ...) # other arguments to be passed to the 'remove_incomplete_periods()' function ('min_obs', 'max_skew', 'remove_mode')
                            # or to the 'E_diel()' function ('use_shifted_time', 'light_ON')
  {
  # Retrieve the pot id
  pot <- unique(dat$idPot)
  
  # Get period information (diel cycles, time fraction, day/night)
  dat <- get_cycles(dat)
  
  for (idperiod in 1:nrow(col.per))
    {
    # Select data of the relevant period
    t1 <- col.per[idperiod, paste(Time_var, 1, sep = "")]
    t2 <- col.per[idperiod, paste(Time_var, 2, sep = "")]
    t1_min <- 24 * 60 * (t1 - if (is.null(list(...)$light_ON)) 0 else list(...)$light_ON)
    t2_min <- 24 * 60 * (t2 - if (is.null(list(...)$light_ON)) 0 else list(...)$light_ON)

    if (idperiod == 1)
      {
      period <- dat$Exact_Time_min >= t1_min & dat$Exact_Time_min <= t2_min
      #period <- dat[, Time_var] >= t1 & dat[, Time_var] <= t2
      }
    else
      {
      period <- dat$Exact_Time_min > t1_min & dat$Exact_Time_min <= t2_min
      #period <- dat[, Time_var] > t1 & dat[, Time_var] <= t2
      }
    
    # create the dataframe
    xy <- dat[period, c(Time_var, E_var, "cycle", "time_fraction_min", "DN")]
    colnames(xy)[1:2] <- c("x", "y")
    
    # Remove data of incomplete periods that could generate bias
    args <- list(xy = xy, min_obs = list(...)$min_obs, max_skew = list(...)$max_skew, remove_mode = list(...)$remove_mode)
    args <- Filter(Negate(is.null), args)
    xy <- do.call(remove_incomplete_periods, args)

    # Plot the data used for the fit
    if (length(which(!is.na(xy$y))) > 0 & add_points)
      {
      N_per_DN <- aggregate(list(N = xy$y), by = list(DN = xy$DN), FUN = function (x) length(na.omit(x)))
      if (nrow(N_per_DN) == 2 & min(N_per_DN$N) > 0) # at least 10 observations in the day and in the night
        {
        points(xy$x, xy$y, col = col.per$col[idperiod],
               pch = if (!export_PPTX) 1 else 16,
               cex = if (!export_PPTX) 1 else 0.5)
        }
      }

    
    # Set the fit type and the acclimation range, if required
    fit.type <- "standard"
    if (use_acclim)
      {
      t1_acclim <- col.per[idperiod, paste(Time_var, "1_acclim", sep = "")]
      if (!is.na(t1_acclim)) # note that if NA is specified for the acclimation range of the current 'idperiod' (initial, steady conditions), then the standard fit has been run
        {
        t2_acclim <- col.per[idperiod, paste(Time_var, "2_acclim", sep = "")]
        acclim_range <- c(t1_acclim, t2_acclim)
        fit.type <- "acclim" # no need to specify c("acclimating", "acclimated") here because there is one fit for both periods (so the fitted parameters are the same)
        }
      }
       
    # Retrieve fitted parameters
    c(diel_trend, E_mean, A_SQW, A1, phi1, A2, phi2, acclim_slope, E_EON_before_mod) := results[results$idPot == pot & results$idPeriod == paste(col.per$idPeriod[idperiod], if (fit.type == "standard") "" else " acclimating", sep = ""),
                                                                                                c("diel_trend", "E_mean", "A_SQW", "A1", "phi1", "A2", "phi2", "acclim_slope", "E_EON_before_mod")]
    
    # Potentially get additional arguments to be passed to 'E_diel()'
    args <- list(E_mean, A_SQW, A1, phi1, A2, phi2, use_shifted_time = list(...)$use_shifted_time, light_ON = list(...)$light_ON)
    args <- Filter(Negate(is.null), args)

    # Add the fitted line
    if (fit.type == "standard")
      {
      if (!acclim_only)
        {
        ALPHA <- 1 #if (!export_PPTX) 1 else 0.5
        LWD <- 1 #if (!export_PPTX) 1 else 3
        if (idperiod == 1)
          {
          curve(do.call(E_diel, append(list(TIME = x), args)) + diel_trend*x, t1, t2, add = T, #col = col.per$col[idperiod])
                col = rgb(t(col2rgb(col.per$col[idperiod]))/255, alpha = ALPHA), lwd = LWD)
          #curve(E_diel(TIME = x, E_mean, A_SQW, A1, phi1, A2, phi2, ...) + diel_trend*x, t1, t2, add = T, col = col.per$col[idperiod])
          }
        else # the conditions have changed, and t1 has not constrained the fit (see the 'run_fit()' function)...
          {
          # ... so we plot the new curve just after 'SQW_shift'...
          curve(do.call(E_diel, append(list(TIME = x), args)) + diel_trend*x, t1 + SQW_shift/60/24 + 1e-10, t2, add = T, #col = col.per$col[idperiod])
                col = rgb(t(col2rgb(col.per$col[idperiod]))/255, alpha = ALPHA), lwd = LWD)
          #curve(E_diel(TIME = x, E_mean, A_SQW, A1, phi1, A2, phi2, ...) + diel_trend*x, t1 + SQW_shift/60/24 + 1e-10, t2, add = T, col = col.per$col[idperiod])
          
          # ... and for visual display we linearly interpolate between the end of the previous night and the new day conditions
          segments(t1, E_EON_before_mod + diel_trend*t1, t1 + SQW_shift/60/24 + 1e-10, do.call(E_diel, append(list(TIME = t1 + SQW_shift/60/24 + 1e-10), args)) + diel_trend*(t1 + SQW_shift/60/24 + 1e-10), lty = 2, #col = col.per$col[idperiod])
                   col = rgb(t(col2rgb(col.per$col[idperiod]))/255, alpha = ALPHA), lwd = LWD)
          #segments(t1, E_EON_before_mod + diel_trend*t1, t1 + SQW_shift/60/24 + 1e-10, E_diel(t1 + SQW_shift/60/24 + 1e-10, E_mean, A_SQW, A1, phi1, A2, phi2, ...) + diel_trend*(t1 + SQW_shift/60/24 + 1e-10), col = col.per$col[idperiod], lty = 2)
          }
        }
      }
    
    else
      {
      ALPHA <- if (!export_PPTX) 1 else 0.5
      LWD <- if (!export_PPTX) 1 else 4
      if (is.na(acclim_slope) ) # may occur if 'lm_SQW_SIN2_acclim()' has been fitted using only "acclimated" data but no "acclimating" data within 'acclim_range'
        {
        acclim_slope <- 0
        t1 <- t2_acclim
        E_EON_before_mod <- NA
        ALPHA <- 1
        LWD <- 1
        }
      args <- append(args, list(acclim_slope = acclim_slope, acclim_range = acclim_range))
      if (idperiod == 1)
        {
        curve(do.call(E_diel_acclim, append(list(TIME = x), args)) + diel_trend*x, t1, t2, add = T, #col = col.per$col[idperiod])
              col = rgb(t(col2rgb(col.per$col[idperiod]))/255, alpha = ALPHA), lwd = LWD)
        }
      else # the conditions have changed, and t1 has not constrained the fit (see the 'run_fit()' function)...
        {
        # ... so we plot the new curve just after 'SQW_shift'...
        curve(do.call(E_diel_acclim, append(list(TIME = x), args)) + diel_trend*x, t1 + SQW_shift/60/24 + 1e-10, t2, add = T, #col = col.per$col[idperiod])
              col = rgb(t(col2rgb(col.per$col[idperiod]))/255, alpha = ALPHA), lwd = LWD)
          
        # ... and for visual display we linearly interpolate between the end of the previous night and the new day conditions
        LWD <- if (LWD != 1) LWD/2 else LWD
        segments(t1, E_EON_before_mod + diel_trend*t1, t1 + SQW_shift/60/24 + 1e-10, do.call(E_diel_acclim, append(list(TIME = t1 + SQW_shift/60/24 + 1e-10), args)) + diel_trend*(t1 + SQW_shift/60/24 + 1e-10), lty = 2, #col = col.per$col[idperiod])
                 col = rgb(t(col2rgb(col.per$col[idperiod]))/255, alpha = ALPHA), lwd = LWD)
        }
      }

    }
  }

