################################################################################
#                                                                              #
#                             PhenoLeaks - OUTLIERS                            #
#                                                                              #
#                Core functions to automatically detect outliers               #
#                      in the processed transpiration data                     #
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
for (pkg in c("dplyr", "broom")) #"rstatix"
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
#    Outliers based on the variability of each treatment at each time point    #
#------------------------------------------------------------------------------#

# This function detects outliers based on the transpiration difference between adjacent time points,
# which are obtained using different lags (here 30 min, 1 h, 1.5 h and 2h by default).
# Each difference is compared to that of other plants of the same treatment (i.e. genotype, watering regime)
# for the same time point.
#
# Here is the rationale.
# Let's consider a point at time t. For each lag, it will be compared twice (with point t-lag and point t+lag).
# For instance, let's consider the simple case where the point t is an isolated outlier above the baseline.
# Then a strong transpiration difference will be detected for both comparisons, positive for t vs. t-lag and negative for t+lag vs. t.
# So t will be penalized twice, while t-lag and t+lag will be penalized only once.
# This procedure is repeated across different lags because outliers may not be isolated,
# and because locally missing data may prevent us to compute the transpiration difference for some points.
# Then the sum of the "penalties" gives a likelihood that the point is an outlier.

# The implementation is as follows. For each specified lag:
#   - the transpiration difference are computed using the current lag
#   - a linear model per treatment is fitted on the transpiration differences
#     against time and plant as explanatory (qualitative) variables
#   - for each point, the standardized residual (positive or negative) is subtracted to (or summed with)
#     that of the other point that have been used for computing the difference
#     (because each point is used twice, except at the extrema).
# Then, the absolute values of the standardized residual differences are summed for the different lags.
# The sum of the "penalties" is reported in the column 'sum_out', which is appended to the original dataframe.

detect_outliers <- function (all_dat, # the whole dataframe
                             E_var = "E_clean", # the column name of the transpiration axis
                             Trt_var, # the column name(s) of the statistical treatment(s), e.g. "idGenotype" or c("idGenotype", "idWatering")
                             lags = 1:4) # the different lags (vector of positive integers) to be applied for computing the difference between two points
  {
  n <- nrow(all_dat)
  new_pot <- which(diff(all_dat$Exact_Time_min) < 0) + 1
  all_dat$TimeFactor <- as.factor(all_dat$Exact_Time_min)
  
  for (lag in lags)
    {
    # Compute the difference between two points that are separated by (default):
    #   - 0 point (lag = 1, consecutive, here 30 min)
    #   - 1 point (lag = 2, here 1 h)
    #   - 2 points (lag = 3, here 1.5 h)
    #   - 3 points (lag = 4, here 2 h)
    all_dat$diff <- c(rep(NA, lag), diff(all_dat[, E_var], lag = lag))
    ini_new_pot <- vector()
    for (l in 0:(lag-1)) { ini_new_pot <- c(ini_new_pot, new_pot + l) }
    all_dat$diff[ini_new_pot] <- NA

    # Linear models for each lag (takes a few seconds per lag)
    #require(dplyr)
    #require(broom)
    model_lag <- all_dat %>% 
      group_by_at(Trt_var) %>%
      do(augment(lm(diff ~ TimeFactor + idPot, data = ., na.action = "na.exclude"), data = .))
    model_lag <- model_lag[order(model_lag$idObs), ]

    # Get the difference between the standardized residual of the current point and that of the other point
    # (except for NA values including the end of each pot, where we assign 0 - meaning a point preceding or following a missing value
    # has less chance to be detected as an outlier)
    model_lag$.std.resid[is.na(model_lag$.std.resid)] <- 0
    # The line below works well for small islands of "stable" outliers, but not for groups of outliers that steadily increase or decrease and tend to cancel each other
    #all_dat[ncol(all_dat)+1] <- c(-diff(model_lag$.std.resid, lag = lag), model_lag$.std.resid[(n-(lag-1)):n])
    # Instead, get the difference or the sum depending on the sign so that penalties always accumulate:
    plus_or_minus <- sign(model_lag$.std.resid[1:(n-lag)]) / sign(model_lag$.std.resid[(1+lag):n])
    plus_or_minus[plus_or_minus %in% c(-Inf, Inf, NaN)] <- 0
    cumul_penalties <- model_lag$.std.resid[1:(n-lag)] + plus_or_minus * model_lag$.std.resid[(1+lag):n]
    all_dat[ncol(all_dat)+1] <- c(cumul_penalties, model_lag$.std.resid[(n-(lag-1)):n])
    }

  # Sum the penalties
  sum_out <- rep(0, n)
  for (l in 1:length(lags)) { sum_out <- sum_out + abs(all_dat[, ncol(all_dat) - length(lags) + l]) }
  all_dat$sum_out <- sum_out
  
  # Clean and return the augmented 'all_dat'
  all_dat <- all_dat[, -c(which(colnames(all_dat) %in% c("TimeFactor", "diff")), ncol(all_dat) - 1:length(lags))]
  return (all_dat)
  }



#------------------------------------------------------------------------------#
#      Identification of the points at the day/night transitions and after     #
#------------------------------------------------------------------------------#

identify_transitions <- function (all_dat) # the whole dataframe
  {
  all_dat$dark2light <- round(all_dat$Exact_Time_min/60/24) == all_dat$Exact_Time_min/60/24
  all_dat$light2dark <- round(all_dat$Exact_Time_min/60/24 - Sko_Per/24) == all_dat$Exact_Time_min/60/24 - Sko_Per/24
  all_dat$transitions <- all_dat$dark2light | all_dat$light2dark
  
  all_dat$after_dark2light <- round((all_dat$Exact_Time_min - Rapid_dur)/60/24) == (all_dat$Exact_Time_min - Rapid_dur)/60/24
  all_dat$after_light2dark <- round((all_dat$Exact_Time_min - Rapid_dur)/60/24 - Sko_Per/24) == (all_dat$Exact_Time_min - Rapid_dur)/60/24 - Sko_Per/24
  all_dat$after_transitions <- all_dat$after_dark2light | all_dat$after_light2dark

  return(all_dat)
  }



#------------------------------------------------------------------------------#
#      Re-estimation of transpiration values at the day/night transitions      #
#------------------------------------------------------------------------------#

# This function re-estimates transpiration values at the transitions (start of the day or night, end of the day or night)
# in the whole dataset. It extrapolates (backwards or forwards) a linear model of transpiration against time fitted over
# the user-defined interval 'interval_est'. The data are taken:
#   - after the start of the day period to estimate transpiration 30 min after the start of the day (E_SOD)
#   - after the start of the night period to estimate transpiration 30 min after the start of the night (E_SON)
#   - before the end of the day period to estimate transpiration at the end of the day (E_EOD), i.e. the day-to-night transition
#   - before the end of the night period to estimate transpiration at the end of the night (E_EON), i.e. the night-to-day transition

estim_transitions <- function (all_dat, # the whole dataframe
                               Time_var = "Time", # the column name of the time axis
                               E_var = "E_clean", # the column name of the transpiration axis
                               time_start = min(all_dat[, Time_var], na.rm = T), # time at which outlier detection starts
                               time_end = max(all_dat[, Time_var], na.rm = T), # time at which outlier detection ends
                               interval_est = 0.15, # (d)
                               light_ON = if (Time_var == "decimalDay") Time_ON0 - floor(Time_ON0) else 0, # time of the diel cycle when the light switches on
                               allow_extrapolation = F, # set to TRUE to allow extrapolation, i.e. estimate transitions beyond the period covered by the measurements (NOT TESTED)
                               Trt_var) # only if allow_extrapolation = TRUE - the column name(s) of the statistical treatment(s), e.g. "idGenotype" or c("idGenotype", "idWatering")
  {
  # Prepare the dataframe
  transition_estimates <- data.frame(idPot = NULL, Time_var = NULL, E_estimate = NULL)
  
  # Set the time interval so that it starts and ends when the light switches on
  t1 <- if (is.integer(time_start - light_ON)) time_start else floor(time_start) + light_ON
  t2 <- if (is.integer(time_end - light_ON)) time_end else ceiling(time_end) + light_ON - (Pho_Per+Sko_Per)/24

  # Iterate for each pot and each period
  for (pot in unique(all_dat$idPot))
    {
    dat <- all_dat[all_dat$idPot == pot, ]
    for (period in t1:t2)
      {
      estimate_at <- c(period + Rapid_dur/60/24,
                       period + Pho_Per/24,
                       period + Pho_Per/24 + Rapid_dur/60/24,
                       period + (Pho_Per + Sko_Per)/24)
      
      # Backward estimate of transpiration at the start of the day (30 min after the night-to-day transition)
      E_SOD_est <- NA
      dat_SOD <- dat[dat[, Time_var] > estimate_at[1] & dat[, Time_var] < (estimate_at[1] + interval_est), ]
      if (length(which(!is.na(dat_SOD[, E_var]))) >= 3)
        {
        E_SOD_est <- as.numeric(predict(lm(get(E_var) ~ get(Time_var), dat_SOD), setNames(data.frame(estimate_at[1]), Time_var)))
        }
      
      # Forward estimate of transpiration at the end of the day
      E_EOD_est <- NA
      dat_EOD <- dat[dat[, Time_var] > (estimate_at[2] - interval_est) & dat[, Time_var] < estimate_at[2], ]
      if (length(which(!is.na(dat_EOD[, E_var]))) >= 3)
        {
        E_EOD_est <- as.numeric(predict(lm(get(E_var) ~ get(Time_var), dat_EOD), setNames(data.frame(estimate_at[2]), Time_var)))
        }
      
      # Backward estimate of transpiration at the start of the night (30 min after the day-to-night transition)
      E_SON_est <- NA
      dat_SON <- dat[dat[, Time_var] > estimate_at[3] & dat[, Time_var] < (estimate_at[3] + interval_est), ]
      if (length(which(!is.na(dat_SON[, E_var]))) >= 3)
        {
        E_SON_est <- as.numeric(predict(lm(get(E_var) ~ get(Time_var), dat_SON), setNames(data.frame(estimate_at[3]), Time_var)))
        }
      
      # Forward estimate of transpiration at the end of the night
      E_EON_est <- NA
      dat_EON <- dat[dat[, Time_var] > (estimate_at[4] - interval_est) & dat[, Time_var] < estimate_at[4], ]
      if (length(which(!is.na(dat_EON[, E_var]))) >= 3)
        {
        E_EON_est <- as.numeric(predict(lm(get(E_var) ~ get(Time_var), dat_EON), setNames(data.frame(estimate_at[4]), Time_var)))
        }
      
      transition_estimates <- rbind(transition_estimates,
                                    data.frame(idPot = rep(pot, 4),
                                               Time_var = estimate_at,
                                               E_estimate = c(E_SOD_est, E_EOD_est, E_SON_est, E_EON_est)))
      }
    
    rm(dat)
    }
  
  # Rename time
  names(transition_estimates)[2] <- Time_var
  
  # Keep only the selected period
  transition_estimates <- transition_estimates[transition_estimates[, Time_var] >= time_start & transition_estimates[, Time_var] <= time_end, ]
  
  # Merge
  transition_estimates$Exact_Time_min <- as.integer(round((transition_estimates[, Time_var] - light_ON) * 60 * 24))
  transition_estimates <- transition_estimates[-2]
  all_dat <- identify_transitions(all_dat)
  all_dat <- merge(all_dat, transition_estimates, all = T)
  
  # Check that no new lines have been created.
  # This would occur if some time points have been deleted for some pots instead of having set the transpiration value to NA.
  # Step#01 normally avoids this problem by generating regular time points for all pots over the same time interval.
  # Early call of 'format_time()' at Step#02 also normally corrects this problem.
  # But this would also occur if time_start < Time_start_exp or time_end > Time_end_exp (i.e. extrapolation).
  
    ## If extrapolation is deemed acceptable, then the dataframe should be rearranged at user's own risks (NOT TESTED):
    if (length(which(is.na(all_dat$idObs))) > 0 & allow_extrapolation)
      {
      all_dat <- all_dat[order(all_dat$idPot, all_dat$Exact_Time_min), ]
      for (i in which(is.na(all_dat$idObs)))
        {
        # Note that 'meanVPD' is not calculated for the extrapolated times,
        # so the function will need to be run on both "E_clean" and "E_clean_per_kPa", and the results merged, if both values are required
        pot <- all_dat$idPot[i]
        for (trt in Trt_var)
          {
          id_trt <- na.omit(unique(all_dat[all_dat$idPot == pot, trt]))[1]
          all_dat[i, trt] <- id_trt
          }
        all_dat[i, Time_var] <- all_dat$Exact_Time_min[i] / 60 / 24 + light_ON
        all_dat$outlier[i] <- F
        }
      row.names(all_dat) <- 1:nrow(all_dat)
      all_dat$idObs <- 1:nrow(all_dat)
      all_dat <- identify_transitions(all_dat)
      }
    
    ## Check that there is no line with empty 'idobs'
    if (length(which(is.na(all_dat$idObs))) > 0) stop ("Inconsistent number of transition estimates ---> debug")
    if (length(which(diff(all_dat$idObs) != 1)) > 0) stop ("Inconsistent number of transition estimates ---> debug")
  
    ## Check that the number of estimates is consistent
    n_estimates <- length(which(!is.na(transition_estimates$E_estimate)))
    
    if (n_estimates != length(which(!is.na(all_dat$E_estimate)))) stop ("Inconsistent number of transition estimates ---> debug")
    
    if (n_estimates != length(which(!is.na(all_dat$E_estimate[all_dat$after_dark2light]))) +
        length(which(!is.na(all_dat$E_estimate[all_dat$light2dark]))) +
        length(which(!is.na(all_dat$E_estimate[all_dat$after_light2dark]))) +
        length(which(!is.na(all_dat$E_estimate[all_dat$dark2light])))) stop ("Inconsistent number of transition estimates ---> debug")
  
  # Return
  return (all_dat)
  }




#~~~~~~~~~~~~~~~~~~~~~~        OLD IMPLEMENTATION        ~~~~~~~~~~~~~~~~~~~~~~#

#------------------------------------------------------------------------------#
#    Outliers based on the variability of each treatment at each time point    #
#------------------------------------------------------------------------------#
#
# # This function detects outliers based on the transpiration difference between adjacent time points,
# # which are obtained using different lags (here 30 min, 1 h and 1.5 h).
# # Each difference is compared to that of other plants of the same treatment (i.e. genotype, watering regime)
# # for the same time point.
# # Here is the rationale.
# # Let's consider a point at time t. For each lag, it will be compared twice (with point t-lag and point t+lag).
# # For instance, let's take the minimal lag (two consecutive points). If point t is an isolated outlier,
# # then a strong transpiration difference will be detected for both comparisons, so it will take 2 penalties.
# # By contrast, the previous and next points will each receive only 1 penalty.
# # This procedure is repeated across different lags because outliers may not be isolated,
# # and because locally missing data may prevent us to compute the transpiration difference for some points.
# # Then the sum of the penalties gives a likelihood that the point is an outlier.
# #
# # The implementation is as follows. First, a linear model per treatment is fitted on the transpiration differences,
# # obtained using the different lags, against time and plant as explanatory (qualitative) variables.
# #
# # Then one of two methods may be used, for each lag and for each time point within each treatment:
# #   - method "std.resid" compare the standardized residuals of the model to a fixed threshold ('thr').
# #     It is possible (and advised) to also provide a higher threshold ('thr_after_transitions')
# #     for the points just after the day/night transitions that are intrinsically more noisy.
# #   - method "IQR" looks at the standardized residuals that are outside the 1.5 interquartile range
# #     (i.e. outliers in classical R boxplots)
# #
# # If a specific time point is outside the relevant limits for a given lag, it takes 1 penalty, otherwise it takes 0.
# # It takes 1 penalty more if the lag point is also outside the limits.
# # Note that intermediate penalties are also given when the difference cannot be computed
# # due to missing value before/after (see the function for details).
# # The sum of the penalties is reported in the column 'sum_out', which is added to the original dataframe.
# # Details of the individual components for the different lags are also added.
# detect_outliers <- function (all_dat, # the whole dataframe
#                              treatment, # the column name(s) of the statistical treatment(s), e.g. "idGenotype" or c("idGenotype", "idWatering")
#                              Time_var = if ("decimalDay" %in% names(all_dat)) "decimalDay" else "Time", # the column name of the time axis
#                              E_var = "E_clean", # the column name of the transpiration axis
#                              method, # either "std.resid" or "IQR"
#                              thr = 1, # threshold above which the standardized residual is suspicious (method "std.resid")
#                              thr_after_transitions = thr + 2) # same parameter for data just after the day/night transitions (method "std.resid")
#   {
#   # Compute the difference between two points that are separated by:
#   #   - 0 point (lag = 1, consecutive, here 30 min)
#   #   - 1 point (lag = 2, here 1 h)
#   #   - 2 points (lag = 3, here 1.5 h)
#   new_pot <- which(diff(all_dat[, Time_var]) < 0) + 1
#   all_dat$diff1 <- c(NA, diff(all_dat[, E_var], lag = 1))
#   all_dat$diff1[new_pot] <- NA
#   all_dat$diff2 <- c(NA, NA, diff(all_dat[, E_var], lag = 2))
#   all_dat$diff2[c(new_pot, new_pot+1)] <- NA
#   all_dat$diff3 <- c(NA, NA, NA, diff(all_dat[, E_var], lag = 3))
#   all_dat$diff3[c(new_pot, new_pot+1, new_pot+2)] <- NA
# 
#   # Linear models for each lag (takes a few seconds)
#   #require(dplyr)
#   #require(broom)
#   all_dat$TimeFactor <- as.factor(all_dat[, Time_var])
#   model_lag1 <- all_dat %>% 
#     group_by_at(treatment) %>%
#     do(augment(lm(diff1 ~ TimeFactor + idPot, data = ., na.action = "na.exclude"), data = .))
#   model_lag1 <- model_lag1[order(model_lag1$idObs), ]
#   model_lag2 <- all_dat %>% 
#     group_by_at(treatment) %>%
#     do(augment(lm(diff2 ~ TimeFactor + idPot, data = ., na.action = "na.exclude"), data = .))
#   model_lag2 <- model_lag2[order(model_lag2$idObs), ]
#   model_lag3 <- all_dat %>% 
#     group_by_at(treatment) %>%
#     do(augment(lm(diff3 ~ TimeFactor + idPot, data = ., na.action = "na.exclude"), data = .))
#   model_lag3 <- model_lag3[order(model_lag3$idObs), ]
# 
#   if (method == "std.resid")
#     {
#     is_after_trans <- which(all_dat$after_transitions)
# 
#     all_dat$lag1_out <- as.numeric(abs(model_lag1$.std.resid) > thr)
#     all_dat$lag1_out[is_after_trans] <- as.numeric(abs(model_lag1$.std.resid) > thr_after_transitions)[is_after_trans]
# 
#     all_dat$lag2_out <- as.numeric(abs(model_lag2$.std.resid) > thr)
#     all_dat$lag2_out[is_after_trans] <- as.numeric(abs(model_lag2$.std.resid) > thr_after_transitions)[is_after_trans]
#     all_dat$lag2_out[is_after_trans+1] <- as.numeric(abs(model_lag2$.std.resid) > thr_after_transitions)[is_after_trans+1]
# 
#     all_dat$lag3_out <- as.numeric(abs(model_lag3$.std.resid) > thr)
#     all_dat$lag3_out[is_after_trans] <- as.numeric(abs(model_lag3$.std.resid) > thr_after_transitions)[is_after_trans]
#     all_dat$lag3_out[is_after_trans+1] <- as.numeric(abs(model_lag3$.std.resid) > thr_after_transitions)[is_after_trans+1]
#     all_dat$lag3_out[is_after_trans+2] <- as.numeric(abs(model_lag3$.std.resid) > thr_after_transitions)[is_after_trans+2]
#     }
#   
#   else if (method == "IQR")
#     {
#     #require(rstatix)
#     # This part takes a few minutes
#     lag1_out <- model_lag1 %>%
#       group_by_at(c(treatment, "TimeFactor")) %>%
#       identify_outliers(".std.resid")
#     all_dat$lag1_out <- 0
#     all_dat$lag1_out[all_dat$idObs %in% lag1_out$idObs] <- 1
#     
#     lag2_out <- model_lag2 %>%
#       group_by_at(c(treatment, "TimeFactor")) %>%
#       identify_outliers(".std.resid")
#     all_dat$lag2_out <- 0
#     all_dat$lag2_out[all_dat$idObs %in% lag2_out$idObs] <- 1
#     
#     lag3_out <- model_lag3 %>%
#       group_by_at(c(treatment, "TimeFactor")) %>%
#       identify_outliers(".std.resid")
#     all_dat$lag3_out <- 0
#     all_dat$lag3_out[all_dat$idObs %in% lag3_out$idObs] <- 1
#     }
#  
#   else
#     {
#     stop("The 'method' argument should be either \"std.resid\" or \"IQR\".")
#     }
#   
#   all_dat$lag1_out[is.na(all_dat$diff1)] <- 1
#   all_dat$lag1_out_next <- c(all_dat$lag1_out[2:nrow(all_dat)], 1)
#   all_dat$lag2_out[is.na(all_dat$diff2)] <- 0.5
#   all_dat$lag2_out_next <- c(all_dat$lag2_out[3:nrow(all_dat)], 0.5, 0.5)
#   all_dat$lag3_out[is.na(all_dat$diff3)] <- 0.5
#   all_dat$lag3_out_next <- c(all_dat$lag3_out[4:nrow(all_dat)], 0.5, 0.5, 0.5)
#   all_dat$sum_out <- all_dat$lag1_out + all_dat$lag1_out_next + all_dat$lag2_out + all_dat$lag2_out_next + all_dat$lag3_out + all_dat$lag3_out_next
# 
#   return (all_dat)
#   }
# 
# 
# 
#------------------------------------------------------------------------------#
#                  Outliers based on day and night thresholds                  #
#------------------------------------------------------------------------------#
# 
#
# # The functions below detect outliers quite differently,
# # based on distinct thresholds that are set for day and night kinetics.
# #
# # They are normally used in conjunction with 'detect.outliers()' to ascertain that some points are outliers.
# # It is not recommended to use these functions as the sole criterion for detecting outliers because:
# #  - they rely on quite strong assumptions regarding the shape of the kinetics
# #  - not all cases have been implemented (e.g. 2+ consecutive outliers that are not around the transitions)
# #  - thresholds are unique for the whole dataset (while they should be adapted to each genotype/environment).
# 
# 
# # The function for the night
# find_outlier_night <- function (x, y, # the time and transpiration data
#                                 rel_threshold = 5, # relative threshold, compared to the normalized local slope: (local slope - avg) / avg
#                                 abs_threshold = 0.01/(0.5/24)) # absolute threshold for a transpiration difference between two consecutive points to be considered as unrealistic (because normalization is barely reliable when average slope is close to zero)
#   {
#   # Initialize vector of outliers (related to x and y)
#   out <- NULL
# 
#   # Remove NAs for outlier detection
#   NAs <- which(is.na(y))
#   ifelse (length(NAs) == 0, xout <- x, xout <- x[-NAs])
#   ifelse (length(NAs) == 0, yout <- y, yout <- y[-NAs])
# 
#   # Compute the average local slope by removing the start and the end (in case there are some outliers around)
#   avg <- mean((diff(yout)/diff(xout))[-(c(1:3, (length(yout)-2):length(yout)))])
# 
#   # Identify potential 'positive' (unrealistic opening) and 'negative' (unrealistic closure) outliers
#   pos <- which((diff(yout)/diff(xout)-avg)/avg > rel_threshold & abs(diff(yout)/diff(xout)) > abs_threshold)
#   neg <- which((diff(yout)/diff(xout)-avg)/avg < -rel_threshold & abs(diff(yout)/diff(xout)) > abs_threshold)
# 
#   # Initialize vector of outliers (related to xout and yout)
#   outdiff <- NULL
# 
#   # Run detection
#   if (length(pos) > 0)
#     {
#     # First, detect isolated outliers (cannot detect first point, second point if excessive closure, last point, or a group of 2+ consecutive outliers)
#     for (p in pos)
#       {
#       if (p == 1 | (p+1) %in% neg) { outdiff <- c(outdiff, p) }
#       else if ((p-1) %in% neg & (p-1) != 1) { outdiff <- c(outdiff, p-1) }
#       }
#     # Second, detect outliers showing excessive closure at the light-to-dark transition
#     if (pos[1] == 2)
#       {
#       if (length(pos) == 1)
#         {
#         outdiff <- c(outdiff, 2-1)
#         }
#       else
#         {
#         if (sum(diff(pos)) == length(pos)-1)
#           {
#           outdiff <- c(outdiff, pos-1)
#           }
#         else
#           {
#           outdiff <- c(outdiff, (2:pos[min(which(diff(pos) != 1))])-1)
#           }
#         }
#       }
#     # Third, detect outliers showing excessive reopening (or closure if avg < 0) at the dark-to-light transition
#     last_pos <- pos[length(pos)]
#     if (last_pos == length(yout) - 1)
#       {
#       if (length(pos) == 1)
#         {
#         outdiff <- c(outdiff, last_pos)
#         }
#       else
#         {
#         if (sum(diff(pos)) == length(pos)-1)
#           {
#           outdiff <- c(outdiff, pos)
#           }
#         else
#           {
#           first_consecutive_pos <- rev(pos)[min(which(rev(diff(pos)) != 1))]
#           outdiff <- c(outdiff, first_consecutive_pos:last_pos)
#           }
#         }
#       }
#     }
# 
#   if (length(neg) > 0)
#     {
#     # Fourth, detect outliers showing excessive closure (or reopening if avg < 0) at the dark-to-light transition
#     last_neg <- neg[length(neg)]
#     if (last_neg == length(yout) - 1)
#       {
#       if (length(neg) == 1)
#         {
#         outdiff <- c(outdiff, last_neg)
#         }
#       else
#         {
#         if (sum(diff(neg)) == length(neg)-1)
#           {
#           outdiff <- c(outdiff, neg)
#           }
#         else
#           {
#           first_consecutive_neg <- rev(neg)[min(which(rev(diff(neg)) != 1))]
#           outdiff <- c(outdiff, first_consecutive_neg:last_neg)
#           }
#         }
#       }
#     }
# 
#   # Assign the outliers to the correct vector and visualise
#   if (!is.null(outdiff))
#     {
#     if (length(NAs) == 0)
#       {
#       out <- outdiff + 1
#       }
#     else
#       {
#       for (o in 1:length(outdiff))
#         {
#         nb_NA_before <- length(NAs <= outdiff[o]+1)
#         out <- outdiff + 1 + nb_NA_before
#         }
#       }
#     }
# 
#   # Returns the indices
#   return (out)
#   }
# 
# 
# # The function for the day
# find_outlier_day <- function (x, y, # the time and transpiration data
#                               abs_threshold = 0.05/(0.5/24)) # absolute threshold for a transpiration difference between two consecutive points to be considered as unrealistic if the next pair of points shows the opposite behaviour
#   {
#   # Initialize vector of outliers (related to x and y)
#   out <- NULL
# 
#   # Remove NAs for outlier detection
#   NAs <- which(is.na(y))
#   ifelse (length(NAs) == 0, xout <- x, xout <- x[-NAs])
#   ifelse (length(NAs) == 0, yout <- y, yout <- y[-NAs])
# 
#   # Identify potential 'positive' (unrealistic opening) and 'negative' (unrealistic closure) outliers
#   pos <- which(diff(yout)/diff(xout) > abs_threshold)
#   neg <- which(diff(yout)/diff(xout) < -abs_threshold)
#   neutral <- which(diff(yout)/diff(xout) >= -abs_threshold & diff(yout)/diff(xout) <= abs_threshold)
# 
#   # Initialize vector of outliers (related to xout and yout)
#   outdiff <- NULL
# 
#   # Initialize check for last point (should not be considered as pos outlier if the point before is neg)
#   checklast <- TRUE
# 
#   # Run detection
#   if (length(neg) > 0)
#     {
#     # First, detect isolated outliers (cannot detect first point, last point, or a group of 2+ consecutive outliers)
#     for (n in neg)
#       {
#       if (n == 1 | (n+1) %in% pos) { outdiff <- c(outdiff, n) }
#       else if ((n-1) %in% pos) { outdiff <- c(outdiff, n-1) }
#       }
#     if (neg[length(neg)] == length(yout) - 2)
#       {
#       checklast <- FALSE
#       }
#     }
# 
#   if (checklast & length(pos) > 0)
#     {
#     # Second, detect outliers showing excessive reopening at the light-to-dark transition
#     last_pos <- pos[length(pos)]
#     if (last_pos == length(yout) - 1)
#       {
#       if (length(pos) == 1)
#         {
#         outdiff <- c(outdiff, last_pos)
#         }
#       else
#         {
#         if (sum(diff(pos)) == length(pos)-1)
#           {
#           outdiff <- c(outdiff, pos)
#           }
#         else
#           {
#           first_consecutive_pos <- rev(pos)[min(which(rev(diff(pos)) != 1))]
#           outdiff <- c(outdiff, first_consecutive_pos:last_pos)
#           }
#         }
#       }
#     }
# 
#   # Assign the outliers to the correct vector and visualise
#   if (!is.null(outdiff))
#     {
#     if (length(NAs) == 0)
#       {
#       out <- outdiff + 1
#       }
#     else
#       {
#       for (o in 1:length(outdiff))
#         {
#         nb_NA_before <- length(NAs <= outdiff[o]+1)
#         out <- outdiff + 1 + nb_NA_before
#         }
#       }
#     }
# 
#   # Returns the indices
#   return (out)
#   }
# 
# 
# # The function to run day and night detection on the data
# day_night_outliers <- function (all_dat, # the whole dataframe
#                                 Time_var = if ("decimalDay" %in% names(all_dat)) "decimalDay" else "Time", # the column name of the time axis
#                                 E_var = "E_clean", # the column name of the transpiration axis
#                                 time_start = min(all_dat[, Time_var], na.rm = T), # time at which outlier detection starts; should match a time at which light is switched on (start of the day)
#                                 time_end = max(all_dat[, Time_var], na.rm = T), # time at which outlier detection ends
#                                 night_rel_threshold = 5,
#                                 night_abs_threshold = 0.01/(0.5/24),
#                                 day_abs_threshold = 0.05/(0.5/24))
#   {
#   all_dat$outlier_DAY_NIGHT <- F
#   for (pot in unique(all_dat$idPot))
#     {
#     t2 <- time_start
#     dat <- all_dat[all_dat$idPot == pot, ]
#     while (t2 < time_end)
#       {
#       # Select the diel period
#       t1 <- t2
#       t2 <- t1 + Pho_Per + Sko_Per
#       
#       # Day
#       day <- dat[, Time_var] >= t1 & dat[, Time_var] <= t1 + Pho_Per
#       x_day <- dat[day, Time_var]
#       y_day <- dat[day, E_var]
#       if (length(which(!is.na(y_day))) > 1)
#         {
#         out_day <- find_outlier_day(x_day, y_day, day_abs_threshold)
#         if (!is.null(out_day))
#           {
#           all_dat$outlier_DAY_NIGHT[dat$idObs[day][out_day]] <- T
#           }
#         }
#       
#       # Night
#       night <- dat[, Time_var] >= t1 + Pho_Per & dat[, Time_var] <= t2
#       x_night <- dat[night, Time_var]
#       y_night <- dat[night, E_var]
#       if (length(which(!is.na(y_night))) > 1)
#         {
#         out_night <- find_outlier_night(x_night, y_night, night_rel_threshold, night_abs_threshold)
#         if (!is.null(out_night))
#           {
#           all_dat$outlier_DAY_NIGHT[dat$idObs[night][out_night]] <- T
#           }
#         }
#       }
#     }
#   return (all_dat)
#   }
#
#
#------------------------------------------------------------------------------#
#             Calling the functions at Step#02 (example from C2M47)            #
#------------------------------------------------------------------------------#
#
# # Identify transitions (required for method "std.resid")
# df <- identify_transitions(df, Time_var = "decimalDay", light_ON = Time_ON0 - floor(Time_ON0))
# 
# # The method "std.resid" (takes less than 1 min) requires defining two thresholds to select outliers
# # ('thr' for most points, and 'thr_after_transitions' for being more permissive right after the day/night transitions).
# # See 'PhenoLeaks_outliers.R' to read more about the 'detect_outliers()' function.
# df <- detect_outliers(df,
#                       treatment = c("idGenotype", "idWatering"),
#                       Time_var = "decimalDay",
#                       E_var = "E_clean",
#                       method = "std.resid",
#                       thr = 1, thr_after_transitions = 3)
# 
# # Alternatively, one can use the  method "IQR" (takes a few minutes)
# # but it looks less performant in the Arabidopsis dataset:
# #df <- detect_outliers(df,
# #                     treatment = c("idGenotype", "idWatering"),
# #                     Time_var = "decimalDay",
# #                     E_var = "E_clean",
# #                     method = "IQR")
# 
# # Additional outlier detection based on separate day and night thresholds
# # for the difference between two consecutive points (takes a few seconds).
# # See 'PhenoLeaks_outliers.R' to read more about the 'day_night_outliers()' function.
# df <- day_night_outliers(df,
#                          Time_var = "decimalDay",
#                          E_var = "E_clean",
#                          time_start = Time_ON0 - 1,
#                          time_end = Time_end_exp,
#                          night_rel_threshold = 5,
#                          night_abs_threshold = 0.01/(0.5/24),
#                          day_abs_threshold = 0.05/(0.5/24))
# 
# 
# # Visualize individual curves with outliers
# pdf(file.path(getwd(), "Figures", "Step#02c_visualize_potential_outliers_all_pots.pdf"), width = 10, height = 5)
# for (irr in c("WW", "WS"))
#   {
#   for (geno in sort(unique(df$idGenotype[df$idWatering == irr])))
#     {
#     for (pot in sort(unique(df$idPot[df$idWatering == irr & df$idGenotype == geno])))
#       {
#       # Prepare the plot
#       dat <- df[df$idPot == pot, ]
#       prepare_kin(dat, Time_var = "decimalDay", E_var = "E_mmol_per_m2_s",
#                   main = paste(geno, irr, paste("Pot", pot, sep = " "), sep = " - "),
#                   light_ON = Time_ON0 - floor(Time_ON0),
#                   inside = F)
#          
#       # Print data
#       points(E_clean ~ decimalDay, data = dat, pch = 16, cex = 0.5)
#       # which - at this point - is the same as:
#       #points(E_mmol_per_m2_s ~ decimalDay, data = dat[!dat$outlier, ] pch = 16, cex = 0.5)
#       
#       # Print outliers identified manually (whole periods)
#       points(E_mmol_per_m2_s ~ decimalDay, data = dat[dat$outlier, ], pch = 1, cex = 0.5, col = "grey")
#       
#       # Check individual penalties
#       #points(E_mmol_per_m2_s ~ decimalDay, data = dat[!dat$outlier & dat$lag3_out, ], pch = 16, col = "darkolivegreen2", cex = 1.4)
#       #points(E_mmol_per_m2_s ~ decimalDay, data = dat[!dat$outlier & dat$lag3_out_next, ], pch = 21, col = "darkolivegreen4", cex = 1.5)    
#       #points(E_mmol_per_m2_s ~ decimalDay, data = dat[!dat$outlier & dat$lag2_out, ], pch = 16, col = "aquamarine2", cex = 0.9)
#       #points(E_mmol_per_m2_s ~ decimalDay, data = dat[!dat$outlier & dat$lag2_out_next, ], pch = 21, col = "aquamarine4", cex = 1)
#       #points(E_mmol_per_m2_s ~ decimalDay, data = dat[!dat$outlier & dat$lag1_out, ], pch = 16, col = "coral2", cex = 0.5)
#       #points(E_mmol_per_m2_s ~ decimalDay, data = dat[!dat$outlier & dat$lag1_out_next, ], pch = 21, col = "coral4", cex = 0.6)
#       
#       # Visualize different thresholds for the sum of penalties
#       points(E_mmol_per_m2_s ~ decimalDay, data = dat[!dat$outlier & dat$sum_out >= 3, ], pch = 16, col = 4, cex = 0.5)
#       points(E_mmol_per_m2_s ~ decimalDay, data = dat[!dat$outlier & dat$sum_out >= 4, ], pch = 16, col = 3, cex = 0.5)
#       points(E_mmol_per_m2_s ~ decimalDay, data = dat[!dat$outlier & dat$sum_out >= 5, ], pch = 16, col = 2, cex = 0.5)
# 
#       # Visualize the points identify as suspicious by the fixed day/night thresholds
#       points(E_mmol_per_m2_s ~ decimalDay, data = dat[!dat$outlier & dat$outlier_DAY_NIGHT, ], pch = 21, col = 5, cex = 0.9)
#       
#       # Add the legend
#       legend("topright", ncol = 4, pch = 1, col = c(NA, "grey", NA, NA), cex = 0.5, bty = "n",
#              legend = c("VISUALLY:", "Outlier", "", "                    "))
#       legend("topright", inset = c(0, 0.025), ncol = 4, pch = 16, col = c(NA, 4:2), cex = 0.5, bty = "n",
#              legend = c("STD.RESID:", "3 penalties", "4 penalties", "5+ penalties"))
#       legend("topright", inset = c(0, 0.05), ncol = 4, pch = 21, col = c(NA, 5, NA, NA), cex = 0.5, pt.cex = 0.9, bty = "n",
#              legend = c("DAY/NIGHT:", "diff > thr", "", ""))
#             
#       rm(dat)
#       }
#     }
#   }
# 
# dev.off()
# 
# 
# # Look at the kinetics and define the rule for setting points as outliers
# df$rule_out <- df$sum_out >= 5 | (df$sum_out >= 3 & df$outlier_DAY_NIGHT)
#
#~~~~~~~~~~~~~~~~~~~        END OF OLD IMPLEMENTATION       ~~~~~~~~~~~~~~~~~~~#
