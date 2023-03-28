################################################################################
#                                                                              #
#                             PhenoLeaks - GENERIC                             #
#                                                                              #
#                          Generic functions required                          #
#                       to run PhenoLeaks core functions                       #
#                                                                              #
#                             Florent Pantin, 2022                             #
#                                                                              #
################################################################################




#------------------------------------------------------------------------------#
#                       Operator for multiple assignments                      #
#------------------------------------------------------------------------------#


# Function to define an operator ':=' allowing us to retrieve more than one value
# from a function output without having to use 'unlist'.
# See http://code.google.com/p/miscell/source/browse/rvalues/rvalues.r
':=' = function(lhs, rhs)
  {
  frame = parent.frame()
  lhs = as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs = lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs = list(rhs)
  if (length(lhs) > length(rhs))
    rhs = c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL))
  }



#------------------------------------------------------------------------------#
#                               Manage constants                               #
#------------------------------------------------------------------------------#

# Lock the constants of an experiment in the global R environment
# and return them as a named list (for console display)
lock_constants <- function (e) # the environment of the caller function, typically passed as 'environment()'
  {
  cst <- ls(e)
  constants <- list()
  for (i in 1:length(cst))
    {
    constants[[cst[i]]] <- get(cst[i], e)
    if (exists(cst[i], globalenv())) unlockBinding(cst[i], globalenv())
    assign(cst[i], constants[[cst[i]]], envir = globalenv())
    lockBinding(cst[i], globalenv())
    }
  return (constants)
  }

