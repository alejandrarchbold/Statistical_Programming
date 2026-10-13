# Statistical Programming
# Practical 1: Sequence duels and perils of independence

# Group 3: Alejandra Campo Archbold (s3023700), Weihang Gu (s2932201), Weichen Zhao (s2929481)

# Contributions:

# - Alejandra Campo Archbold:
# - Weihang Gu: 
# - Weichen Zhao:

# GitHub URL project: https://github.com/alejandrarchbold/Statistical_Programming.git

#-------------------------------------------------------------------------------
# Overview:

#
#
#

#-------------------------------------------------------------------------------


#setwd("put/your/local/repo/location/here")

firsthit <- function(x, a) {
  # Description:
  # Input: x is a vector of 0s and 1s where 1=, a pattern
  
  k <- length(a)
  n <- length(x)
  
  if (k == 0) {
    stop("Add at least one element.")
  }
  
  if (n < k) {
    return(NA)
  }
  
  for (m in k:n) {
    w <- x[(m - k + 1):m]
    
    if (all(w == a)) {
      return(m)
    }
  }
  return(NA)
}
