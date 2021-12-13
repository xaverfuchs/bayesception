# Code referring to the two-tap model with equal variances
# Author: Dr. Xaver Fuchs, xaver.fuchs@plus.ac.at
# Formulae based on: Goldreich & Tong, 2013, Frontiers in Psychology

##formula A08
computePrior <- function(x1, x2, sigma_v, time_t) {
  Prior <- 1/(sqrt(2*pi)*sigma_v*time_t) * exp(-(x2-x1)^2/(2*(sigma_v*time_t)^2))
  return(Prior)
}

##formula A11
computeLikelihood <- function(x1m, x2m, x1, x2, sigma_s) {
  p_x1m_given_x1 <- 1/(sqrt(2*pi)*sigma_s) * exp(-(x1m-x1)^2/(2*sigma_s^2))
  p_x2m_given_x2 <- 1/(sqrt(2*pi)*sigma_s) * exp(-(x2m-x2)^2/(2*sigma_s^2))
  Likelihood <- p_x1m_given_x1*p_x2m_given_x2
  return(Likelihood)
}

##formula A13
computePosterior <- function(x1m, x2m, x1, x2, sigma_s, sigma_v, time_t) {
  Posterior <- exp(-(((x1m-x1)^2+(x2m-x2)^2)/(2*sigma_s^2) + (x2-x1)^2/(2*(sigma_v*time_t)^2)))
  return(Posterior)
}


##formula A14
compute_x1_star <- function(x1m, x2m, sigma_s, sigma_v, time_t) {
  x1_star <- x1m * (((sigma_v*time_t)^2+sigma_s^2) / ((sigma_v*time_t)^2+2*sigma_s^2)) + x2m * (sigma_s^2 / ((sigma_v*time_t)^2+2*sigma_s^2))
  
  return(x1_star)
}

compute_x2_star <- function(x1m, x2m, sigma_s, sigma_v, time_t) {
  x2_star <- x1m * (sigma_s^2 / ((sigma_v*time_t)^2+2*sigma_s^2)) + x2m * (((sigma_v*time_t)^2+sigma_s^2) / ((sigma_v*time_t)^2+2*sigma_s^2)) 
  
  return(x2_star)
}

compute_common_sigma_square <- function(sigma_s, sigma_v, time_t) {
  sigma_square <- sigma_s^2 * (sigma_s^2 + (sigma_v*time_t)^2) / (2*sigma_s^2 + (sigma_v*time_t)^2)
  return(sigma_square)
}

compute_correlation <- function(sigma_s, sigma_v, time_t) {
  distr_correlation <- sigma_s^2 / (sigma_s^2 + (sigma_v*time_t)^2)
  return(distr_correlation)
}


computePosterior2 <- function(x1m, x2m, x1, x2, sigma_s, sigma_v, time_t, returnList=F) {
  #first compute the necessary variables
  x1_star <- compute_x1_star(x1m = x1m, x2m = x2m, sigma_s = sigma_s, sigma_v = sigma_v, time_t = time_t)
  x2_star <- compute_x2_star(x1m = x1m, x2m = x2m, sigma_s = sigma_s, sigma_v = sigma_v, time_t = time_t)
  common_sigma_square <- compute_common_sigma_square(sigma_s = sigma_s, sigma_v = sigma_v, time_t = time_t)
  correlation <- compute_correlation (sigma_s = sigma_s, sigma_v = sigma_v, time_t = time_t)
  
  #now use formula 14 to compute posterior 
  Posterior <- exp(  ( -  1/(2*(1-correlation^2)) * 
                         ( ((x1-x1_star)^2) + (x2-x2_star)^2 - 2*correlation * 
                             (x1-x1_star)*(x2-x2_star) )/common_sigma_square ) )
  if (returnList==T) {
    #extend to be a list
    returnList <- list( x1_star= x1_star, x2_star=x2_star, 
                        common_sigma_square=common_sigma_square, correlation=correlation, 
                        Posterior = Posterior )
    return(returnList)
    
  } else {
    return(Posterior)
  }
}
