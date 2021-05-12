############################ Generalized Linear Models ###########################
############################ Homework 3 Problem 3 Code ###########################

# By: Diana Gerardo

# Libraries ----------------------------------------------------------------------
library(MASS)

# Problem 3 ----------------------------------------------------------------------

# Consider the data set from: 
# http://www.stat.columbia.edu/~gelman/book/data/fabric.asc
# on the incidence of faults in the manufacturing of rolls of fabric.
# 1st column contains the length of each roll, 2nd contains the number of faults.

fabric <- read.csv("~/Generalized_Linear_Models/HW2/fabric.csv")
head(fabric)

# Part (a) -----------------------------------------------------------------------

# Fit a Poisson GLM, with logarithmic link.
# glm()'s Poisson default link is the logarithmic
fabric_glm_pois <- glm(formula = faults ~ length,
                       family = "poisson", data = fabric)

# Increasing the fabric length by 1 unit 
# multiplies the avg number of faults by
exp(fabric_glm_pois$coefficients[2])

# The expected fabric length of faults 
# found in 32 rolls of fabric is
exp(fabric_glm_pois$coefficients[1] + fabric_glm_pois$coefficients[2])

summary(fabric_glm_pois)
# From the summary, the coefficients are highly significant.
# Dispersion parameter for the poisson family taken to be 1.
# We also see that the residual deviance is greater than the 
# degrees of freedom, so that we have over-dispersion. 
# This means that there is extra variance not accounted for 
# by the model or by the error structure.

# Part (b) -----------------------------------------------------------------------

# Fit a QuasiPoisson GLM, with logarithmic link
fabric_glm_quasi <- glm(formula = faults ~ length,
                        family = "quasipoisson", data = fabric)

summary(fabric_glm_quasi)
# Estimates are the same as previous.
# Dispersion parameter taken to be 2.121965

# Part (c) -----------------------------------------------------------------------

point_estimate <- function(x,betas){
  Xb <- t(x)%*%betas
  
  return(Xb)
}

confidence_interval <- function(x, betas, inv_fisher){
  # fixed at the 95% confidence interval
  Xb <- t(x)%*%betas
  XJX <- t(x)%*%inv_fisher%*%x
  low <- Xb - 1.96*sqrt(XJX)
  up <- Xb + 1.96*sqrt(XJX)
  interval <- c(low,up)
  
  return(interval)
}

# point and interval estimates at
x <- rbind(1,500) 
xo <-rbind(1,995)


########### Poisson GLM ########### 
J.inv_pois <- vcov(fabric_glm_pois) 
betas_pois <- fabric_glm_pois$coefficients

pois_point_x <- point_estimate(x,betas_pois)
pois_point_xo <- point_estimate(xo,betas_pois)

pois_interval_x <- confidence_interval(x, betas_pois, J.inv_pois)
pois_interval_xo <- confidence_interval(xo,betas_pois, J.inv_pois)

######## Quasi Poisson GLM ########
J.inv_quasi <- vcov(fabric_glm_quasi) 
betas_quasi <- fabric_glm_quasi$coefficients

quasi_point_x <- point_estimate(x,betas_quasi)
quasi_point_xo <- point_estimate(xo,betas_quasi)

quasi_interval_x <- confidence_interval(x, betas_quasi, J.inv_quasi)
quasi_interval_xo <- confidence_interval(xo,betas_quasi, J.inv_quasi)
