#Diana Gerardo
#AMS 274 Homework 3
#Problem 2 part c

library(MCMCpack)
library(MASS)
library(mvtnorm)
library(readxl)

options(scipen = 9999)

# the data
beetle_dat <- read_excel("~/Generalized_Linear_Models/HW3/beetle_dat.xlsx")

y <- beetle_dat$num.killed 
m <- beetle_dat$beetles

x <- beetle_dat$logdose # little x, feature
X <- as.matrix(x)
X <- cbind(rep(1, nrow(X)), X) # big X, covariates

# The parametric link funcition given for part c is not a
# common link function built in glm(), so I manually retrieve the 
# fitted values and variance-covariance matrix to perform M-H

# Functions-----------------------------------------------------------------

modified_logit_link <- function(covariates, beta, alpha){
  eta <- covariates%*%beta
  value <- rep(0,length(covariates[1,]))
  
  for (i in 1:length(covariates[1,])) {
    value[i] <- exp(alpha*eta[i]) / (1 + exp(eta[i]))^alpha
  }
  return(value)
}

modlog_func <- function(x,b,a){
  eta <- b[1]+b[2]*x
  numer <- exp(eta)
  denom <- 1 + exp(eta)
  val <- (numer/denom)^a
}

modlog_func2 <- function(x,b1,b2,a){
  eta <- b1+b2*x
  numer <- exp(eta)
  denom <- 1 + exp(eta)
  val <- (numer/denom)^a 
}

modlogit_devres <- function(covariates, obs_denom, obs_numer, beta_hat, alpha_hat){
  p_hat <- modified_logit_link(covariates, beta_hat, alpha_hat)
  m <- obs_denom
  y_hat <- p_hat*m
  y_obs <- obs_numer
  n <- length(y_obs)
  
  term1 <- 2*abs(y_obs*log(y_obs/y_hat))
  term2 <- (m - y_obs)*log(1E-10+(m - y_obs)/(m - y_hat))
  deviance <- sign(y_obs - y_hat) * sqrt(term1 + term2)
  
  return(data.frame(p_hat, deviance))
}

# --------------------------------------------------------------------------

# given values
beta_hat <- c(-113.625, 62.5)
alpha_hat <- .279

modlogit_dat <- modlogit_devres(X, m, y, beta_hat, alpha_hat)
fitted_values <- modlogit_dat$p_hat

# obtaining vcov manually
W <- beetle_dat$beetles*fitted_values*(1 - fitted_values)
W <- diag(as.numeric(W))
infMatrix <- t(X) %*% W %*% X 
vcov_by_hand <- solve(infMatrix)

# Metropolis-Hastings Algorithm---------------------------------------

prior.beta <- function(x){return(log(1))}
prior.alpha <- function(x){dgamma(x,2,3)}

posterior <- function(beta,alpha){
  priors <- prior.beta(beta) + prior.alpha(alpha)
  term1 <- y*log(modlog_func(x,beta,alpha))
  term2 <- (m-y)*log(1 - modlog_func(x,beta,alpha))
  sum <- sum(term1+term2)
  value <- exp(priors + sum)
  
  return(value)
  
}

#initialize 
iters <- 425000
iter.beta <- matrix(NA, nrow=2, ncol=iters)
iter.beta[,1] <- c(-113.625, 62.5)
iter.alpha <- matrix(NA, nrow=1, ncol=iters)
iter.alpha[1] <- .279

new.beta <- matrix(NA, nrow=2, ncol=iters)
new.alpha <- matrix(NA, nrow=1, ncol=iters)
 

J <- vcov_by_hand
d <- 43 #tuning parameter

set.seed(88800)

#MCMC
for(i in 2:iters)
{
  new.beta[,i] <- rmvnorm(1, iter.beta[,i-1], d*J)
  new.alpha[i] <- rgamma(1, iter.alpha[i-1], 1)
  
  numer <- posterior(new.beta[,i],new.alpha[i])
  denom <- posterior(iter.beta[,i-1], iter.alpha[i-1])
  
  q <- min(1, numer/denom)
  
  while (q != 'NaN') {
    if(runif(1) < q){
      iter.beta[,i] <- new.beta[,i]
      iter.alpha[i] <- new.alpha[i]
      break
    }
    
    else{
      iter.beta[,i] <- iter.beta[,i-1]
      iter.alpha[i] <- iter.alpha[i-1]
      break
    }
  }
  
  if (q == 'NaN'){
    iter.beta[,i] <- iter.beta[,i-1]
    iter.alpha[i] <- iter.alpha[i-1]
  }
  
}

# -----------------------------------------------------------------------

## rename simulated samples
b1.new <- iter.beta[1,]
b2.new <- iter.beta[2,]
alpha.new <- iter.alpha

# ## Discarding burn-in samples
# b1.burnin = 100
# b2.burnin = 100
# alpha.burnin = 100
# 
# b1.new <- b1.new[-(1:b1.burnin)]
# b2.new <- b2.new[-(1:b2.burnin)]
# alpha.new<- alpha.new[-(1:alpha.burnin)]

## These tend to change from sample to sample because of randomization
b1.thin = 133
b2.thin = 133
alpha.thin = 133

## Thinning the chains
b1.new <- b1.new[seq(1, length(b1.new), b1.thin)]
b2.new <- b2.new[seq(1, length(b2.new), b2.thin)]
alpha.new <- alpha.new[seq(1, length(alpha.new), alpha.thin)]

## trace plots of beta0, beta1, and alpha 
plot.ts(b1.new, ylab=expression(beta[1]))
plot.ts(b2.new, ylab=expression(beta[2]))
plot.ts(as.numeric(alpha.new), ylab=expression(alpha))

## mle estimates compared to simulation mean estimates
iter.beta[,1]; mean(b1.new); mean(b2.new) 
iter.alpha[1]; mean(alpha.new)

## Density plots of estimated parameters
plot(density(b1.new), main=expression(beta[1]),xlab=expression(beta[1]))
polygon(density(b1.new), col="coral3", border="brown")

plot(density(b2.new), main=expression(beta[2]), xlab=expression(beta[2]))
polygon(density(b2.new), col="coral3", border="brown")

plot(density(alpha.new), main=expression(alpha), xlab=expression(alpha))
polygon(density(alpha.new), col="coral3", border="brown")

# ----------------------------------------------------------------------------

## median lethal dose function
LD_50 <- function(b1.new,b2.new,alpha.new){
  term1 <- exp(b1.new)*2^(1/alpha.new)
  term2 <- exp(b1.new)
  the_log <- -1*log(term1-term2)
  value <- the_log/b2.new
  
  return(value)
}

## median lethal does density
LD_values <- LD_50(b1.new, b2.new, alpha.new)
plot(density(LD_values), main=expression("LD"[50]))
polygon(density(LD_values), col="coral3", border="brown")
mean(LD_values)

## dose-repsonse curve
x.grid <- seq(1.65, 1.9, length=100)
quantile <- matrix(NA,100,3)

for(i in 1:100){
  new.pi <- modlog_func2(x.grid[i],b1.new, b2.new, alpha.new)
  quantile[i,] <- quantile(new.pi,c(0.025,0.5,0.975))
}

plot(x.grid,quantile[,1], type = "l", lty=2, lwd=2,
     xlab="log dose", ylab="# killed / #beetles")
lines(x.grid,quantile[,2], type = "l", lwd=2, col="red")
lines(x.grid, quantile[,3], type="l", lty=2, lwd=2)
points((y/m)~x,pch=19, col="grey2") # little x, without intercept column

# -----------------------------------------------------------------------------

## Residuals vs Fitted values plot
## note use little x, without intercept column 
mean_param <- c(mean(b1.new), mean(b2.new), mean(alpha.new))
residuals.new <- (y/m) - modlog_func2(x,mean_param[1],mean_param[2],mean_param[3])
plot(residuals.new ~ modlog_func2(x,mean_param[1],mean_param[2], mean_param[3]),
     ylim=c(-0.3,0.3), ylab="residuals", xlab="fitted values",
     main="ModLogit Deviance Residuals")
abline(h=0)
abline(h=-0.2, col="red")
abline(h=0.2, col="red")

## Drawing from the Posterior Predictive Distribution
predictive <- function(x, b1, b2, a, weights){
  num_obs <- length(weights)
  iters <- length(b1)
  predictive_draws <- matrix(0, nrow=iters, ncol=num_obs)
  
  for(i in 1:iters){
    pi <- modlog_func2(x,b1[i],b2[i],a[i])
    
    for (j in 1:num_obs) {
      predictive_draws[i,j] <- rbinom(n = 1, size = weights[j], prob = pi[j])
      
    }
  }
  return(predictive_draws)
}

## each row is a predictive distribution
pred_draws <- predictive(x, b1.new, b2.new, alpha.new, m)
dim(pred_draws)
head(pred_draws)

## Gelfand and Gosh Criterion (posterior predictive quad loss)
gg <- function(obs_var, preds){
  y <- obs_var
  
  means <- apply(preds, 2, mean)
  func_of_mean <- (y-means)^2
  
  variance <- apply(preds, 2, var)
  
  summation <- sum(func_of_mean) + sum(variance)
  
  return(summation)
}

gg(y, pred_draws)
# [1] 116.1154
