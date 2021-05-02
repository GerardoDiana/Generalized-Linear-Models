#Diana Gerardo
#AMS 274 Homework 3
#Problem 2

library(MCMCpack)
library(MASS)
library(mvtnorm)

library(readxl)
beetle_dat <- read_excel("~/Generalized_Linear_Models/HW3/beetle_dat.xlsx")
# View(beetle_dat)


########################################### Part (a) ################################################
y <- beetle_dat$num.killed 
m <- beetle_dat$beetles
x <- beetle_dat$logdose

# Complementary log-log (CLL) Bin GLM
CLL <- glm((num.killed/beetles) ~ logdose, weights = beetles, family = binomial(link="cloglog"), 
           data = beetle_dat)
 

#initialize 
iter.beta <- matrix(NA, nrow=2, ncol=10000)
iter.beta[,1] <- c(0,0)

#flat prior log(1), normal prior dnorm dnorm(x[1],0,1,log = T)+dnorm(x[2],0,1,log=T)
prior <- function(x){
  return(log(1))
}

posterior <- function(beta)
{
  return(exp(prior(beta) + sum( y*log(1-exp(-exp(beta[1]+beta[2]*x))) + (m-y)*(-exp(beta[1]+beta[2]*x))) ))
}

new.beta <- matrix(NA, nrow=2, ncol=10000)
new.beta[,1] <- c(0,0) 

J <- vcov(CLL)
d <- 2 #tuning parameter

#MCMC
for(i in 2:10000)
{
  new.beta[,i] <- rmvnorm(1, iter.beta[,i-1], d*J)
  
  q <- min(1, posterior(new.beta[,i])/posterior(iter.beta[,i-1]))
  
  if(runif(1) < q)
  {iter.beta[,i] <- new.beta[,i]}
  
  else
  {iter.beta[,i] <- iter.beta[,i-1]}
   
}

#Remove first 1000 iters
T <- 1000
iter.beta<- iter.beta[,-(1:T)] 


#trace plots of beta1 and beta2
plot.ts(iter.beta[1,], ylab=expression(beta[1]))
plot.ts(iter.beta[2,], ylab=expression(beta[2]))
CLL
mean(iter.beta[2,])
plot(density(iter.beta[1,]), main=expression(beta[1]),
     xlab=expression(beta[1]))
polygon(density(iter.beta[1,]), col="thistle", border="purple")

plot(density(iter.beta[2,]), main=expression(beta[2]), 
     xlab=expression(beta[2]))
polygon(density(iter.beta[2,]), col="thistle", border="purple")


#LD50 values and density plot
LD50 <- (-0.366513-iter.beta[1,])/iter.beta[2,]
plot(density(LD50), main=expression("LD"[50]))
polygon(density(LD50), col="thistle", border="purple")

#evaluating pi(x_i) for all x.grid
pi <- function(x,beta)
{return(1-exp(-exp(beta[1,]+beta[2,]*x)))}

x.grid <- seq(1.65, 1.9, length=100)
quantile <- matrix(NA,100,3)
for(i in 1:100){
  new.pi <- pi(x.grid[i],iter.beta)
  quantile[i,] <- quantile(new.pi,c(0.025,0.5,0.975))
}

# dose-response graph
plot(x.grid,quantile[,1], type = "l", lty=2, xlab="log dose", ylab="# killed/ # beetles", lwd=2)
lines(x.grid,quantile[,2], type = "l", col="purple", lwd=2)
lines(x.grid, quantile[,3], type="l", lty=2, lwd=2)
points((beetle_dat$num.killed/beetle_dat$beetles)~beetle_dat$logdose, pch=19, col="grey2")

pi.e <- function(x,b1,b2)
{return(1-exp(-exp(b1+b2*x)))}

cll.res <- (y/m)-pi.e(x,-40,22)
plot(cll.res~pi.e(x,-40,22), ylim=c(-0.25,0.25),
     ylab="residuals", xlab="fitted values", main="CLL Deviance Residuals")
abline(h=0)
abline(h=0.2, col="red")
abline(h=-0.2, col="red")
## quad loss
sum(cll.res^2)/8

########################################### Part (b) ################################################

y <- beetle_dat$num.killed 
m <- beetle_dat$beetles
x <- beetle_dat$logdose

# Complementary log-log (CLL) Bin GLM
logit <- glm((num.killed/beetles) ~ logdose, weights = beetles, family = binomial(link="logit"), 
           data = beetle_dat)
 

#initialize 
iter.beta <- matrix(NA, nrow=2, ncol=10000)
iter.beta[,1] <- c(0,0)

#flat prior log(1)
prior <- function(x){
  return(log(1))
}

posterior <- function(beta)
{
  return(exp(prior(beta) + sum( y*log(exp(beta[1]+beta[2]*x)/(1+exp(beta[1]+beta[2]*x))) + 
                                  (m-y)*log(1-exp(beta[1]+beta[2]*x)/(1+exp(beta[1]+beta[2]*x)))) ))
}

new.beta <- matrix(NA, nrow=2, ncol=10000)
new.beta[,1] <- c(0,0) 

J <- vcov(logit)
d <- 2 #tuning parameter

#MCMC
for(i in 1:10000)
{
  new.beta[,i] <- rmvnorm(1, iter.beta[,i], d*J)
  
  q <- min(1, posterior(new.beta[,i])/posterior(iter.beta[,i]))
  
  if(runif(1) < q)
  {iter.beta[,i+1] <- new.beta[,i]}
  
  else
  {iter.beta[,i+1] <- iter.beta[,i]}
   
}

#Remove first 1000 iters
T <- 1000
iter.beta<- iter.beta[,-(1:T)] 


#trace plots of beta1 and beta2
plot.ts(iter.beta[1,], ylab=expression(beta[1]))
plot.ts(iter.beta[2,], ylab=expression(beta[2]))
logit

plot(density(iter.beta[1,]), main=expression(beta[1]),
     xlab=expression(beta[1]))
polygon(density(iter.beta[1,]), col="lightblue", border="blue")

plot(density(iter.beta[2,]), main=expression(beta[2]), 
     xlab=expression(beta[2]))
polygon(density(iter.beta[2,]), col="lightblue", border="blue")

#LD50 values and density plot
LD50 <- (-iter.beta[1,])/iter.beta[2,]
plot(density(LD50), main=expression("LD"[50]))
polygon(density(LD50), col="lightblue", border="blue")

#evaluating pi(x_i) for all x.grid
pi <- function(x,beta)
{return(exp(beta[1,]+beta[2,]*x)/(1+exp(beta[1,]+beta[2,]*x)))}

x.grid <- seq(1.65, 1.9, length=100)
quantile <- matrix(NA,100,3)
for(i in 1:100){
  new.pi <- pi(x.grid[i],iter.beta)
  quantile[i,] <- quantile(new.pi,c(0.025,0.5,0.975))
}

# dose-response graph
plot(x.grid,quantile[,1], type = "l", lty=2, lwd=2, xlab="log dose", ylab="# killed / #beetles")
lines(x.grid,quantile[,2], type = "l", lwd=2, col="blue")
lines(x.grid, quantile[,3], type="l", lty=2, lwd=2)
points((beetle_dat$num.killed/beetle_dat$beetles)~beetle_dat$logdose, pch=19, col="grey2")

pi.est <- function(x,b1,b2)
{return(exp(b1+b2*x)/(1+exp(b1+b2*x)))}

logires <- (y/m) - pi.est(x,-60,35)
plot(logires~pi.est(x,-60,35), ylim=c(-0.5,0.5),ylab="residuals",xlab="fitted values",
     main="Logit Deviance Residuals")
abline(h=0)
abline(h=0.2, col="red")
abline(h=-0.2, col="red")

##quad loss
sum(logires^2)/8
