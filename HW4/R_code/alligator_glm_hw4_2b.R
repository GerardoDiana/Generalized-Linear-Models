############################ Generalized Linear Models ###########################
########################### Homework 4 Problem 2b Code ###########################

# By: Diana Gerardo

rm(list=ls()) # clear out global environment
options(scipen = 999)

# Libraries ----------------------------------------------------------------------

libraries <- c('MASS','stargazer','mvtnorm', 'nnet','foreign', 'readxl')
sapply(libraries, require, character.only = TRUE)

alligators <- read_excel("~/Generalized_Linear_Models/HW4/alligators.xlsx")

# Defining X, y, z ---------------------------------------------------------------
y_col_choice <- function(x, letter){
  column <- rep(0, length(x))
  
  for(i in 1:length(x)){
    if(alligators$choice[i] == letter) column[i] <- 1
    else {column[i] <- 0}
  }
  return(column)
}

z_col_gender <- function(x){
  column <- rep(0, length(x))
  
  for(i in 1:length(x)){
    if(alligators$gender[i] == "male" ) column[i] <- 1
    else {column[i] <- 0}
  }
  return(column)
}

X <- alligators$length

y.I <- y_col_choice(X,'I') 
y.O <- y_col_choice(X,'O')
y.F <- y_col_choice(X,'F')

Z <- z_col_gender(X) # binary version of gender

# MLE ----------------------------------------------------------------------------
multi <- multinom(alligators$choice ~ alligators$length + alligators$gender,
                  Hess = TRUE)
summary(multi)

# Prior, Pi, Posterior Functions -------------------------------------------------
flat.prior <- function(para){ return(log(1)) }

pi.I <- function(para,x,z){
  Xb.1 <- para[1]+para[2]*x + para[3]*z
  Xb.2 <- para[4]+para[5]*x + para[6]*z
  numer <- exp(Xb.1)
  denom <- 1+exp(Xb.1)+exp(Xb.2)
  
  return(numer/denom)
}

pi.O <- function(para,x,z){
  Xb.1 <- para[1]+para[2]*x + para[3]*z
  Xb.2 <- para[4]+para[5]*x + para[6]*z
  numer <- exp(Xb.2)
  denom <- 1+exp(Xb.1)+exp(Xb.2)
  
  return(numer/denom)
}

pi.F <- function(para,x,z){
  Xb.1 <- para[1]+para[2]*x + para[3]*z
  Xb.2 <- para[4]+para[5]*x + para[6]*z
  numer <- 1
  denom <- 1+exp(Xb.1)+exp(Xb.2)
  
  return(numer/denom)
}

posterior <- function(para,y.I,y.O,y.F,x,z){
  I.term <- y.I * log(pi.I(para,x,z))
  O.term <- y.O * log(pi.O(para,x,z))
  F.term <- y.F * log(pi.F(para,x,z))
  inside.sum <- I.term + O.term + F.term
  inside.exp <- flat.prior(para) +  sum(inside.sum)
  
  return( exp(inside.exp) ) 
}

# M-H Algorithm ------------------------------------------------------------------
set.seed(2000)

iters <- 200000
param <- matrix(NA, nrow=6, ncol=iters)
param[,1] <- rep(0,6)

new.param <- matrix(NA, nrow=6, ncol=iters)
new.param[,1] <- rep(0,6)

J <- vcov(multi)
d <- .25


for(i in 2:iters)
{
  new.param[,i] <- rmvnorm(1, param[,i-1], d*J)
  
  q <- min(1,posterior(new.param[,i],y.I,y.O,y.F,X,Z)/posterior(param[,i-1],
                                                              y.I,y.O,y.F,X,Z))
  
  if(runif(1) < q)
  {param[,i] <- new.param[,i]}
  
  else
  {param[,i] <- param[,i-1]}
  
}

# Discard Initial Simulations ----------------------------------------------------
burn <- function(samps, discard){
  T <- discard
  return(samps[,-(1:T)])
}

samples <- burn(param,iters/2)
dim(samples)

# Check Convergence --------------------------------------------------------------
traceplot <- function(data, ylab_list){
  for(i in 1:length(data[,1])){
    plot.ts(data[i,], ylab = ylab_list[i], xlab = 'index')
  }
}

lab.list <- c(expression(alpha[1]),expression(beta[1]), expression(lambda[1]),
              expression(alpha[2]),expression(beta[2]), expression(lambda[2]))
traceplot(samples,lab.list)

# 95% Confidence Interval of sampels
confidence_interval <- apply(samples,1,quantile, c(0.025,0.5,0.975))
colnames(confidence_interval) <- c("alpha1","beta1","lambda1",
                                   "alpha2","beta2","lambda2")
confidence_interval

# Density Plots of samples
density_plots <- function(samps, xlab_list){
  for(i in 1:length(samps[,1])){
    plot(density(samps[i,]),xlab=xlab_list[i],
         main=xlab_list[i]) 
  }
}

density_plots(samples,lab.list)


# Estimated and Interval Plots ---------------------------------------------------
Pi.I <- function(para,x,z){
  Xb.1 <- para[1,]+para[2,]*x + para[3,]*z
  Xb.2 <- para[4,]+para[5,]*x + para[6,]*z
  numer <- exp(Xb.1)
  denom <- 1+exp(Xb.1)+exp(Xb.2)
  
  return(numer/denom)
}

Pi.O <- function(para,x,z){
  Xb.1 <- para[1,]+para[2,]*x + para[3,]*z
  Xb.2 <- para[4,]+para[5,]*x + para[6,]*z
  numer <- exp(Xb.2)
  denom <- 1+exp(Xb.1)+exp(Xb.2)
  
  return(numer/denom)
}

Pi.F <- function(para,x,z){
  Xb.1 <- para[1,]+para[2,]*x + para[3,]*z
  Xb.2 <- para[4,]+para[5,]*x + para[6,]*z
  numer <- 1
  denom <- 1+exp(Xb.1)+exp(Xb.2)
  
  return(numer/denom)
}

quant <- function(x, z, samps, choice){
  q <- matrix(NA,length(x),3)
  
  if(choice == "I"){
    for(i in 1:length(x)){
      new.pi<- Pi.I(samps,x[i],z[i])
      q[i,] <- quantile(new.pi,c(0.025,0.5,0.975))
    }
    return(q)
  }
  
  if(choice == "O"){
    for(i in 1:length(x)){
      new.pi<- Pi.O(samps,x[i],z[i])
      q[i,] <- quantile(new.pi,c(0.025,0.5,0.975))
    }
    return(q)
  }
  
  if(choice == "F"){
    for(i in 1:length(x)){
      new.pi<- Pi.F(samps,x[i],z[i])
      q[i,] <- quantile(new.pi,c(0.025,0.5,0.975))
    }
    return(q)
  }
}

x.grid <- seq(0, 4, length=100)
z.grid <- rep(x=1,times=100) # x=1 male, x=0 female

I.quant <- quant(x.grid, z.grid, samples, "I")
O.quant <- quant(x.grid, z.grid, samples, "O")
F.quant <- quant(x.grid, z.grid, samples, "F")

plot(x.grid,I.quant[,1], type = "l", xlab="Length of alligator (meters)", 
     ylab="predicted probability",ylim=c(0,1), col="white",
     main=expression(hat(pi)[j](x)~" (Male)"))

polygon(c(rev(x.grid), x.grid), 
        c(rev(I.quant[ ,3]), I.quant[ ,1]), 
        col =rgb(1,0,0,alpha=0.3), border = NA)

polygon(c(rev(x.grid), x.grid), 
        c(rev(O.quant[ ,3]), O.quant[ ,1]), 
        col =rgb(0,1,0,alpha=0.3), border = NA)

polygon(c(rev(x.grid), x.grid), 
        c(rev(F.quant[ ,3]), F.quant[ ,1]), 
        col =rgb(0,0,1,alpha=0.3), border = NA)

lines(x.grid, I.quant[,2], col="red")

lines(x.grid, O.quant[,2], col="darkgreen")

lines(x.grid, F.quant[,2], col="blue")

legend(-0.1,.55,legend=c('Invertebrates','Other','Fish'),lty=1,bty='n',col=c("red","darkgreen","blue"))
