############################ Generalized Linear Models ###########################
########################### Homework 4 Problem 3b Code ###########################

# By: Diana Gerardo

#options(scipen = 9999999)
rm(list=ls()) # clear out global environment


# Libraries ----------------------------------------------------------------------

libraries <- c('MASS','mvtnorm','readxl')
sapply(libraries, require, character.only = TRUE)

mice <- read_excel("~/Generalized_Linear_Models/HW4/mice.xlsx")

# Adding columns -----------------------------------------------------------------
# prop represents a proportion
mice$prop1 <- mice$Dead/mice$num.subj 
mice$alive <- mice$num.subj - mice$Dead
mice$prop2 <- mice$Malform/mice$alive

# Frequentist --------------------------------------------------------------------
M1 <- glm(prop1 ~ Concentration, weights=num.subj, 
          family=binomial(link=logit), data=mice)

M2 <- glm(prop2 ~ Concentration, weights=alive,
          family=binomial(link=logit), data=mice)

new.concentration <- seq(0, 500, len=100)

p1 <- predict(M1,newdata=list(Concentration=new.concentration),type="response")
p2 <- predict(M2,newdata=list(Concentration=new.concentration),type="response")

pi1 <- p1
pi2 <- p2 * (1-pi1)
pi3 <- 1 - pi1 - pi2

plot(mice$Concentration, seq(0,1,length=length(mice$Concentration)),type='n',
     fg='grey', xlab='Concentration (mg/kg per day)', ylab='response porportion', 
     las=1, main=expression(hat(pi)[j](x)),ylim=c(0,1))

points(mice$Concentration, mice$prop1, pch=19, col='red')
lines(new.concentration, pi1, col='red')

points(mice$Concentration, mice$prop2*(1-mice$prop1), pch=19, col='darkgreen')
lines(new.concentration, pi2, col='darkgreen')

points(mice$Concentration, 1-mice$prop1-mice$prop2*(1-mice$prop1),
       pch=19, col='blue')
lines(new.concentration, pi3, col='blue')

legend("left",legend=c('1: Dead','2: Malformation','3: Normal'),
       lty=1,bty='n',col=c("red","darkgreen","blue"))

# Bayesian -----------------------------------------------------------------------

### Relabing variables as described in hw ###
y1 <- mice$Dead
y2 <- mice$Malform
x <- mice$Concentration
m <- mice$num.subj


##### Prior, Pi, and Posterior Functions #####
flat.prior <- function (para){ return(log(1)) }

posterior.1 <- function(beta){
  Xb <- beta[1]+beta[2]*x
  numer <- exp(Xb)
  denom <- 1+exp(Xb)
  log.1 <- log(numer/denom)
  log.2 <- log(1 - numer/denom)
  inside.sum <- y1*log.1 + (m-y1)*log.2
  posterior <- exp(flat.prior(beta) + sum(inside.sum))
  
  return(posterior)
}

pi_1 <- function(alpha1,beta1){
  Xb <- alpha1 + beta1*x
  numer <- exp(Xb)
  denom <- 1 + exp(Xb)
  val <- numer/denom
  
  return(val)
} 

pi_2 <- function(beta,alpha1,beta1){
  Xb <- beta[1]+beta[2]*x
  numer <- exp(Xb)
  denom <- 1 + exp(Xb)
  val <- (numer/denom) * 1/(1+exp(alpha1 + beta1*x))
  
  return(val)
}

pi_3 <-function(beta,alpha1,beta1){
  return(1 - pi_1(alpha1,beta1) - pi_2(beta,alpha1,beta1) )
}

posterior.2 <- function(beta,alpha1,beta1){
  numer.1 <- pi_2(beta,alpha1,beta1)
  numer.2 <- pi_3(beta,alpha1,beta1)
  denom <- pi_2(beta,alpha1,beta1)+pi_3(beta,alpha1,beta1)
  log.1 <- log(numer.1/denom)
  log.2 <- log(numer.2/denom)
  inside.sum <- y2*log.1 + (m-y1-y2)*log.2
  posterior <- exp(flat.prior(beta) + sum(inside.sum))
  
  return(posterior)
}

########### M-H Algorithm Functions ############
M.H.1 <- function(iters, glm, tune){
  # initialize parameters
  iter.beta <- matrix(NA, nrow=2, ncol=iters)
  iter.beta[,1] <- c(0,0)
  
  new.beta <- matrix(NA, nrow=2, ncol=iters)
  new.beta[,1] <- c(1,1) 
  
  J <- vcov(glm)
  d <- tune
  
  # M-H process
  for(i in 2:iters){
    new.beta[,i] <- rmvnorm(1, iter.beta[,i-1], d*J)
    
    q <- min(1, posterior.1(new.beta[,i])/posterior.1(iter.beta[,i-1]))
    
    if(runif(1) < q){iter.beta[,i] <- new.beta[,i]}
    
    else {iter.beta[,i] <- iter.beta[,i-1]}
  }
  
  # Discard initial simulations
  T <- iters/2
  iter.beta <- iter.beta[,-(1:T)]
  
  return(iter.beta)
}

M.H.2 <- function(iters, glm, tune, alpha1, beta1){
  # initialize parameters
  iter.beta <- matrix(NA, nrow=2, ncol=iters)
  iter.beta[,1] <- c(0,0)
  
  new.beta <- matrix(NA, nrow=2, ncol=iters)
  new.beta[,1] <- c(1,1) 
  
  J <- vcov(glm)
  d <- tune
  
  for(i in 2:iters){
    new.beta[,i] <- rmvnorm(1, iter.beta[,i-1], d*J)
    
    top <- posterior.2(new.beta[,i],alpha1,beta1)
    bottom <- posterior.2(iter.beta[,i-1],alpha1,beta1)
    
    q <- min(1, top/bottom)
    
    if(runif(1) < q){iter.beta[,i] <- new.beta[,i]}
    
    else {iter.beta[,i] <- iter.beta[,i-1]}
  }
  
  # Discard initial simulations
  T <- iters/2
  iter.beta <- iter.beta[,-(1:T)]
  
  return(iter.beta)
}


######### Sampling From The Posterior ##########
mle.alpha1 <- M1$coefficients[1] 
mle.beta1 <- M1$coefficients[2]

set.seed(80000)
sample.1 <- M.H.1(iters = 20*10^4, glm = M1, tune = 39)

sample.2 <- M.H.2(iters = 20*10^4, glm = M2, tune = 4, mle.alpha1, mle.beta1)

######### Check For Healthy Convergence ##########
samps <- list(sample.1,sample.2)
labels <- list(c(expression(alpha[1]), expression(beta[1])),
               c(expression(alpha[2]), expression(beta[2])) )

trace_plots <- function(samples,para_list){
  for(i in 1:length(samples)){
    plot.ts(samples[[i]][1,], ylab=para_list[[i]][1], xlab="index")
    plot.ts(samples[[i]][2,], ylab=para_list[[i]][2], xlab="index")
  }
}

trace_plots(samps, labels)

plot(density(samps[[2]][1,]), main=expression(alpha[2]~" Distribution"),
     xlab=expression(alpha[2]))
plot(density(samps[[2]][2,]), main=expression(beta[2]~" Distribution"),
     xlab=expression(beta[2]))

############ 95% Confidence Interval #############
confidence_interval <- function(samples, para_vector_names){
  interval <- apply(samples,1,quantile,c(0.025,0.5,0.975))
  colnames(interval) <- para_vector_names
  interval<- rbind(interval, apply(samples,1,sd))
  rownames(interval) <- c('2.5%','50%','97.5%','sd')
  
  return(interval)
}

table1 <- confidence_interval(samps[[1]], c('alpha1',"beta1"))
table1

table2 <- confidence_interval(samps[[2]], c('alpha2','beta2'))
table2

############ Estimated Response Curve #############

PI.1 <- function(x,samps1){
  Xb <- samps1[1,] + samps1[2,]*x
  numer <- exp(Xb)
  denom <- 1 + exp(Xb)
  val <- numer/denom
  
  return(val)
}

PI.2 <- function(x,samps1,samps2){
  Xb <- samps2[1,]+samps2[2,]*x
  numer <- exp(Xb)
  denom <- 1 + exp(Xb)
  val <- (numer/denom) * 1/(1+exp(samps1[1,] + samps1[2,]*x))
  
  return(val)
}

PI.3 <- function(x,samps1,samps2){
  return(1-PI.1(x,samps1)-PI.2(x,samps1,samps2))
}

quant <- function(x, samps1, samps2, choice){
  q <- matrix(NA,length(x),3)
  
  if(choice == 1){
    for(i in 1:length(x)){
      new.pi<- PI.1(x[i], samps1)
      q[i,] <- quantile(new.pi,c(0.025,0.5,0.975))
    }
    return(q)
  }
  
  if(choice == 2){
    for(i in 1:length(x)){
      new.pi<- PI.2(x[i], samps1, samps2)
      q[i,] <- quantile(new.pi,c(0.025,0.5,0.975))
    }
    return(q)
  }
  
  if(choice == 3){
    for(i in 1:length(x)){
      new.pi<- PI.3(x[i], samps1, samps2)
      q[i,] <- quantile(new.pi,c(0.025,0.5,0.975))
    }
    return(q)
  }
}



# dose-response graph
x.grid <- seq(0, 500, length=100)
q1 <- quant(x.grid,samps[[1]],samps[[2]],1)
q2 <- quant(x.grid,samps[[1]],samps[[2]],2)
q3 <- quant(x.grid,samps[[1]],samps[[2]],3)

plot(x.grid,q1[,1], type = "l", xlab="Concentration (mg/kg per day)", 
     ylab="response porportion",ylim=c(0,1), col="white",
     main=expression(hat(pi)[j](x)))

polygon(c(rev(x.grid), x.grid), 
        c(rev(q1[ ,3]), q1[ ,1]), 
        col =rgb(1,0,0,alpha=0.3), border = NA)

polygon(c(rev(x.grid), x.grid), 
        c(rev(q2[ ,3]), q2[ ,1]), 
        col =rgb(0,1,0,alpha=0.3), border = NA)

polygon(c(rev(x.grid), x.grid), 
        c(rev(q3[ ,3]), q3[ ,1]), 
        col =rgb(0,0,1,alpha=0.3), border = NA)


lines(x.grid,q1[,2], col="red")
points((mice$Dead/mice$num.subj)~ mice$Concentration, pch=19, col="red")

lines(x.grid,q2[,2], col="darkgreen")
points(mice$Concentration, mice$prop2*(1-mice$prop1), pch=19, col='darkgreen')

lines(x.grid, q3[,2], col="blue")
points(mice$Concentration, 1-mice$prop1-mice$prop2*(1-mice$prop1), pch=19, col='blue')

legend("left",legend=c('1: Dead','2: Malformation','3: Normal'),lty=1,bty='n',col=c("red","darkgreen","blue"))


