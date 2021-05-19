#Diana Gerardo
#AMS 274 Homework 3
#Problem 1

library(readxl)
beetle_dat <- read_excel("~/Generalized_Linear_Models/HW3/beetle_dat.xlsx")
# View(beetle_dat)

############################### Part (a) ##################################
y <- (beetle_dat$num.killed/beetle_dat$beetles)

# LOGIT Bin GLM
logit <- glm((num.killed/beetles) ~ logdose, weights = beetles ,family = binomial(link="logit"),data = beetle_dat)

# PROBIT Bin GLM
probit <- glm((num.killed/beetles) ~ logdose, weights = beetles, family = binomial(link="probit"), 
              data = beetle_dat)
# Complementary log-log (CLL) Bin GLM
CLL <- glm((num.killed/beetles) ~ logdose, weights = beetles, family = binomial(link="cloglog"), 
           data = beetle_dat)


#Deviance Residuals
logit.resid <- residuals(logit,type = "deviance")
probit.resid  <- residuals(probit, type = "deviance")
CLL.resid <- residuals(CLL, type = "deviance")


#Residual Deviance
l.resdev <- sum(logit.resid^2)
p.resdev <- sum(probit.resid^2)
C.resdev <- sum(CLL.resid^2)


#Mean Squared Residuals
l.msr <- sqrt(mean(logit.resid^2))
p.msr <- sqrt(mean(probit.resid^2))
C.msr <- sqrt(mean(CLL.resid^2))
  

#Residual plots
plot(logit.resid ~ logit$fitted.values, ylab="residuals", xlab="fitted values", 
     main="Model 1 Deviance Residuals")
abline(h=0)
abline(h=1, col="red")
abline(h=-1, col="red")


plot(probit.resid ~ probit$fitted.values, ylab="residuals", xlab="fitted values", 
     main = "Model 2 Deviance Residuals")
abline(h=0)
abline(h=1, col="red")
abline(h=-1, col="red")

plot(CLL.resid ~ CLL$fitted.values,ylab="residuals", ylim=c(-1,1.5),
     xlab="fitted values",main = "Model 3 Deviance Residuals")
abline(h=0)
abline(h=1, col="red")
abline(h=-1, col="red")


#observed against fitted values
plot(logit$fitted.values ~ y , 
     ylab="fitted values",xlab="# killed / # beetles",data=beetle_dat, col="blue")
points(probit$fitted.values ~ y , data=beetle_dat, 
       col="green3", pch = 2)
points(CLL$fitted.values ~ y , data=beetle_dat, 
       col="purple", pch = 3)
abline(0,1)
legend(0.1,1, legend = c("logit","probit","CLL"),col=c("blue","green3","purple"),
       pch=1:3, box.lty = 0)

#Pred_pi (fitted values)
logit.pred_pi <- predict(logit,type="response")
probit.pred_pi <- predict(probit,type="response")
CLL.pred_pi <- predict(CLL,type="response")


## grid predictions
xx <- seq(1.65, 1.9, length=100)
l.pred <- predict(logit, newdata = list(logdose=xx), type="response")
p.pred <- predict(probit, newdata = list(logdose=xx), type="response")
C.pred <- predict(CLL, newdata = list(logdose=xx), type="response")

plot(xx,l.pred,type='n',fg='grey',bty='n',
     ylab="# killed / # beetles",xlab="log dose")
lines(xx, l.pred,col="blue", lwd = 2)
lines(xx,p.pred, col="green3", lwd= 2)
lines(xx,C.pred, col="purple", lwd = 2)
points(beetle_dat$logdose,y,pch=19,col='grey30')
legend(1.65,1, legend = c("logit","probit","CLL","Data"),col=c("blue","green3","purple","grey30"),
       lty=c(1,1,1,0), pch=c(NA,NA,NA,19), box.lty = 0)



########################################### Part (c) ###############################################

inv.mod.logit <- function(x, b, a) {
  eta <- c(x %*% b)
  exp(a*eta) / (1+exp(eta))^ a
}
b.hat <- c(-113.625, 62.5)
a.hat <- .279
mod_logit.pred <- inv.mod.logit(cbind(1,xx), b.hat, a.hat)


plot(xx,l.pred,type='n',fg='grey',bty='n',
     ylab="# killed / # beetles",xlab="log dose")
lines(xx, l.pred, col = "blue", lwd = 2)
lines(xx,p.pred, col="green3", lwd= 2)
lines(xx,C.pred, col="purple", lwd = 2)
lines(xx,mod_logit.pred, col="red", lwd=2)
points(beetle_dat$logdose,y,pch=19,col='grey30')
legend(1.65,1, legend = c("logit","probit","CLL","modlogit","Data"),col=c("blue","green3","purple","red","grey30"),
       lty=c(1,1,1,1,0), pch=c(NA,NA,NA,NA,19), box.lty = 0)


p.hat <- inv.mod.logit(cbind(1,beetle_dat$logdose), b.hat, a.hat)
m <- beetle_dat$beetles
y.hat <- p.hat * m
y.obs <- beetle_dat$num.killed
n <- length(y.obs)
dev.new <- sign(y.obs-y.hat) * sqrt(2*abs(y.obs*log(y.obs/y.hat) 
                                          + (m-y.obs)*log(1E-10+(m-y.obs)/(m-y.hat))))
plot(dev.new~p.hat,ylab="residuals",xlab="fitted values", ylim=c(-1,1),
     main="Model 4 Deviance Residuals")
abline(h=0)
abline(h=-1, col="red")
abline(h=1,col="red")



########################################### Part (d) ###############################################

logchoose <- function(m,y) lgamma(m+1) - lgamma(y+1) - lgamma(m-y+1)


K <- 3 # number of parameters
aic <- c(sapply(models,function(m) m$aic),-2 * sum(logchoose(m,y.obs) + y.obs*log(p.hat) + (m-y.obs)*log(1-p.hat)) + 2*K)
names(aic)[4] <- "logit2"
bic <- aic + c(2,2,2,3) * (log(n)-2)

aic
bic

beetle_dat$prob <- beetle_dat$num.killed/beetle_dat$beetles

link <- c("logit","probit","cloglog")
models <- lapply(as.list(link), function(lnk){
  glm(prob~logdose,
      weights=beetles,
      family=binomial(link=lnk),
      data=beetle_dat)
})
names(models) <- sapply(models,function(m) m$family$link)
