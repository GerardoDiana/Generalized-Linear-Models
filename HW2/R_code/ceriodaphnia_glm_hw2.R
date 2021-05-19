############################ Generalized Linear Models ###########################
############################ Homework 3 Problem 4 Code ###########################

# By: Diana Gerardo

# Libraries ----------------------------------------------------------------------
library(MASS)
library(ggplot2)

# Problem 4 ----------------------------------------------------------------------

# The data. 
# The first column includes the number of organisms.
# The second the concentration of jet fuel.
# The third is the strain of the organism (binary)
cerio <- read.table("~/Generalized_Linear_Models/HW2/ceriodaphnia.csv", 
                           quote="\"", comment.char="")
colnames(cerio) <- c("num_organisms", "jet_fuel", "strain")
strain_color <- ifelse(cerio[,3] == 1, rgb(0,0,1,.5),rgb(1,0.2,0,.5))
head(cerio) 

# Renaming of variables ----------------------------------------------------------
response <- cerio$num_organisms
predictor <- cerio$jet_fuel

sqrt.response <- sqrt(cerio$num_organisms)
sqrt.predictor <- sqrt(cerio$jet_fuel)

log.response <- log(cerio$num_organisms)
log.predictor <- log(cerio$jet_fuel + 1)

response_dat <- rbind(response, sqrt.response, log.response)
predictor_dat <- rbind(predictor, sqrt.predictor, log.predictor)

response_dentitles <- rbind('Density of Ceriodaphnia',
                            'Density of Sqrt Ceriodaphnia',
                            'Density of Log Ceriodaphnia')
predictor_dentitles <- rbind('Density of Jet Fuel',
                             'Density of Sqrt Jet Fuel',
                             'Density of Log Jet Fuel+1')
response_names <- rbind('Ceriodaphnia','Sqrt Ceriodaphnia','Log Ceriodaphnia')
predictor_names <- rbind('Jet Fuel', 'Sqrt Jet Fuel', 'Log Jet Fuel+1')

# Function Helpers ---------------------------------------------------------------
density_plots <- function(y, x, y_titles, x_titles){
  n <- length(y[,1])
  
  for (i in 1:n){
    
    plot(density(y[i,]), main = y_titles[i,])
    plot(density(x[i,]), main = x_titles[i,])
    
  }
}

scatter_plots <- function(y,x, y_names, x_names){
  n <- length(y[,1])
  
  for (i in 1:n){
    
    for (j in 1:n){
      rounded = round(cor(y[i,], x[j,]),5)
      scatter.smooth(jitter(x[j,]),jitter(y[i,]), col = strain_color,pch=19,
                     main = paste(y_names[i],x_names[j],
                                  sep = ' vs '),
                     ylab = 'ceriodaphnias', xlab = "jet fuel")
      legend('topright', box.lty=0, inset=.02,
             legend = paste('r = ', rounded, sep = ''))
    }
  }
}

gg_density <- function(y,x, y_names, x_names){
  n <- length(y[,1])
  
  for (i in 1:n){
    for(j in 1:n){
      print(ggplot(cerio, aes (x = x[j,], y = y[i,], colour = factor(strain))) +
        stat_density2d () + expand_limits(y = c(-0.5, 12), x = c(-0.5,2.5)) +
        ggtitle(paste(y_names[i], x_names[j], sep = ' vs ')) + 
        xlab(x_names[j]) + ylab(y_names[i]))
    }
  }
}

print_summary <- function(model_list){
  for(i in 1:length(model_list)){
    print(model_list[[i]])
  }
}

model_criterion <- function(model_list, model_names){
  AIC <- c()
  BIC <- c()
  Res.Dev <- c()
  Chi.sq <- c()
  
  for(i in 1:length(model_list)){
    AIC[i] <- round(model_list[[i]]$aic,5)
    Res.Dev[i] <- round(model_list[[i]]$deviance,5)
    BIC[i] <- round(model_list[[i]]$aic - 2*3 + 3*log(length(response)),5)
    Chi.sq[i] <- round(sum(resid(model_list[[i]], type = "pearson")^2),5)
  }
  
  df <- cbind(model_names,AIC,Res.Dev,BIC,Chi.sq)
  return(df)
}

residual_plots <- function(model_list, model_names){
  for(i in 1:length(model_list)){
    pearson <- resid(model_list[[i]], type = "pearson")
    plot(pearson, ylim = c(-4.3,4.3),
         main = paste(model_names[i], 'Pearson Residuals', sep = ' '))
    abline(h=2);abline(h=-2)
    
    deviance <- resid(model_list[[i]], type = "deviance")
    plot(deviance,ylim = c(-4.3,4.3),
         main = paste(model_names[i], 'Deviance Residuals', sep = ' '))
    abline(h=2);abline(h=-2)
  }
}

# Exploring choice of transformation --------------------------------------------- 
density_plots(response_dat, predictor_dat,
              response_dentitles, predictor_dentitles)

scatter_plots(response_dat, predictor_dat,
              response_names, predictor_names)

gg_density(response_dat, predictor_dat,
           response_names, predictor_names)

ggplot(cerio, aes (x = predictor, y = response, colour = factor(strain))) +
  stat_density2d () + expand_limits(y = c(-20, 130), x = c(-0.5,2.5)) +
  ggtitle('Ceriodaphnia vs Jet Fuel') + xlab('jet fuel') + ylab('ceriodaphnia')

# Models -------------------------------------------------------------------------
M1 <- glm(num_organisms ~ jet_fuel + strain, data = cerio,
          family = poisson(link = 'log'))
M2 <- glm(num_organisms ~ jet_fuel + strain, data = cerio,
          family = poisson(link = 'sqrt')) 
M3 <- glm(num_organisms ~ log(jet_fuel+1) + strain, data = cerio,
          family = poisson(link = 'log'))
M4 <- glm(num_organisms ~ jet_fuel + strain + jet_fuel:strain, data = cerio,
          family = poisson(link = 'log'))
          
model <- list(M1,M2,M3,M4)
mod_names <- rbind('M1','M2','M3','M4')

print_summary(model)
model_criterion(model, mod_names)
residual_plots(model, mod_names)

# Estimated Plots ----------------------------------------------------------------

plot(jitter(log(response)) ~ jitter(predictor), col = strain_color, pch=19)
# strain = 0
abline(model[[1]]$coefficients[1],
       model[[1]]$coefficients[2],
       col ="darkorange", lwd=2) 
#strain = 1
abline(model[[1]]$coefficients[1] + model[[1]]$coefficients[3],
       model[[1]]$coefficients[2], col = "darkblue", lwd=2)

plot(jitter(log(response)) ~ jitter(predictor), col = strain_color, pch=19)
#strain=0
abline(model[[4]]$coefficients[1],
       model[[4]]$coefficients[2],
       col ="darkorange", lwd=2) 
#strain=1
abline(model[[4]]$coefficients[1] + model[[4]]$coefficients[3]
       + model[[4]]$coefficients[4],
       model[[4]]$coefficients[2],
       col = "blue", lwd=2)

