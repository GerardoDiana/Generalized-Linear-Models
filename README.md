# Generalized-Linear-Models

This repository contains graduate coursework related to (Bayesian) generalized linear models.

## HW2 ##
[Write-Up](https://github.com/GerardoDiana/Generalized-Linear-Models/blob/master/HW2/DG_GLM_HW2_WRITEUP.pdf) description: 
1. Derived the deviance from a Poisson glm
2. Expressed the Gamma distribution as a member of the exponential dispersion family. Derived deviance from a Gamma glm
3. Model Implementation of Poisson GLM and Quasi Poisson GLM on the fabric dataset.
4. Built a Poisson GLM to study the effect of the covariates (jet fuel concentration and organism strain) on the number of Ceriodaphnia organisms

My [R scripts](https://github.com/GerardoDiana/Generalized-Linear-Models/tree/master/HW2/R_code) for 3 and 4. Datasets can be found [here](https://github.com/GerardoDiana/Generalized-Linear-Models/tree/master/HW2) 

## HW3 ##
Note: complete write-up is under construction

Write-Up description: 
1. Using R, I fitted 4 Binomial GLMs to the beetles mortality data corresponding to 4 link functions: logit, probit, complementary log-log, and a more general (parameteric) link function which I refer to as the modified logit link function. I also obtained residuals, fitted values, and estimated dose-response curve for all models. I perform a model comparison with AIC and BIC.
2. Bayesian analysis of the beetle mortality data. Considering each of the Binomial GLM Models (minus the probit model), I designed and implemented a MCMC method to sample from the posterior distribution.  I plot the point and interval estimates of the dose-response curve for each model. To compare models I use the Gelfand&Gosh criterion (posterior predictive loss)
3. Expressed that the Inverse Gaussian distribution is a member of the exponential dispersion family. Obtained the scaled deviance for the comparison of the full inverse Gaussian model with the inverse Gaussian GLM.

My R scripts for 1 and 2:
1. [prob1_274HW3.R](https://github.com/GerardoDiana/Generalized-Linear-Models/blob/master/HW3/R_code/prob1_274HW3.R): with all 4 models
2. [prob2_274HW3.R](https://github.com/GerardoDiana/Generalized-Linear-Models/blob/master/HW3/R_code/prob2_274HW3.R) : cloglog and logit model, [prob2_c_274HW3.R](https://github.com/GerardoDiana/Generalized-Linear-Models/blob/master/HW3/R_code/prob2_c_274HW3.R): modified logit model

My Python scripts for 1 and 2:
1. [beetle_hw3_part1.ipynb](https://github.com/GerardoDiana/Generalized-Linear-Models/blob/master/HW3/Python_code/beetle_hw3_part1.ipynb): with all 4 models
2. [beetle_hw3_part2_a.ipynb](https://github.com/GerardoDiana/Generalized-Linear-Models/blob/master/HW3/Python_code/beetle_hw3_part2_a.ipynb): cloglog model, [beetle_hw3_part2_b.ipynb](https://github.com/GerardoDiana/Generalized-Linear-Models/blob/master/HW3/Python_code/beetle_hw3_part2_b.ipynb): logit model, [beetles_hw3_part2_c.ipynb](https://github.com/GerardoDiana/Generalized-Linear-Models/blob/master/HW3/Python_code/beetles_hw3_part2_c.ipynb): modified logit function using Pyro

## HW4 ##
[Write-Up](https://github.com/GerardoDiana/Generalized-Linear-Models/blob/master/HW4/DG_GLM_HW4_WRITEUP.pdf) description:
1. refer to the pdf write-up. (non-coding problem)
2. Using the alligator data, I developed a Bayesian multinomial regression model, using the baseline-category logits formulation with "fish" as the baseline category, to estimate the response probabilities as a function of length. And later extend the model to describe the effects of both length and gender on food choice.
3. Using the mice data, I fitted two Binomials GLMS from multinomial data and obtained MLE estimates. I then implement a Bayesian version of the model.

The data can be found [here](https://github.com/GerardoDiana/Generalized-Linear-Models/tree/master/HW4) and my R scripts can be found [here](https://github.com/GerardoDiana/Generalized-Linear-Models/tree/master/HW4/R_code)