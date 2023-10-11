# Environmental-Effects-on-the-Spatio-Temporal-Variability-of-Sardine-Distribution

R Code corresponding to the manuscript: Environmental Effects on the Spatio-Temporal Variability of Sardine Distribution along the Portuguese Continental Coast.
This repository provides the R Code for four main approaches:

- Time Lagged Covariates using weighted averages

Script_1_Create_kernel_covarites: Having as input an original dataset (with the response, the date of measuring, the sampled locations, and the covariates observed during an interval of time lags - the provided example has 3 covariates measured in 29 days including the day when the response was measured), we transform the time lagged covariates in kernel weighted covariates defining the maximum distance between the mode and the minimum of the time interval, that is, the maximum for the parameter 'l' in Equation (3) from the paper. The parameter 'l' can take different values for each covariate. The result is the completed dataset (the input dataset plus the kernel covariates).


- Preliminary evaluation of kernel covariates

Script_2_Evaluate_kernel_covariates: Using as input a completed dataset with the kernel covariates and the two reponses (a binary covariate representing species presence/absence and a continuous covariate representing the estimated species biomass if the species is present), we create a ranking of kernel covariates for each response process. At the same time, we provide the optimized number of knots and spline type according to the AIC metric for each kernel covariate. 


- Barrier Model

Script_3_Run_BModel: Having as input the completed dataset and using the ranking information obtained throught the previous step, we run the Barrier model using flat priors and the R-INLA package. The result is an INLA object that contains the posterior marginal distributions for both fixed and hyperparameters, the posterior distributions of the fitted values for both responses and for the random effects.


- Occupancy areas definition

Script_4_Occupancy_areas: Having as input a INLA object obtained in the previous step, we compute the posterior median of the species biomass for the mesh points and the years of study and the uncertainty linked with the posterior distribution of the interest (biomass) process. The estimation of the uncertainty requires the definition of the 'scale' parameter which is 4/5 by default. Then, the posterior median of the species biomass and the regarding uncertainty are used to classify the study region according to species occupancy. To perform this, a threshold for uncertainty is required. The default is 0.5. Despite of this, we advocate for decisions driven by biological insights and expert judgment to ensure credible and meaningful outcomes.




References:

Bakka, H., Vanhatalo, J., Illian, J. B., Simpson, D., and Rue, H. (2019). Non-stationary gaussian models with physical barriers. Spatial Statistics, 29:268–288.

Krainski, E., Gómez-Rubio, V., Bakka, H., Lenzi, A., Castro-Camilo, D., Simpson, D., Lindgren, F., and Rue, H. (2019). Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA. Chapman & Hall/CRC Press, Boca Raton, FL.

Lindgren, F., Rue, H., and Lindstrom, J. (2011). An explicit link between gaussian fields and gaussian markov random fields: the stochastic partial differential equation approach. Statistical Methodology, 73(4):423–498.

