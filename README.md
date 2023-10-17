# Environmental-Effects-on-the-Spatio-Temporal-Variability-of-Sardine-Distribution

R Code corresponding to the manuscript: Environmental Effects on the Spatio-Temporal Variability of Sardine Distribution along the Portuguese Continental Coast.
This repository provides the R Code for four main approaches:

- Time Lagged Covariates using weighted averages

*Script_1_Create_kernel_covarites*: Having as input an original dataset (with the response, the date of measuring, the sampled locations, and the covariates observed during an interval of time lags - the provided example has 3 covariates measured in 29 days including the day when the response was measured),

| Longitude | Latitude | Date | Response | X1_before_0days | ... | X1_before_28days | X2_before_0days | ... | X2_before_28days | X3_before_0days | ... | X3_before_28days | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 406.78 | 4300.87 | 2017-05-15 | 1345.76 | 11.98 | ...  | 13.78 | 1.01 | ... | 1.21 | 0.00 | ... | 0.00 |
| 412.45 | 4210.75 | 2017-05-15 | 0 | 12.41 | ... | 12.54 | 1.32 | ... | 1.67 | 0.01 | ... | 0.03 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 409.63 | 4249.98 | 2022-05-10 | 23.90 | 12.76 | ... | 12.98 | ... | 2.15 | ... | 1.99 | 0.65 | ... | 0.54 |
| 410.99 | 4227.65 | 2022-05-10 | 2.01 | 14.41 | ... | 13.11 | 3.01 | ... | 2.34 | 1.17 | ... | 1.12 |

we transform the time-lagged covariates into kernel weighted covariates defining the maximum distance between the mode and the minimum of the time interval, that is, the maximum for the parameter 'l' in Equation (3) from the paper. The parameter 'l' can take different values for each covariate. The result is the completed dataset (the input dataset plus the kernel covariates).


- Preliminary evaluation of kernel covariates

*Script_2_Evaluate_kernel_covariates*: Using as input the completed dataset with the kernel covariates and the two responses (a binary covariate representing species presence/absence and a continuous covariate representing the estimated species biomass if the species is present), 

| Longitude | Latitude | Date | Response | X1_before_0days |  ... | X3_before_28days | X1_center_1_lag_1 | ... | X2_center_27_lag_1 | Presence | Pos_Response |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | 
| 406.78 | 4300.87 | 2017-05-15 | 1345.76 | 11.98 | ...  | 0.00 | 12.87 | ... | 0.02 | 1 | 1345.76 |
| 412.45 | 4210.75 | 2017-05-15 | 0 | 12.41 | ... | 0.03 | 14.12 | ... | 0.03 | 0 | NA |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 409.63 | 4249.98 | 2022-05-10 | 23.90 | 12.76 | ... | 0.54 | 13.98 | ... | 0.02 | 1 | 23.90 |
| 410.99 | 4227.65 | 2022-05-10 | 2.01 | 14.41 | ... | 1.12 | 16.01 | ... | 1.23 | 1 | 2.01 |

we create a ranking of kernel covariates for each response process. At the same time, we provide the optimized number of knots and spline type according to the AIC metric for each kernel covariate. 


- Barrier Model

*Script_3_Run_BModel*: Having as input the completed dataset and using the ranking information obtained through the previous step, we run the Barrier model using flat priors and the R-INLA package. The result is an INLA object that contains the posterior marginal distributions for both fixed and hyperparameters, the posterior distributions of the fitted values for both responses and the random effects.


- Occupancy areas definition

*Script_4_Occupancy_areas*: Having as input an INLA object obtained in the previous step, we compute the posterior median of the species biomass for the mesh points and the years of study, and the uncertainty linked with the posterior distribution of the interest (biomass) process. The estimation of the uncertainty requires the definition of the 'scale' parameter which is 4/5 by default. Then, the posterior median of the species biomass and the uncertainty are used to classify the study region according to species occupancy. To perform this, a threshold for uncertainty is required. The default is 0.5. Despite of this, we advocate for decisions driven by biological insights and expert judgment to ensure credible and meaningful outcomes.




References:

Bakka, H., Vanhatalo, J., Illian, J. B., Simpson, D., and Rue, H. (2019). Non-stationary gaussian models with physical barriers. Spatial Statistics, 29:268–288.

Krainski, E., Gómez-Rubio, V., Bakka, H., Lenzi, A., Castro-Camilo, D., Simpson, D., Lindgren, F., and Rue, H. (2019). Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA. Chapman & Hall/CRC Press, Boca Raton, FL.

Lindgren, F., Rue, H., and Lindstrom, J. (2011). An explicit link between gaussian fields and gaussian markov random fields: the stochastic partial differential equation approach. Statistical Methodology, 73(4):423–498.

