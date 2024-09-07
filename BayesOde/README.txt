Introduction

Generating parameters and data

First, we generate 100 sets of parameters from the prior distribution. In Simulation I, we use the model-based approach to generate 100 sets of data. For Simulation II and subsequent comparative studies, including sensitivity analyses, the data is generated using the Gillespie algorithm. The code and corresponding datasets for data generation are located in the relevant folders.

Code for this part:
set_para:  Generating parameters from a priori
generateddata_model: Generating data from our model
generateddata_gillespie: Generating data from Gillespie algorithm

Simulation I
Verify algorithm assuming no missing data 

simulation1: 
main program----load datasets, perform parameter estimation, and compare the results with the true parameters

logL_nomissing: 
functions----log-likelihood without missing data 

accept_rate_com: 
function----calculate and output the standard deviation of the proposed distribution, adjusting based on the appropriate acceptance rate

Gibbs: 
function----perform a Gibbs cycle

check_convergence:
function----Determine whether the MCMC chain has converged and output point estimator, whether the confidence interval covers the true value, confidence interval

Simulation II
Verify algorithm assuming missing data 

simulation2
main program----load datasets, perform parameter estimation, and compare the results with the true parameters

logL_missing: 
functions----log-likelihood with missing data 

add_data_new:
function----expand the amount of data, inputs data with missing values, and outputs the expanded data

reset_yobsv：
function----arrange the inserted observations back into their corresponding positions

accept_rate_missing：
function----calculate and output the standard deviation of the proposed distribution, adjusting based on the appropriate acceptance rate

Gibbs_missing_not0beta_NotEqualLambda：
function----perform a Gibbs cycle

check_convergence_missing_not0beta_NotEqualLambda：
function----determine whether the MCMC chain has converged and output point estimator, whether the confidence interval covers the true value, confidence interval

IC：
function----calculation of the DIC value

maxlikeli:
function----calculate the likelihood of the phantom parameters of the mcmc chain and output the maximum likelihood

Sensitivity analysis
Change the initial total number of cells in simulation2 to test thethe sensitivity of the algorithm to the total number of initial cells

Realdata analysis

realdata_differentParas: 
main program----assuming that the parameters of the four groups are different, four separate models are used for each set of data, and the point estimates of the parameters and DIC values are calculated for each model. Gibbs_missing_not0beta_NotEqualLambda and check_convergence_missing_not0beta_NotEqualLambda are replaced accordingly.

realdata_sameParas:
main program----assuming that the parameters of the four groups are the same.
Four groups data are as an input, and DIC values are calculated for each model.

Comparison study

NLS: estimate parameters using nonlinear least square algorithm and RK4

smoothing method: estimate parameters using smoothing method

discrete model: estimate parameters using discrete-time Markov chain model