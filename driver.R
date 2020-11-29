# #### This is the overall driver function
# #### The first part fits the model to data from Eatern and Western Cambodia. You may skip this part if you 
# #### already have parameter sets. The second creates the model output and plots using the parameter sets
# #### NOTE: THIS FULL CODE WILL TAKE SEVERAL HOURS-1 DAY, RUNNING 8 PARALLEL CLUSTERS

# #######################################################################################
# ########################## FITTING MODEL PARAMETERS ###############################
# #######################################################################################

# ##### Output 1: east_pars_0_0.csv  (parameters sets for Eastern Cambodia)
# ##### Output 2: west_pars_0_0.csv (parameter sets for Western Cambodia)
# source('model_fitting.R')
# source('likelihoods.R')
# source('IMIS_mod_par.R')
# library(deSolve)
# library(ggplot2)
# library(logitnorm)
# library(IMIS)
# library(dplyr)
# library(parallel)

# ##### Fitting model - east
# source('malaria_IMIS_east.R')

# #### Fitting model - west
# source('malaria_IMIS_west.R')

# #######################################################################################
# ########################## PLOTTING FIGURES FOR PAPER ###############################
# #######################################################################################

# # Outputs here include: CSV files for each scenario and outcome (upper, median, mean, and lower for each year)
# # PDF images of each scenario and outcome

library(deSolve)
library(ggplot2)
library(dplyr)
library(reshape2)
library(logitnorm)
library(tidyverse)

source('plotting_function_2.R') # actually creates plots, based on output CSV files

# Parameters shared across regions
b = 1/10 # relative infectiousness asymptomatic infections
w =  1/365 # waning immunity
w_t = 1/20 # duration protected by treatment drug
r_a = 1/60 # cure rate asymptomatic
r_s_0 <- 1/10 # cure rate symptomatic if no relapse
s_a <- 0.5 # proportion receiving non-appropriate treatment (mostly artemisinin monotherapy) who succeed
					# reduced to account for artemisinin resistance in Cambodia
p_h <- 1/10 # progression from latency to infectiousness (human)
p_m <- 1/14 # progression from latency to infectiousness (mosquito)
mu_0 = 1/(30*365) # duration in population at risk at baseline
mu_m <- 1/8 # mosquito mortality
a_0 = 1e-6  # acquired resistance (on a given drug) # 1 in 1,000,000
 s_p <- 0.5  # proportion with multiple pfpm2 copies who are successfully treated with DP
s_m <- 0.75  # proportion with multiple pfmdr1 copies who are successfully treated with ASMQ
ppq_init <- 0.01
mosq_c = 0.5 # P(infection | bite) human -> mosquito
mosq_b = 0.5 # # P(infection | bite) mosquito -> human
mosq_a = 1/3  # human biting rate
	

# # # ######################### Baseline plots #################################
# source('model_scenarios.R') # version of the model that includes interventions starting 2020
# s_both <- 1 - (1- s_m)*(1- s_p)

# # ####################  EAST
# source('running_model_final_east.R')
# region = 2 # East
# mdr_2000 = 0.1 
# final_pars <- read.csv('east_pars_0_0.csv')
# variable_inputs_all <- subset(final_pars,select=-X)
# #variable_inputs_all <- variable_inputs_all[1:20,] # TEMPORARY

# fixed_inputs <- c(region=region, b=b, w=w, r_a=r_a, mu_0=mu_0, r_s_0=r_s_0,
				# s_a=s_a, p_h = p_h, p_m = p_m, mu_m = mu_m, w_t=w_t, a_0=a_0, 
				# s_p=s_p, s_m=s_m, mdr_2000=mdr_2000, ppq_init=ppq_init, s_both=s_both,
				# mosq_a = mosq_a, mosq_b = mosq_b, mosq_c = mosq_c) 
# fixed_inputs <- fixed_inputs[setdiff(names(fixed_inputs),colnames(final_pars))] # if input varies, don't keep fixed

# running_model_final_east(scenario="baseline",which_proph="NONE", fitness_sens_ana=F, fixed_inputs,variable_inputs_all) #no intervention scenario
# running_model_final_east(scenario="baseline",which_proph="PPQ", fitness_sens_ana=F, fixed_inputs, variable_inputs_all) #chemoprophylaxis with PPQ scenario
# running_model_final_east(scenario="baseline",which_proph="TRT", fitness_sens_ana=F, fixed_inputs, variable_inputs_all) # triple ACT treatment scenario

# #############  WEST

# source('running_model_final_west.R')
# region = 1 # West
# mdr_2000 = 0.15 
# final_pars <- read.csv('west_pars_0_0.csv')
# variable_inputs_all <- subset(final_pars,select=-X)
# #variable_inputs_all <- variable_inputs_all[1:20,] # TEMPORARY

# fixed_inputs <- c(region=region, b=b, w=w, r_a=r_a, mu_0=mu_0, r_s_0=r_s_0,
				# s_a=s_a, p_h = p_h, p_m = p_m, mu_m = mu_m, w_t=w_t, a_0=a_0, 
				# s_p=s_p, s_m=s_m, mdr_2000=mdr_2000, ppq_init=ppq_init, s_both=s_both,
				# mosq_a = mosq_a, mosq_b = mosq_b, mosq_c = mosq_c) 
# fixed_inputs <- fixed_inputs[setdiff(names(fixed_inputs),colnames(final_pars))] # if input varies, don't keep fixed

# running_model_final_west(scenario="baseline",which_proph="NONE", fitness_sens_ana=F, fixed_inputs, variable_inputs_all) #no intervention
# running_model_final_west(scenario="baseline",which_proph="PPQ", fitness_sens_ana=F, fixed_inputs, variable_inputs_all) #chemoprophylaxis with pPQ
# running_model_final_west(scenario="baseline",which_proph="TRT", fitness_sens_ana=F, fixed_inputs, variable_inputs_all) # triple ACT treatment

# # ######### Plotting
# plotting_function(scenario="baseline")

# ##################################### Short baseline (for plots showing model fit) #################################
# source('model_scenarios.R')
# s_both <- 1 - (1- s_m)*(1- s_p)

# ####################  EAST
# source('running_model_final_east.R')
# region = 2 # East
# mdr_2000 = 0.1 
# final_pars <- read.csv('east_pars_0_0.csv')
# variable_inputs_all <- subset(final_pars,select=-X)
# variable_inputs_all <- variable_inputs_all[1:100,] 

# fixed_inputs <- c(region=region, b=b, w=w, r_a=r_a, mu_0=mu_0, r_s_0=r_s_0,
				# s_a=s_a, p_h = p_h, p_m = p_m, mu_m = mu_m, w_t=w_t, a_0=a_0, 
				# s_p=s_p, s_m=s_m, mdr_2000=mdr_2000, ppq_init=ppq_init, s_both=s_both,
				# mosq_a = mosq_a, mosq_b = mosq_b, mosq_c = mosq_c) 
# fixed_inputs <- fixed_inputs[setdiff(names(fixed_inputs),colnames(final_pars))] # if input varies, don't keep fixed

# running_model_final_east(scenario="short_baseline",which_proph="NONE", fitness_sens_ana=F, fixed_inputs,variable_inputs_all)
# running_model_final_east(scenario="short_baseline",which_proph="PPQ",fitness_sens_ana=F, fixed_inputs, variable_inputs_all)
# running_model_final_east(scenario="short_baseline",which_proph="TRT",fitness_sens_ana=F, fixed_inputs, variable_inputs_all) 

# #############  WEST

# source('running_model_final_west.R')
# region = 1 # West
# mdr_2000 = 0.15 
# final_pars <- read.csv('west_pars_0_0.csv')
# variable_inputs_all <- subset(final_pars,select=-X)
# variable_inputs_all <- variable_inputs_all[1:100,]  

# fixed_inputs <- c(region=region, b=b, w=w, r_a=r_a, mu_0=mu_0, r_s_0=r_s_0,
				# s_a=s_a, p_h = p_h, p_m = p_m, mu_m = mu_m, w_t=w_t, a_0=a_0, 
				# s_p=s_p, s_m=s_m, mdr_2000=mdr_2000, ppq_init=ppq_init, s_both=s_both,
				# mosq_a = mosq_a, mosq_b = mosq_b, mosq_c = mosq_c) 
# fixed_inputs <- fixed_inputs[setdiff(names(fixed_inputs),colnames(final_pars))] # if input varies, don't keep fixed

# running_model_final_west(scenario="short_baseline",which_proph="NONE", fitness_sens_ana=F, fixed_inputs, variable_inputs_all)
# running_model_final_west(scenario="short_baseline",which_proph="PPQ",fitness_sens_ana=F, fixed_inputs, variable_inputs_all)
# running_model_final_west(scenario="short_baseline",which_proph="TRT",fitness_sens_ana=F, fixed_inputs, variable_inputs_all) 

# ######## Plotting
# #### No need here

# # ######################### Sensitivity analysis: changing f_mp only #################################
# ########## (fitness of strains resistance to both MQ and PPQ)
# source('model_scenarios.R')
# s_both <- 1 - (1- s_m)*(1- s_p)

# ####################  EAST
# source('running_model_final_east.R')
# region = 2 # East
# mdr_2000 = 0.1 
# final_pars <- read.csv('east_pars_0_0.csv')
# variable_inputs_all <- subset(final_pars,select=-X)
 # variable_inputs_all <- variable_inputs_all[1:200,] 

# fixed_inputs <- c(region=region, b=b, w=w, r_a=r_a, mu_0=mu_0, r_s_0=r_s_0,
				# s_a=s_a, p_h = p_h, p_m = p_m, mu_m = mu_m, w_t=w_t, a_0=a_0, 
				# s_p=s_p, s_m=s_m, mdr_2000=mdr_2000, ppq_init=ppq_init, s_both=s_both) 
# fixed_inputs <- fixed_inputs[setdiff(names(fixed_inputs),colnames(final_pars))] # if input varies, don't keep fixed

# running_model_final_east(scenario="f",which_proph="NONE", fitness_sens_ana=T, fixed_inputs, variable_inputs_all)
# running_model_final_east(scenario="f",which_proph="PPQ", fitness_sens_ana=T, fixed_inputs, variable_inputs_all)
# running_model_final_east(scenario="f",which_proph="TRT", fitness_sens_ana=T, fixed_inputs, variable_inputs_all) 

# #############  WEST

# source('running_model_final_west.R')
# region = 1 # West
# mdr_2000 = 0.15 # NEW
# final_pars <- read.csv('west_pars_0_0.csv')
# variable_inputs_all <- subset(final_pars,select=-X)
 # variable_inputs_all <- variable_inputs_all[1:200,] 

# fixed_inputs <- c(region=region, b=b, w=w, r_a=r_a, mu_0=mu_0, r_s_0=r_s_0,
				# s_a=s_a, p_h = p_h, p_m = p_m, mu_m = mu_m, w_t=w_t, a_0=a_0, 
				# s_p=s_p, s_m=s_m, mdr_2000=mdr_2000, ppq_init=ppq_init, s_both=s_both) 
# fixed_inputs <- fixed_inputs[setdiff(names(fixed_inputs),colnames(final_pars))] # if input varies, don't keep fixed

# running_model_final_west(scenario="f",which_proph="NONE", fitness_sens_ana=T, fixed_inputs, variable_inputs_all)
# running_model_final_west(scenario="f",which_proph="PPQ", fitness_sens_ana=T, fixed_inputs, variable_inputs_all)
# running_model_final_west(scenario="f",which_proph="TRT", fitness_sens_ana=T, fixed_inputs, variable_inputs_all) 

# # ######### Plotting
 # plotting_function(scenario="f")

# # ######################### Sensitivity analysis: changing f_mp and s_both #################################
########## (fitness of strains resistance to both MQ and PPQ) and (successful treatment if resistant to both)######

source('model_scenarios.R')
s_both <- max((s_m), (s_p)) 

###################  EAST
source('running_model_final_east.R')
region = 2 # East
mdr_2000 = 0.1 
final_pars <- read.csv('east_pars_0_0.csv')
variable_inputs_all <- subset(final_pars,select=-X)
variable_inputs_all <- variable_inputs_all[1:200,] 

fixed_inputs <- c(region=region, b=b, w=w, r_a=r_a, mu_0=mu_0, r_s_0=r_s_0,
				s_a=s_a, p_h = p_h, p_m = p_m, mu_m = mu_m, w_t=w_t, a_0=a_0, 
				s_p=s_p, s_m=s_m, mdr_2000=mdr_2000, ppq_init=ppq_init, s_both=s_both) 
fixed_inputs <- fixed_inputs[setdiff(names(fixed_inputs),colnames(final_pars))] # if input varies, don't keep fixed

running_model_final_east(scenario="f_s",which_proph="NONE", fitness_sens_ana=T, fixed_inputs, variable_inputs_all)
running_model_final_east(scenario="f_s",which_proph="PPQ", fitness_sens_ana=T, fixed_inputs, variable_inputs_all)
running_model_final_east(scenario="f_s",which_proph="TRT", fitness_sens_ana=T, fixed_inputs, variable_inputs_all) # joint treatment

############  WEST

source('running_model_final_west.R')
region = 1 # West
mdr_2000 = 0.15 # NEW
final_pars <- read.csv('west_pars_0_0.csv')
variable_inputs_all <- subset(final_pars,select=-X)
variable_inputs_all <- variable_inputs_all[1:200,] 


fixed_inputs <- c(region=region, b=b, w=w, r_a=r_a, mu_0=mu_0, r_s_0=r_s_0,
				s_a=s_a, p_h = p_h, p_m = p_m, mu_m = mu_m, w_t=w_t, a_0=a_0, 
				s_p=s_p, s_m=s_m, mdr_2000=mdr_2000, ppq_init=ppq_init, s_both=s_both) 
fixed_inputs <- fixed_inputs[setdiff(names(fixed_inputs),colnames(final_pars))] # if input varies, don't keep fixed

running_model_final_west(scenario="f_s",which_proph="NONE", fitness_sens_ana=T, fixed_inputs, variable_inputs_all)
running_model_final_west(scenario="f_s",which_proph="PPQ", fitness_sens_ana=T, fixed_inputs, variable_inputs_all)
running_model_final_west(scenario="f_s",which_proph="TRT", fitness_sens_ana=T, fixed_inputs, variable_inputs_all) # joint treatment

####### Plotting
plotting_function(scenario="f_s")