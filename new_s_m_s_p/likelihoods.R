# This file contains likelihoods for malaria reported cases, prevalence, and resistance for both Eastern and
# Western Cambodia.
# It is called by the IMIS functions

################# Monthly reported cases (max and min per year) #########################

## Western Cambodia

calc_monthly_west_llik <- function(results) {
	
	# Years. Typically the peak occurs in the fall and the low point in the spring
	dates <- seq(2004,2019,1)  
	dates_spring <- dates + 0.25
	dates_fall <- dates + 0.75

	# These data were extracted from Maude (2014) plots
	monthly_max_obs <- c(5160,2460,4800,2160,2200,3080,1580,1820,780,340) #2004-2013 peaks
	monthly_min_obs <- c(1660, 940,1820, 900,700, 1160,640, 780, 600, 220 ) #2004-2013 troughs
	
	# These data are estimated based on the World Malaria Report for more recent years
	monthly_max_est <- c(834,792,518,926,708,386) #2014-2019 peaks
	monthly_min_est <- c(361,343,224,400,307,167) #2014-2019 troughs
	monthly_max_obs <- c(monthly_max_obs, monthly_max_est)
	monthly_min_obs <- c(monthly_min_obs, monthly_min_est)
	
	# extracting model estimates
	model_spring <- results[results$year %in% dates_spring,]$rep + 0.0001 # add .0001 to avoid problems when =0
	model_fall <- results[results$year %in% dates_fall,]$rep + 0.0001
	
	# Log-likelihood. The /10 is introduced to allow for greater uncertainty in model outputs
	llik <- sum(dnbinom(monthly_min_obs, mu= model_spring, size = model_spring/10, log=T)) + 
		sum(dnbinom(monthly_max_obs, mu= model_fall, size = model_fall/10, log=T))		
	return(llik)	
}

## Eastern Cambodia
calc_monthly_east_llik <- function(results) {
	
	# Years. Typically the peak occurs in the fall and the low point in the spring
	dates <- seq(2004,2019,1) 
	dates_spring <- dates + 0.25
	dates_fall <- dates + 0.75

	# These data were extracted from Maude (2014) plots
	monthly_max_obs <- c(1451,1398,2193,1149,1369,1857,2065,1436,790,503) #2004-2013 peaks
	monthly_min_obs <- c(227,202,439,414,299,633,311,304,278,180) #2004-2013 troughs
	
	# These data are estimated based on the World Malaria Report for more recent years
	monthly_max_est <- c(906,861,563,1006,770,419) #2014-2019 peaks
	monthly_min_est <- c(220,209,137,245,187,102) #2014-2019 troughs
	monthly_max_obs <- c(monthly_max_obs, monthly_max_est)
	monthly_min_obs <- c(monthly_min_obs, monthly_min_est)
	
	# extracting model estimates
	model_spring <- results[results$year %in% dates_spring,]$rep + 0.0001 # add .0001 to avoid problems when =0
	model_fall <- results[results$year %in% dates_fall,]$rep + 0.0001

	 # Log-likelihood. The /10 is introduced to allow for greater uncertainty in model outputs
		llik <- sum(dnbinom(monthly_min_obs, mu= model_spring, size = model_spring/10, log=T)) + 
		sum(dnbinom(monthly_max_obs, mu= model_fall, size = model_fall/10, log=T))	
	return(llik)	
}

################# Malaria prevalence from surveys #########################

## Western Cambodia
calc_prev_west_llik <- function(results) {
	
	# Years surveys took place. All surveys occurred in the fall, around peak malaria season
	dates_prev <- c(2004.75,2007.75,2010.75,2013.75) 
	
	gen_slide_N <-c(2811,2943,4586,4114) # Number with slides tested
	gen_slide_N <-gen_slide_N/5 # The /5 is to allow for greater uncertainty in model outputs
	gen_slide_N_adj <- round(gen_slide_N)
	gen_slide_pf <- c(45,11,8,1) # Number with slides positive for P falciparum
	gen_slide_pf <- round(gen_slide_pf/5) 
	
	# Extracing model PCR prevalence and converting to slide prevalence
	# add .0001 to avoid problems when =0
	 model_pcr_prev <- results[results$year %in% dates_prev,]$prev_pcr 
	tmp <- (log((model_pcr_prev+ 0.0001)/(1-model_pcr_prev + 0.0001)) - 0.954)/0.868
	model_slide_prev <- exp(tmp)/(1+exp(tmp))
		
	# log likelihood
	llik <- sum(dbinom(gen_slide_pf, size= gen_slide_N_adj, prob= model_slide_prev, log=T)) 
	return(llik)	
}

## Eastern Cambodia
calc_prev_east_llik <- function(results) {
	
	# Years surveys took place. All surveys occurred in the fall, around peak malaria season
	dates_prev <- c(2004.75,2007.75,2010.75,2013.75) 
	
	gen_slide_N <- c(2885,3330,4505,4482) # Number with slides tested
		gen_slide_N <-gen_slide_N/5 # The /5 is to allow for greater uncertainty in model outputs
	gen_slide_N_adj <- round(gen_slide_N) 
	gen_slide_pf <-c(128,114,25,2) # Number positive for P falciparum
	gen_slide_pf <- round(gen_slide_pf/5)
	
	# Extracing model PCR prevalence and converting to slide prevalence
	# add .0001 to avoid problems when =0
	 model_pcr_prev <- results[results$year %in% dates_prev,]$prev_pcr 
	tmp <- tmp <- (log((model_pcr_prev+ 0.0001)/(1-model_pcr_prev + 0.0001)) - 0.954)/0.868
	model_slide_prev <- exp(tmp)/(1+exp(tmp))
		
	# log likelihood
	llik <- sum(dbinom(gen_slide_pf, size= gen_slide_N_adj, prob= model_slide_prev, log=T)) 
	return(llik)
	
}

########################### Resistance likelihoods ###############################

############# Multiple copies of pfmdr1 (genotypic MQ resistance) ##################################

### Western Cambodia
calc_pfmdr_west_llik <- function(results) {
	
	# Year collected. Assumed to be in the middle of the year(s) data were collected
	dates <- c(2008.5, 2009.5, 2010.5, 2011.5, 2012.5, 2013.5, 2011.5, 2012.5, 2013.5, 2008.5, 
		2009.5, 2010.5, 2011.5, 2009.5, 2010.5, 2011.5, 2012.5, 2013.5, 2011.5, 2012.5, 2011.5, 
		2012.5, 2013.5, 2010.5, 2011.5, 2012.5, 2013.5, 2008, 2005, 2003, 2014.5, 2014.5, 2016.5, 
		2017.5, 2014.5, 2015.5, 2016.5, 2018.5, 2015.5, 2016.5, 2018.5, 2015.5, 2016.5, 2015.5, 2016.5, 2017.5,
	    2019.5, 2019.5, 2019.5) 
	
	# Number of samples
	N_tested <- c(38, 91, 60, 110, 49, 64, 81, 72, 32, 21, 58, 37, 52, 163, 71, 57, 8, 42, 100, 
		15, 67, 34, 20, 53, 94, 22, 39, 40, 273, 224, 47, 8, 3, 56, 53, 48, 29, 43, 51, 40, 63, 53, 4, 96, 33, 200,
	     22, 13, 16) 
	N_tested <- round(N_tested/5)  # The /5 is to allow for greater uncertainty in model outputs
	
	# Number of samples with multiple copies of pfmdr1 (genotypic MQ resistance)
	N_multicopies <-c(19, 27, 20, 37, 5, 7, 18, 5, 4, 2, 11, 5, 11, 66, 25, 20, 2, 5, 54,
		 10, 8, 2, 3, 16, 37, 2, 3, 2, 144, 42, 5, 0, 1, 0, 0, 2, 0, 0, 1, 2, 5, 3, 0, 9, 16, 0,
		 4, 0, 2) 
		N_multicopies <- round(N_multicopies/5)
	 
		# Model output
	  model_pfmdr_prev <- results[results$year %in% dates,]$prev_pfmdr 
	  model_pfmdr_prev[model_pfmdr_prev < 0.0001] <- 0.0001 # to prevent issues when = 0 or 100%
	  model_pfmdr_prev[model_pfmdr_prev > 0.9999] <- 0.9999
	
	# Log likelihood
	llik <- sum(dbinom(N_multicopies, size= N_tested, prob= model_pfmdr_prev, log=T)) 
	return(llik)
	
}

## Eastern Cambodia
calc_pfmdr_east_llik <- function(results) {
	
	 # Year collected. Assumed to be in the middle of the year(s) data were collected
	dates <- c(2010.5, 2011.5, 2012.5, 2013.5, 2010.5, 2011.5, 2012.5, 2013.5, 
		 2017.5, 2015.5, 2016.5, 2016.5,  2017.5, 2018.5, 2017.5, 2017.5, 2018.5,  2016.5, 2017.5, 2015.5, 2016.5, 2018.5,
		2019.5, 2019.5) 

	# Number of samples
	N_tested <- c(54, 100, 52, 34, 49, 72, 23, 13,
		247, 52, 44, 55, 
		59, 53, 68, 60, 44, 31, 149, 48, 14, 24,
		8,52 ) 
	N_tested <- round(N_tested/5) # The /5 is to allow for greater uncertainty in model outputs
	
	# Number of samples with multiple copies of pfmdr1 (genotypic MQ resistance)
	N_multicopies <-c(0, 9, 0, 1, 0, 1, 0, 0, 1,  1, 15, 17,  0, 4, 1, 2, 0, 0, 2, 0, 0, 1,
	   0, 26 ) 
	N_multicopies <- round(N_multicopies/5)  
	 
	# Model output
	 model_pfmdr_prev <- results[results$year %in% dates,]$prev_pfmdr
	model_pfmdr_prev[model_pfmdr_prev < 0.0001] <- 0.0001 # to prevent issues when = 0 or 100%
	model_pfmdr_prev[model_pfmdr_prev > 0.9999] <- 0.9999
	
	# Log likelihood
	llik <- sum(dbinom(N_multicopies, size= N_tested, prob= model_pfmdr_prev, log=T)) 
	return(llik)
	
}

############## Multiple copies of pfpm2 (genotypic PPQ resistance) #########################

calc_pfpm2_west_llik <- function(results) {
	
	# Year collected. Assumed to be in the middle of the year(s) data were collected
	dates <- c(2009, 2009, 2011, 2011, 2011, 2013, 2013, 2013, 2013, 2013, 
		2013, 2015, 2015, 2011.5, 2012.5, 2013.5, 2010.5, 2011.5, 2012.5, 2013.5, 
		2015.5, 2017.5, 2015.5, 2016.5, 2015.5, 2016.5, 2018.5, 2015.5, 2016.5, 
		2018.5, 2016.5, 2015.5, 2015.5, 2016.5, 2017.5,
		2019.5, 2019.5, 2019.5) 
	
	# Number of samples
	N_tested <- c(68, 30, 21, 73, 34, 39, 23, 39, 22, 17, 16, 57, 57, 67, 34,
		 20, 53, 94, 22, 39, 14, 56, 53, 3, 48, 29, 43, 51, 41, 63, 4, 94, 129, 35, 200,
		 22, 13, 16) 
	N_tested <- round(N_tested/5)  # The /5 is to allow for greater uncertainty in model outputs
	
	# Number of samples with multiple copies of pfpm2 (genotypic PPQ resistance)
	N_multicopies <-c(19, 0, 13, 28, 3, 31, 13, 31, 1, 8, 3, 52, 52, 1, 1, 8, 
		13, 56, 14, 32, 10, 35, 43, 1, 41, 7, 7, 46, 22, 29, 3, 65, 90, 25, 165,
		8, 1, 2 ) 
		 N_multicopies <- round(N_multicopies/5)
	
		# Model output
	 model_pfpm_prev <- results[results$year %in% dates,]$prev_pfpm
	model_pfpm_prev[model_pfpm_prev < 0.0001] <- 0.0001 # to prevent issues when = 0 or 100%
	model_pfpm_prev[model_pfpm_prev > 0.9999] <- 0.9999
	
		# Log likelihood	
	llik <- sum(dbinom(N_multicopies, size= N_tested, prob= model_pfpm_prev, log=T)) 
	return(llik)
	
}

## Eastern Cambodia
calc_pfpm2_east_llik <- function(results) {
	 
	 # Year collected. Assumed to be in the middle of the year(s) data were collected
	dates <- c(2011, 2011, 2013, 2013, 2013, 2015, 2015, 2015, 2010.5, 
		2011.5, 2012.5, 2013.5, 2017.5, 2015.5, 2016.5, 2016.5, 2017.5, 
		2018.5, 2017.5, 2016.5, 2017.5, 2018.5, 2017.5, 2015.5, 2016.5, 2018.5,
		2019.5,2019.5) 

	# Number of samples
	N_tested <- c(30, 51, 31, 38, 22, 88, 55, 46, 49, 72, 23, 13, 247, 52, 45, 53, 
		59, 53, 68, 31, 60, 44, 149, 48, 15, 24,
			8, 52) 
	N_tested <- round(N_tested/5)   # The /5 is to allow for greater uncertainty in model outputs
	
	# Number of samples with multiple copies of pfpm2 (genotypic PPQ resistance)
	N_multicopies <-c(0, 1, 1, 3, 1, 40, 14, 29, 0, 1, 0, 2, 153, 26, 7, 30, 39, 8, 43, 4, 34, 7, 101, 36, 9, 6,
	1, 2) 
	N_multicopies <- round(N_multicopies/5)
	
		# Model output
	 model_pfpm_prev <- results[results$year %in% dates,]$prev_pfpm
	model_pfpm_prev[model_pfpm_prev < 0.0001] <- 0.0001 # to prevent issues when = 0 or 100%
	model_pfpm_prev[model_pfpm_prev > 0.9999] <- 0.9999
	
	# Log likelihood
	llik <- sum(dbinom(N_multicopies, size= N_tested, prob= model_pfpm_prev, log=T)) 
	return(llik)
	
}

################ Treatment efficacy studies  ###########################
########## ASMQ - late treatment failures - phenotypic MQ resistance ###############

## Western Cambodia
calc_mtes_west_llik <- function(results) {
	
	# Year collected. Assumed to be in the middle of the year(s) data were collected
	dates <- c(2001.5, 2002.5, 2002.5, 2002.5, 2003.5, 2003.5, 2003.5, 2004.5, 2004.5, 2004.5, 
		2011.5, 2014.5, 2015.5, 2017.5, 2010.5, 2016.5, 2017.5, 2015.5, 2014.5, 2016.5, 2003, 2007,
		 2019.5, 2019.5, 2019.5, 2018.5) # NEW

	# Number of individuals tested
	N_tested <- c(50, 70, 67, 29, 91, 90, 52, 81, 81, 85, 
	     28, 138, 59, 51, 45, 46, 55, 60, 59, 60, 162, 143,
		12, 22, 16, 11) 
	N_tested <- round(N_tested/5) # The /5 is to allow for greater uncertainty in model outputs
	
	# Number with late treatment failure
	N_LTF <-c(2, 10, 4, 1, 3, 2, 4, 8, 6, 0, 
	    0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 4, 27,
	    0, 0, 2, 0) 
	N_LTF <- round(N_LTF/5)   
	 
	# Model output
	model_mres_prev <- results[results$year %in% dates,]$prev_mres
	model_mres_prev[model_mres_prev < 0.0001] <- 0.0001 # to prevent issues when = 0 or 100%
	model_mres_prev[model_mres_prev > 0.9999] <- 0.9999
	
		# Log likelihood
	llik <- sum(dbinom(N_LTF, size= N_tested, prob= model_mres_prev, log=T)) 
	return(llik)
	
}

# Eastern Cambodia
calc_mtes_east_llik <- function(results) {
	
	# Year collected. Assumed to be in the middle of the year(s) data were collected
	dates <- c(2017.5, 2014.5, 2015.5, 2016.5, 2018.5, 2018.5, 2018.5, 2018.5, 2018.5,
		2019.5, 2019.5) 
	
	# Number of individuals tested
	N_tested <- c(59, 91, 54, 19, 29, 53, 6, 44, 24,
	  8, 50) 
	  N_tested <- round(N_tested/5)  # The /5 is to allow for greater uncertainty in model outputs
	
	# Number with late treatment failure
	N_LTF <-c(1, 1, 0, 0, 0, 1, 0, 0, 0,
	          0, 0) 
	 N_LTF <- round(N_LTF/5)# note the /5
	 
	# Model output
	model_mres_prev <- results[results$year %in% dates,]$prev_mres
	model_mres_prev[model_mres_prev < 0.0001] <- 0.0001 # to prevent issues when = 0 or 100%
	model_mres_prev[model_mres_prev > 0.9999] <- 0.9999
	
		# Log likelihood
	llik <- sum(dbinom(N_LTF, size= N_tested, prob= model_mres_prev, log=T)) 
	return(llik)
	
}

########## DP - late treatment failures - phenotypic PPQ resistance ###############

## Western Cambodia
calc_ptes_west_llik <- function(results) {
	
	# Year collected. Assumed to be in the middle of the year(s) data were collected
	dates <- c(2010.5, 2012.5, 2010.5, 2010.5, 2011.5, 2012.5, 2015.5, 2013.5, 2011.5,
		 2013.5, 2014.5, 2013.5, 2013.5, 2008.5, 2009.5, 2008.5, 2009.5, 2013.5)

	# Number of individuals tested
	N_tested <- c(28, 39, 11, 56, 42, 22, 53, 17, 59, 16, 40, 81, 63, 47, 37, 74, 60, 89)
	N_tested <- round(N_tested/5) # The /5 is to allow for greater uncertainty in model outputs
	
	# Number with late treatment failure
	N_LTF <-c(7, 12, 2, 6, 7, 0, 29, 1, 2, 1, 25, 37, 10, 5, 3, 1, 0, 42)
	N_LTF <- round(N_LTF/5) 
	
		# Model output
	 model_pres_prev <- results[results$year %in% dates,]$prev_pres
	model_pres_prev[model_pres_prev < 0.0001] <- 0.0001  # to prevent issues when = 0 or 100%
	model_pres_prev[model_pres_prev > 0.9999] <- 0.9999
	
		# Log likelihood
	llik <- sum(dbinom(N_LTF, size= N_tested, prob= model_pres_prev, log=T)) 
	return(llik)
	
}

## Eastern Cambodia
calc_ptes_east_llik <- function(results) {
	
	# Year collected. Assumed to be in the middle of the year(s) data were collected
	dates <- c(2012.5, 2011.5, 2013.5, 2014.5, 2014.5, 2015.5, 2016.5, 2010.5, 2013.5, 
		2016.5, 2013.5)
	
	# Number of individuals tested
	N_tested <- c(60, 55, 22, 40, 32, 60, 45, 59, 60, 58, 60)
		N_tested <- round(N_tested/5)   # The /5 is to allow for greater uncertainty in model outputs
	
	# Number with late treatment failure
	N_LTF <-c(0, 2, 0, 4, 11, 5, 12, 0, 1, 8, 1)
	N_LTF <- round(N_LTF/5) 

	# Model output
	model_pres_prev <- results[results$year %in% dates,]$prev_pres
	model_pres_prev[model_pres_prev < 0.0001] <- 0.0001 # to prevent issues when = 0 or 100%
	model_pres_prev[model_pres_prev > 0.9999] <- 0.9999
		
			# Log likelihood
	llik <- sum(dbinom(N_LTF, size= N_tested, prob= model_pres_prev, log=T)) 
	return(llik)
	
}

