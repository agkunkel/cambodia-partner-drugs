#### This file contains code needed to fit the parameteres for Western Cambodia

options(error=browser)

# This function returns the overall log likelihood of the model based on a set or sets of input parameters
# inpars are those parameters that are allowed to vary; fixed parameters are set within the model
llikelihood <- function(inpars) {
	
	# it may take as input a single parameter set (vector), or multiple (matrix)
	if (is.vector(inpars)) {
		rows = 1
		inpars <- as.matrix(t(inpars))
	} else {
		rows = dim(inpars)[1]
	}
	
	# Fixed parameters (based on units of days)
	# See Appendix of paper for more details
	fixed_inputs <- NULL
	region = 1 # West 
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
	a_0 = 1e-6 # acquired resistance (on a given drug) # 1 in 1,000,000
	 s_p <- 0.625 # proportion with multiple pfpm2 copies who are successfully treated with DP
	 s_m <- 0.625 # proportion with multiple pfmdr1 copies who are successfully treated with ASMQ
	 mdr_2000 = 0.15 # proportion with multiple pfmdr1 copies in the year 2000
	 mosq_c = 0.5 # P(infection | bite) human -> mosquito
	mosq_b = 0.5 # P(infection | bite) mosquito -> human
	mosq_a = 1/3 # human biting rate

	# combining inputs
	fixed_inputs_vec <- c(region=region, b=b, w=w, r_a=r_a, mu_0=mu_0, r_s_0=r_s_0,
			s_a=s_a, p_h = p_h, p_m = p_m, mu_m = mu_m, w_t=w_t, a_0=a_0,  s_p=s_p, s_m=s_m,
			mdr_2000=mdr_2000, mosq_a = mosq_a, mosq_b = mosq_b, mosq_c = mosq_c)
				
	# add fixed inputs to each row of variable inputs
	rep.row<-function(x,no){
	   matrix(rep(x,each=no),nrow=no)
	}
	fixed_inputs <- rep.row(fixed_inputs_vec,rows)
	colnames(fixed_inputs) <- names(fixed_inputs_vec)
	pars <- as.data.frame(cbind(fixed_inputs,inpars))
	
	pars$f_b <- pars$f_m*pars$f_p # fitness cost if genotypic resistance to both MQ and PPQ

	# Time periods for fitting, corresponding to different treatment policies
	# Note how time periods differ for Eastern and Western versions of the code
	# for simplicity, treat as if 12 equal months of 30.5 days
	start_time = 0 # start date (days). letting day 0 correspond to jan 1 1995
	end_time_init = (2000-1995)*366
	end_time_init2 = (2008-1995)*366 # WEST
	end_time = (2020-1995)*366 # end date (days)
	times_init <- seq(start_time,end_time_init, by = 30.5)
	times_init2 <- seq(end_time_init, end_time_init2, by = 30.5)
	times <- seq(end_time_init2, end_time, by = 30.5)
	
	# if multiple input parameter sets, will loop through the code below
	llik_vec <- NULL
	if (is.vector(pars)) {
		loops = 1
	} else {
		loops = dim(pars)[1]
	}
	
	# For each parameter set
	for (i in 1:loops) {
	
	# Consider one parameter set at a time
	if (is.vector(pars)) {
		pars_vec <- pars
	} else {
		pars_vec <- pars[i,]
	}
	
	# Conditions to initalize the model (1995-2000)
	yinit <- c(S_m = 0.99, E_mn = 0, E_mm = 0, E_mp = 0, E_mb = 0,
			 I_mn = 0.01, I_mm = 0, I_mp = 0, I_mb = 0,
			 S_h = pars$N_init-100, E_hn = 0, E_hm = 0, E_hp = 0, E_hb = 0, 
			 I_hsn = 100, I_hsm = 0, I_hsp = 0, I_hsb = 0,
			 Rt = 0, R = 0, I_han = 0, I_ham = 0, I_hap = 0, I_hab = 0,
			 New_n = 0, New_m = 0, New_p = 0, New_b = 0)

	
	# Run the model for the first time period (1995-2000)
	possible_error1 <- tryCatch({	
		results_init <- as.data.frame(ode(y=yinit, times=times_init, 
			func=malaria_model, pars_vec, method='lsoda'))
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n");e})
		if(inherits(possible_error1, "error")) {
			llik_vec[i] <- -500000 # if errors produced, assign poor likelihood and skip ahead
		} else {	
	
	# Values from 2000 become new inputs for next time period
	yinit2 <- subset(results_init,select=-time)[dim(results_init)[1],]
	results_init <- results_init[-dim(results_init)[1],] # so don't repeat same date

	# Western Cambodia - introduce MQ resistance in mosquitoes in year 2000
	# Mosquitoes
	yinit2$E_mm <- pars_vec$mdr_2000*yinit2$E_mn
	yinit2$E_mn <- (1-pars_vec$mdr_2000)*yinit2$E_mn
	yinit2$I_mm <- pars_vec$mdr_2000*yinit2$I_mn
	yinit2$I_mn <- (1-pars_vec$mdr_2000)*yinit2$I_mn
	# Humans
	yinit2$E_hm <- pars_vec$mdr_2000*yinit2$E_hn
	yinit2$E_hn <- (1-pars_vec$mdr_2000)*yinit2$E_hn
	yinit2$I_hsm <- pars_vec$mdr_2000*yinit2$I_hsn
	yinit2$I_hsn <- (1-pars_vec$mdr_2000)*yinit2$I_hsn
	yinit2$I_ham <- pars_vec$mdr_2000*yinit2$I_han
	yinit2$I_han <- (1-pars_vec$mdr_2000)*yinit2$I_han
	ynames <- names(yinit2)
	yinit2 <- as.numeric(yinit2)
	names(yinit2) <- ynames
	
	# Run the model for second time period (2000-2008, West)
	possible_error2 <- tryCatch({	
		results_init2 <- as.data.frame(ode(y=yinit2, times=times_init2, 
			func=malaria_model, pars_vec, method='lsoda'))
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n");e})
		if(inherits(possible_error2, "error")) {
			llik_vec[i] <- -500000
		} else {	
	
	# 2008 outputs will be inputs for next run
	yinit3 <- subset(results_init2,select=-time)[dim(results_init2)[1],]
		results_init2 <- results_init2[-dim(results_init2)[1],] # so don't repeat same date


	# Introduce PPQ resistance 
	# Need to consider pre-existing MQ resistance as well
	ppq_init <- 0.01 # Assuming 1% to start
	
	# introducing PPQ resistance to mosquitoes
	yinit3$E_mp <- with(yinit3, ppq_init*E_mn*E_mn/(E_mn+E_mm))
	yinit3$E_mb <- with(yinit3, ppq_init*E_mn*E_mm/(E_mn+E_mm))
	yinit3$E_mn <- with(yinit3, (1-ppq_init)*E_mn)
	yinit3$I_mp <- with(yinit3, ppq_init*I_mn*I_mn/(I_mn+I_mm))
	yinit3$I_mb <- with(yinit3, ppq_init*I_mn*I_mm/(I_mn+I_mm))
	yinit3$I_mn <- with(yinit3, (1-ppq_init)*I_mn)
	
	# introducing PPQ resistance to humans
	yinit3$E_hp <- with(yinit3, ppq_init*E_hn*E_hn/(E_hn+E_hm))
	yinit3$E_hb <- with(yinit3, ppq_init*E_hn*E_hm/(E_hn+E_hm))
	yinit3$E_hn <- with(yinit3, (1-ppq_init)*E_hn)
	yinit3$I_hsp <- with(yinit3, ppq_init*I_hsn*I_hsn/(I_hsn+I_hsm))
	yinit3$I_hsb <- with(yinit3, ppq_init*I_hsn*I_hsm/(I_hsn+I_hsm))
	yinit3$I_hsn <- with(yinit3, (1-ppq_init)*I_hsn)
	yinit3$I_hap <- with(yinit3, ppq_init*I_han*I_han/(I_han+I_ham))
	yinit3$I_hab <- with(yinit3, ppq_init*I_han*I_ham/(I_han+I_ham))
	yinit3$I_han <- with(yinit3, (1-ppq_init)*I_han)
	
	ynames <- names(yinit3)
	yinit3 <- as.numeric(yinit3)
	names(yinit3) <- ynames
	
	
		### Running model for final time period	(2008-2020)
	possible_error <- tryCatch({	
		results <- as.data.frame(ode(y=yinit3, times=times, 
			func=malaria_model, pars_vec, method='lsoda'))
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n");e})
		if(inherits(possible_error, "error")) {
			llik_vec[i] <- -500000
		} else {	

	# Formatting the model output
	results <- rbind(results_init, results_init2, results)
	
	results$year <- results$time/366+1995
	
	# Number of humans
	results$Nh <- with(results,{S_h + E_hn + E_hm + E_hp + E_hb +
			 I_hsn + I_hsm + I_hsp + I_hsb +
			 Rt + R + I_han + I_ham + I_hap + I_hab})
		
	# Number of mosquitoes				 
	results$Nm <- with(results,{S_m + E_mn + E_mm + E_mp + E_mb +
				 I_mn + I_mm + I_mp + I_mb })
					 
	# Symtpomatic and asymptoamtic humans					 
	results$I_s <- with(results,{I_hsn + I_hsm + I_hsp + I_hsb})	
	results$I_a <- with(results,{I_han + I_ham + I_hap + I_hab})	
	
	# PCR prevalence (not just population at risk)
	results$prev_pcr <- with(results,{(I_s+I_a)/2700000}) 
	results$prev_pcr_fg <- with(results,{(I_s+I_a)/Nh})  # PCR prevalence in forest-goers (at-risk population)
	results$prev_mosq <- with(results,{(I_mn + I_mm + I_mp + I_mb)/Nm}) # prevalence in mosquitoes

	# undoing cumulative sum
	firstdiff <- function(x) {
	   shifted <- c(0,x[1:(length(x)-1)])
	   x-shifted
	   }
	
	results$New <- with(results,{New_n + New_m + New_p + New_b}) # incident cases	 per month
	results$rep <- pars_vec$p*firstdiff(results$New) # reported incident cases per month
	
	# assuming resistance is measured among new (incident) cases only
	# measure this by looking at force of infection from mosquitoes to humans
		results$I_m_inc <- with(results,{I_mn + pars_vec$f_m*I_mm + pars_vec$f_p*I_mp + pars_vec$f_b*I_mb })
	results$prev_pfmdr <- with(results,{(pars_vec$f_m*I_mm + pars_vec$f_b*I_mb)/(I_m_inc)})
	 results$prev_pfpm <-  with(results,{(pars_vec$f_p*I_mp + pars_vec$f_b*I_mb)/(I_m_inc)})
	results$prev_mres <- (1-pars_vec$s_m)*results$prev_pfmdr 
	results$prev_pres <- (1-pars_vec$s_p)*results$prev_pfpm

	
	################################ Likelihoods #########################
	## ## see llikelihoods.R
		
	# monthly reported cases
	monthly_west_llik <- calc_monthly_west_llik(results)	
	
	## malaria prevalence
	prev_west_llik <- calc_prev_west_llik(results)
	
	## Markers of resistance
	pfmdr_west_llik <- calc_pfmdr_west_llik(results)
	 pfpm2_west_llik <- calc_pfpm2_west_llik(results)
	
	## Phenotypic resistance
	mtes_west_llik <- calc_mtes_west_llik(results)
	ptes_west_llik <- calc_ptes_west_llik(results)
	
	# overall log likelihood
	  llik_west <-  monthly_west_llik   + prev_west_llik  +
			  pfmdr_west_llik + pfpm2_west_llik +
			 mtes_west_llik + ptes_west_llik #+ pmosq_west_llik
	llik <- llik_west
	if(length(llik)==0) {llik <- NA}
	llik_vec[i] <- llik #  if looping through multiple parameter sets
	} # try-catch
	} # 2nd try-catch
	} # 3rd try-catch
	} # for loop

	return((llik_vec)) # returns log likelihood
		
}

##### This function returns the prior density for a set of input parameters
prior <- function(pars) {
	
	  # Densities for each parameter are added
		density <-  log(dlogitnorm(pars[["a"]], mu=0,sigma=0.7)) #seasonal amplitude
		density <- density + log(dlogitnorm(pars[["p"]], mu=-1.55,sigma=0.3)) #P(symptomatic case notified)
		density <- density + log(dlogitnorm(pars[["f_m"]], mu=3.4,sigma=0.6)) # relative fitness, genotypic MQ resistance
		density <- density + log(dlogitnorm(pars[["f_p"]],  mu=3.4,sigma=0.6)) # relative fitness, genotypic PPQ resistance
		density <- density + dlnorm(pars[["mosq_m"]], meanlog=1,sdlog=1, log=T) # ratio mosquitoes:humans
		density <- density + dlnorm(pars[["N_init"]],  meanlog=12.9, sdlog=0.4, log=T) # Initial population at risk - West
		density <- density + log(dlogitnorm(pars[["beta_min_m"]],  mu=-0.2,sigma=1.02)) #multiplier after decline in man-biting rate
		
		return(exp(density)) # overall density is returned (not in log form)
			
}

# Samples from the prior distribution of each parameter. Parametrs same as in prior function above
sample.prior <- function(n) {
		
	variable_inputs <- NULL
	a <- rlogitnorm(n, mu=0, sigma=0.7)	
	p <- rlogitnorm(n, mu=-1.55,sigma=0.3)
	f_m <- rlogitnorm(n,  mu=3.4,sigma=0.6)
	f_p <-	rlogitnorm(n,  mu=3.4,sigma=0.6)
			mosq_m <- rlnorm(n, meanlog=1,sdlog=1)	
		N_init <- rlnorm(n,  meanlog=12.9, sdlog=0.4) # West: 
		beta_min_m <- rlogitnorm(n,  mu=-0.2,sigma=1.02)
		
	# Combining everything
	variable_inputs <- cbind(a, p, f_m, f_p, mosq_m, N_init, beta_min_m) 
	return(as.data.frame(variable_inputs))
				
}

set.seed(1) # to ensure reproducibility if re-run

# CALLING MODIFIED IMIS FUNCTION
out <- IMIS_mod_par(200,500,50,5) 
write.csv(x=out$resample, file='west_pars_0_0.csv') # write CSV file with final parameter set to file



