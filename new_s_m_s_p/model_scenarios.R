# Model that includes post-2020 interventions
# Used to plot results after model has been fit 

malaria_model <- function(times,yinit,pars){
	with(as.list(c(yinit,pars)), {
		
		# size of human population at risk
		Nh = S_h + E_hn + E_hm + E_hp +  E_hmp + 
			 I_hsn + I_hsm + I_hsp +  I_hsmp +
			 Rt + R + I_han + I_ham + I_hap + I_hamp + 
			 S_pr + # prophylaxis variables not present in fitting version of model
			 E_hp_pr + E_hmp_pr  +
			  I_hsp_pr + I_hsmp_pr  +
			 I_hap_pr + I_hamp_pr  +
			  R_pr + Rt_pr
			 
			  # size of mosquito population
		Nq = S_q + E_qn + E_qm + E_qp + E_qmp + 
			 I_qn + I_qm + I_qp + I_qmp 
	
		year = times/366+1995 
		
		# # proportion receiving appropriate treatment increases from 2004-2010
		if (year < 2004) {
			app <- 0.3
		} else if (year < 2010) {
			app <- 0.3 + (0.8 - 0.3)/(2010-2004)*(year-2004)
		} else {
			app <- 0.8
		}
		
	
		r_sn = r_s_0*app + r_s_0*(1-app)*s_a # successful treatment if no drug resistance
		rf_sn = r_s_0*(1-app)*(1-s_a) # failed treatment if no drug resistance (due to inappropriate treatment)
	   acq_rate = r_s_0*a_0*app/(1-a_0)  # acquire resistance if no drug resistance and appropriate treatment

		
		# Changes in first line drug regimen over time
		if (year < 2000) {# Before the introduction of ACTs
			
			r_sm <- r_sn  # successful treatment
			r_sp <- r_sn
			r_smp <- r_sn
			rf_sm = rf_sn #  failed treatment
			rf_sp = rf_sn 
			rf_smp = rf_sn 
			a_m = 0   #acquired resistance - none
			a_p = 0
			prot_t_m = 1 #post-treatment, protected against MQ and PPQ resistant strains
			prot_t_p = 1
				
		} else if (year < 2008) { # using MQ in whole country
			
			r_sm <- r_s_0*app*s_m + r_s_0*(1-app)*s_a  #successful treatment if genotyipc MQ resistance 
			r_sp <- r_sn
			r_smp <- r_s_0*app*s_m + r_s_0*(1-app)*s_a
			
			rf_sm = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)    # failed treatment - due to resistance or inappropriate treatment
			rf_sp = rf_sn 
			rf_smp = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)  
			
			a_m = acq_rate    #acquired resistance - only on appropriate treatment
			a_p = 0
			if (region ==2) { # Eastern Cambodia, introduce MQ resistance to model in 2010 (see Appendix)
				a_m = 0 } 
			
			prot_t_m = 0 # post-treatment protection against PPQ-resistant but not MQ-resistant strains
			prot_t_p = 1
			
		} else if (year < 2010) { # Policy varies by region
			
			if (region ==1) { # west # using DP in western Cambodia
				r_sm <- r_sn  # successful treatment
				r_sp <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
				r_smp <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
				rf_sm = rf_sn # failed treatment
			   rf_sp =  r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a) 
			   rf_smp = r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a)
				a_m = 0    # acquired resistance
				a_p = acq_rate
				prot_t_m = 1 # post-treatment protection against MQ rsistant strains
				prot_t_p = 0
			} else if (region ==2) { #east # still using ASMQ in eastern Cambodia
				r_sm <- r_s_0*app*s_m + r_s_0*(1-app)*s_a  # successful treatment
				r_sp <- r_sn
				r_smp <- r_s_0*app*s_m + r_s_0*(1-app)*s_a
				rf_sm = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)  # failed treatment
		    	rf_sp = rf_sn 
			    rf_smp = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)  
				 a_m = 0 #  # acquired MQ resistance introduced in model in Eastern Cambodia in 2010
				a_p = 0
				prot_t_m = 0
				prot_t_p = 1
			}	
			
				}  else if (year < 2017) { # Using DP in whole country 
			r_sm <- r_sn # successful treatment
			r_sp <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
			r_smp <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
			rf_sm = rf_sn # failed treatment
			rf_sp = r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a) 
			rf_smp = r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a) 
			a_m = 0  # acquired resistance
			a_p = acq_rate
			prot_t_m = 1 # post-treatment protection
			prot_t_p = 0
			
		} else if (year < 2018) {   # Switch back to ASMQ occurs first in Western Cambodia
			
			if (region ==1) { # west # using MQ in western Cambodia
				r_sm <- r_s_0*app*s_m + r_s_0*(1-app)*s_a # successful treatment
				r_sp <- r_sn
				r_smp <- r_s_0*app*s_m + r_s_0*(1-app)*s_a
				rf_sm = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)  # failed treatment
		    	rf_sp = rf_sn
			    rf_smp = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)  
				a_m = acq_rate # acquired resistance
				a_p = 0
				prot_t_m = 0
				prot_t_p = 1
			} else if (region ==2) { #east # still using DP in eastern Cambodia
				r_sm <- r_sn # successful treatment
				r_sp <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
				r_smp <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
				rf_sm = rf_sn #  failed treatment
			   rf_sp = r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a) 
			   rf_smp = r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a) 
				a_m = 0   # acquired resistance
				a_p = acq_rate
				prot_t_m = 1
				prot_t_p = 0
			}	
			
		} else { # using ASMQ in whole country
			
			r_sm <- r_s_0*app*s_m + r_s_0*(1-app)*s_a #successful treatment
			r_sp <- r_sn
			r_smp <- r_s_0*app*s_m + r_s_0*(1-app)*s_a
			rf_sm = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a) # failed treatment
			rf_sp = rf_sn 
			rf_smp = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a) 
			a_m = acq_rate    # acquired resistance
			a_p = 0
			prot_t_m = 0
			prot_t_p = 1
			
		}
		
		############ TRIPLE ACT SCENARIO
	joint_trt = 0 # this is an indicator for whether triple ACT treatment is applied
	# Triple ACT treatment
		if (which_proph < 2.1 & which_proph > 1.9) { # scenario 2=ppq+mq treatment
			if (year > 2020) { #introduce in 2000
				r_sm <- r_sn # successful treatment
				r_sp <- r_sn
				 r_smp <- r_s_0*app*s_both + r_s_0*(1-app)*s_a 
				rf_sm = rf_sn # failed treatment
			     rf_sp = rf_sn # failed treatment
			   rf_smp = r_s_0*app*(1-s_both) + r_s_0*(1-app)*(1-s_a) 
				a_m = acq_rate    # acquired resistance
				a_p = acq_rate 
				prot_t_m = 1 # treatment is effective against both MQ and PPQ resistant bugs separately
				prot_t_p = 1
				joint_trt = 1
			}
		} 
		
		
		############### PROPHYLAXIS SCEnARIO
		
		pr_start <- 2020 # start of prophylaxis enrollment period
		pr_end <- 2020 + 1/12 # end of prophylaxis enrollment period
		pr_prop <- 0.5 # proportion of the population that gets prophylaxis was 0.5
		pr_r <- 1/(12*30.5) # removal from prophylaxis (1 year) 
		
		pr_s <- 0
		pr_p <- 0 
		pr_sp <- 0
		
		# Piperaquine prophylaxis
		if (which_proph < 1.1 & which_proph > 0.9) { # 1=ppq
				pr_p <- 1 # a prophylaxis indicator
				if (year > pr_start & year < pr_end) {
					pr_s <-  -log(1-pr_prop)/((pr_end-pr_start)*366) #starting prophylaxis (not PPQ resistant)
					pr_sp <- -log(1-pr_prop*s_p)/((pr_end-pr_start)*366) #starting prophylaxis, and effective (PPQ resistant)
					}
		} 

		
			# population exit rates (declining population at risk)
		# input, mu_0, is in terms of days similar to other model inputs
		# 2010-2013: decline from about 7.5% to 5% in both E and W = 0.67
		# 2004 - 2007: decline from about 20% to 15% overall (not given separately) = 0.75
		# 2007-2010: question appears to have changed from work in forest to spend night in forest
		start_decl = 2004 # start earlier
		lambda = log(0.7)/((2013 - 2010)*366) # a compromise
		if (year < start_decl) {
			mu = mu_0
		} else if (year < 2018) {
			mu = mu_0 - lambda # note lambda is negative
		} else { # currently, holding steady from 2018 onwards
			mu = mu_0
		}
		
			# transmission rates/parameters
		# also allowed to decline over time, from 2008-2012, matching timing of ITN scale-up
		if (year < 2008) {
			mosq_a_t = mosq_a
		} else if (year < 2012) {
			mosq_a_t = mosq_a - (mosq_a - mosq_a*beta_min_m)*(year-2008)/(2012-2008)
		} else {
			mosq_a_t = mosq_a*beta_min_m
		}
		beta_q = mosq_a_t*mosq_c
		beta_h = mosq_m*mosq_a_t*mosq_b
		
		
		# trying to force trough near april, peak near september through 
		amp <- (1 - a*sin(2*pi*(times)/366)) 	 # was +45
			
		#  arbitrarily, applying fitness costs and seasonality only one direction (mosquito -> human)	
				
		beta_qn <- beta_q # human infecting mosquito doesn't change over time, let's say
		beta_qm <- beta_qn
		beta_qp <- beta_qn
		beta_qmp <- beta_qn
		
		beta_hn <- beta_h*amp # accounting for seasonality
		beta_hm <- beta_hn*f_m # and fitness costs
		beta_hp <- beta_hn*f_p
		beta_hmp <-  beta_hn*f_mp
		
		I_hn <- I_hsn + b*I_han # infectiousness from symptomatic + asymptomatic
		I_hm <- I_hsm + b*I_ham
		I_hp <- (I_hsp + I_hsp_pr) + b*(I_hap + I_hap_pr)
		I_hmp <- (I_hsmp + I_hsmp_pr) + b*(I_hamp + I_hamp_pr)
		
		I_ha <- I_han + I_ham + I_hap +  I_hamp +  I_hap_pr + I_hamp_pr  # sum of all asymptomatic
		
		####################### MOSQUITO EQUATIONS
		
		# Holding population sizes fixed and letting beta vary for seasonality
		# (Considering only % infected, not raw numbers)
		
		# Susceptible mosquitoes 
		dS_q <- mu_m*(Nq - S_q) +
			- beta_qn*I_hn/Nh*S_q - beta_qm*I_hm/Nh*S_q +
			- beta_qp*I_hp/Nh*S_q - beta_qmp*I_hmp/Nh*S_q 
		
			# Exposed mosquitoes
		dE_qn <- beta_qn*I_hn/Nh*S_q - (mu_m + p_m)*E_qn  # no resistance
		
		dE_qm <- beta_qm*I_hm/Nh*S_q - (mu_m + p_m)*E_qm # genotypic MQ resistance
		
		dE_qp <- beta_qp*I_hp/Nh*S_q - (mu_m + p_m)*E_qp # genotypic PPQ resistance
		
		dE_qmp <- beta_qmp*I_hmp/Nh*S_q - (mu_m + p_m)*E_qmp # genotypic MQ and PPQ resistance 
			
		# Infectious mosquitoes
		dI_qn <- p_m*E_qn - mu_m*I_qn
		
		dI_qm <- p_m*E_qm - mu_m*I_qm
		
		dI_qp <- p_m*E_qp - mu_m*I_qp
			
		dI_qmp <- p_m*E_qmp - mu_m*I_qmp
		

		###################### HUMAN EQUATIONS
		
		# Susceptible humans
		dS_h <- mu_0*Nh +
			- beta_hn*I_qn*S_h - beta_hm*I_qm*S_h +
			- beta_hp*I_qp*S_h - beta_hmp*I_qmp*S_h +
			+ w*R - pr_s*S_h + pr_r*S_pr - mu*S_h 
		
		dS_pr <- pr_s*(S_h + E_hn + E_hm) + # On prophylaxis
		 		 pr_sp*(E_hp + E_hmp) + 
				- (mu + pr_r)*S_pr +
				- pr_p*(1-s_p)*beta_hp*I_qp*S_pr +
				- pr_p*(1-s_p)*beta_hmp*I_qmp*S_pr  + w*R_pr

	# Exposed humans
	
		dE_hn <- beta_hn*I_qn*S_h - (mu + p_h)*E_hn - pr_s*E_hn
		 
		dE_hm <- beta_hm*I_qm*S_h - (mu + p_h)*E_hm - pr_s*E_hm
		
		dE_hp <- beta_hp*I_qp*S_h - (mu + p_h)*E_hp - pr_s*E_hp +
				   pr_r*E_hp_pr 
				 
		dE_hp_pr <- 	(pr_s - pr_sp)*(E_hp)	 - (mu + pr_r + p_h)*E_hp_pr  + # prophylaxis
				pr_p*(1-s_p)*beta_hp*I_qp*S_pr
		
		dE_hmp <- beta_hmp*I_qmp*S_h - (mu + p_h)*E_hmp - pr_s*E_hmp +
				   pr_r*E_hmp_pr
		 
		dE_hmp_pr <- 		(pr_s - pr_sp)*(E_hmp)  - (mu + pr_r + p_h)*E_hmp_pr +
				pr_p*(1-s_p)*beta_hmp*I_qmp*S_pr
	
	# Infectious, symptomatic humans
		
		dI_hsn <-  p_h*E_hn - mu*I_hsn +
					- r_sn*I_hsn + - 0.1*rf_sn*I_hsn +
					- a_m*I_hsn - a_p*I_hsn 
		
		dI_hsm <-  p_h*E_hm - mu*I_hsm +
					- r_sm*I_hsm +  - 0.1*rf_sm*I_hsm +
					- a_p*I_hsm + a_m*I_hsn 
		
		dI_hsp <-  p_h*E_hp - mu*I_hsp +
					- r_sp*I_hsp + - 0.1*rf_sp*I_hsp +
					- a_m*I_hsp + a_p*I_hsn	 +
					pr_r*I_hsp_pr							
		
		dI_hsp_pr <- p_h*E_hp_pr - (mu + pr_r)*I_hsp_pr + 
					- r_sp*I_hsp_pr  + - 0.1*rf_sp* I_hsp_pr # as if treated with DP
					
		dI_hsmp <-  p_h*E_hmp - mu*I_hsmp +
					- r_smp*I_hsmp  + - 0.1*rf_smp*I_hsmp +
					a_m*I_hsp + a_p*I_hsm  +
					pr_r*I_hsmp_pr
		
		dI_hsmp_pr	<- p_h*E_hmp_pr - (mu + pr_r)*I_hsmp_pr + 
					- r_smp*I_hsmp_pr 	+ - 0.1*rf_smp* I_hsmp_pr 	# as if treatment with DP

		# protected by treatment drug 
		
		dRt <- app*(r_sn*I_hsn + r_sm*I_hsm + r_sp*I_hsp + r_smp*I_hsmp) +
			- (1-s_m)*(1-prot_t_m)*beta_hm*I_qm*Rt - (1-s_p)*(1-prot_t_p)*beta_hp*I_qp*Rt +
			- ((1-s_m)*(1-prot_t_m) + (1-s_p)*(1-prot_t_p) + (1-s_both)*joint_trt)*beta_hmp*I_qmp*Rt +
			 - w_t*Rt - mu*Rt  + pr_r*Rt_pr
			
			# added protection by prophylaxis, if applicable
		dRt_pr <- app*(r_sp*I_hsp_pr + r_smp*I_hsmp_pr) +
				- (1-s_p)*(1-prot_t_p)*beta_hp*I_qp*Rt_pr +
				- ((1-s_both)*(1-prot_t_m) + (1-s_p)*(1-prot_t_p) + (1-s_both)*joint_trt)*beta_hmp*I_qmp*Rt_pr  +
				- w_t*Rt_pr - mu*Rt_pr  - pr_r*Rt_pr
			
	
		# partially immune
		
		dR <- (1-app)*(r_sn*I_hsn + r_sm*I_hsm + r_sp*I_hsp + r_smp*I_hsmp ) +
			w_t*Rt + r_a*(I_ha - I_hap_pr - I_hamp_pr)  +
			- (beta_hn*I_qn + beta_hm*I_qm)*R + 
			- (beta_hp*I_qp + beta_hmp*I_qmp)*R +
			- w*R - mu*R - pr_s*R + pr_r*R_pr
			
				# partially immune and protected by prophylaxis
			
		dR_pr <- (1-app)*(r_sp*I_hsp_pr + r_smp*I_hsmp_pr) +
			w_t*Rt_pr +  r_a*(I_hap_pr + I_hamp_pr)  +
				pr_s*(R + I_han + I_ham) + # getting on proph
		 		 pr_sp*(I_hap + I_hamp) + # getting on proph - asymptomatic I's. 
		 		 - (mu + w + pr_r)*R_pr + 
				- pr_p*(1-s_p)*beta_hp*I_qp*R_pr +
				- pr_p*(1-s_p)*beta_hmp*I_qmp*R_pr  
		
		# infectious, asymptomatic 
		
			rc = 0.05 # recombination coefficient
		
		dI_han <- beta_hn*I_qn*R - r_a*I_han - mu*I_han - pr_s*I_han + 0.1*rf_sn*I_hsn +
					+ 0.5*beta_hn*I_qn*(I_ham + I_hap) + (0.5 - rc)*beta_hn*I_qn*I_hamp + 
					- 0.5*(beta_hm*I_qm + beta_hp*I_qp)*I_han - (0.5 + rc)*beta_hmp*I_qmp*I_han + 
					+ rc*beta_hm*I_qm*I_hap + rc*beta_hp*I_qp*I_ham 
		
		dI_ham <- beta_hm*I_qm*R + (1-s_m)*(1-prot_t_m)*beta_hm*I_qm*Rt + 0.1*rf_sm*I_hsm +
			- r_a*I_ham - mu*I_ham - pr_s*I_ham  +
			+ 0.5*beta_hm*I_qm*(I_han + I_hamp) + (0.5 - rc)*beta_hm*I_qm*I_hap +  
			- 0.5*(beta_hn*I_qn + beta_hmp*I_qmp)*I_ham - (0.5 + rc)*beta_hp*I_qp*I_ham  + 
			+ rc*beta_hn*I_qn*I_hamp + rc*beta_hmp*I_qmp*I_han 
		
		dI_hap <- beta_hp*I_qp*R + (1-s_p)*(1-prot_t_p)*beta_hp*I_qp*Rt + 0.1*rf_sp*I_hsp +
			- r_a*I_hap - mu*I_hap - (pr_s)*I_hap +  pr_r*I_hap_pr +
			+ 0.5*beta_hp*I_qp*(I_han + I_hamp) + (0.5 - rc)*beta_hp*I_qp*I_ham + 
			- 0.5*(beta_hn*I_qn + beta_hmp*I_qmp)*I_hap - (0.5 + rc)*beta_hm*I_qm*I_hap + 
			+ rc*beta_hn*I_qn*I_hamp + rc*beta_hmp*I_qmp*I_han 
			
			  dI_hap_pr <- (pr_s - pr_sp)*I_hap  - (mu + pr_r + r_a)*I_hap_pr + 0.1*rf_sp* I_hsp_pr +
			  			pr_p*(1-s_p)*beta_hp*I_qp*R_pr +
			 			 + (1-s_p)*(1-prot_t_p)*beta_hp*I_qp*(Rt_pr) +
			 			 + 0.5*beta_hp*I_qp*I_hamp_pr +  
						- 0.5*beta_hmp*I_qmp*I_hap_pr   
			
		dI_hamp <- beta_hmp*I_qmp*R  +  0.1*rf_smp*I_hsmp +
		         - r_a*I_hamp - mu*I_hamp - (pr_s)*I_hamp  +
				 pr_r*I_hamp_pr +((1-s_m)*(1-prot_t_m) + (1-s_p)*(1-prot_t_p) + (1-s_both)*joint_trt)*beta_hmp*I_qmp*Rt +
			+ 0.5*beta_hmp*I_qmp*(I_ham + I_hap) + (0.5 - rc)*beta_hmp*I_qmp*I_han + 
			- 0.5*(beta_hm*I_qm + beta_hp*I_qp)*I_hamp - (0.5 + rc)*beta_hn*I_qn*I_hamp + 
			+ rc*beta_hm*I_qm*I_hap + rc*beta_hp*I_qp*I_ham 

		    
		   dI_hamp_pr <- (pr_s - pr_sp)*I_hamp  - (mu + pr_r + r_a)*I_hamp_pr  +  0.1*rf_smp* I_hsmp_pr +
				 pr_p*(1-s_p)*beta_hmp*I_qmp*R_pr +
				  ((1-s_both)*(1-prot_t_m) + (1-s_p)*(1-prot_t_p) + (1-s_both)*joint_trt)*beta_hmp*I_qmp*Rt_pr  +
			         - 0.5*beta_hp*I_qp*I_hamp_pr +  
						+ 0.5*beta_hmp*I_qmp*I_hap_pr   
				    
		   	# Parameters for likelihoods (newly infected)
	
		dNew_n <- beta_hn*I_qn*S_h
		dNew_m <- beta_hm*I_qm*S_h
		dNew_p <- beta_hp*I_qp*S_h + pr_p*(1-s_p)*beta_hp*I_qp*S_pr 
		dNew_mp <- beta_hmp*I_qmp*S_h + pr_p*(1-s_p)*beta_hmp*I_qmp*S_pr  
				
		
	return(list(c(dS_q, dE_qn, dE_qm, dE_qp, dE_qmp, 
			 dI_qn, dI_qm, dI_qp,  dI_qmp, 
			 dS_h, dE_hn, dE_hm, dE_hp,  dE_hmp, 
			 dI_hsn, dI_hsm, dI_hsp, dI_hsmp, 
			 dRt, dR, dI_han, dI_ham, dI_hap,  dI_hamp, 
			 dNew_n, dNew_m, dNew_p, dNew_mp, 
			 dS_pr, dE_hp_pr, dE_hmp_pr,
			  dI_hsp_pr, dI_hsmp_pr,
			 dI_hap_pr, dI_hamp_pr,
			  dR_pr, dRt_pr)))}) 
	}
	
