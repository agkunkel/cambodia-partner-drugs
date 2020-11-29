# Model for fitting resistance data (smaller version with no interventions)
# Includes: no resistance, MQ only, PPQ only, MQ+PPQ
# unlike for the "scenarios" model, a lowercase "m" is used for mosquitoes and "b" for resistance to both MQ and PPQ

malaria_model <- function(times,yinit,pars){
	with(as.list(c(yinit,pars)), {
		
	# size of human population at risk
		Nh = S_h + E_hn + E_hm + E_hp + E_hb +
			 I_hsn + I_hsm + I_hsp + I_hsb +
			 Rt + R + I_han + I_ham + I_hap + I_hab
			 
			 # size of mosquito population
		Nm = S_m + E_mn + E_mm + E_mp + E_mb +
			 I_mn + I_mm + I_mp + I_mb  
	
		# year
		year = times/366+1995 
		
		# proportion receiving appropriate treatment increases from 2004-2010
		if (year < 2004) {
			app <- 0.3
		} else if (year < 2010) {
			app <- 0.3 + (0.8 - 0.3)/(2010-2004)*(year-2004)
		} else {
			app <- 0.8
		}
		
		r_sn = r_s_0*app + r_s_0*(1-app)*s_a # successful treatment if no drug resistance
		rf_sn = r_s_0*(1-app)*(1-s_a) # failed treatment if no drug resistance (due to inappropriate treatment)
		acq_rate = r_s_0*a_0*app/(1-a_0) # acquire resistance if no drug resistance and appropriate treatment
		
		
		# Changes in first line drug regimen over time
		if (year < 2000) { # Before the introduction of ACTs
			
			r_sm <- r_sn # successful treatment
			r_sp <- r_sn
			r_sb <- r_sn
			rf_sm = rf_sn # failed treatment
			rf_sp = rf_sn 
			rf_sb = rf_sn 
			a_m = 0   #acquired resistance - none
			a_p = 0
			prot_t_m = 1 #post-treatment, protected against MQ and PPQ resistant strains
			prot_t_p = 1
					
		} else if (year < 2008) { # using ASMQ in whole country
			
			r_sm <- r_s_0*app*s_m + r_s_0*(1-app)*s_a  #successful treatment if genotyipc MQ resistance 
			r_sp <- r_sn
			r_sb <- r_s_0*app*s_m + r_s_0*(1-app)*s_a
			
			rf_sm = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)   # failed treatment - due to resistance or inappropriate treatment
			rf_sp = rf_sn 
			rf_sb = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)  
			
			a_m = acq_rate    #acquired resistance - only on appropriate treatment
			if (region ==2) { # Eastern Cambodia, introduce MQ resistance to model in 2010 (see Appendix)
				a_m = 0 } 	
			a_p = 0
			
			prot_t_m = 0 # post-treatment protection against PPQ-resistant but not MQ-resistant strains
			prot_t_p = 1
			
		} else if (year < 2010) { # Policy varies by region
			
			if (region ==1) { # west # using DP in western Cambodia
				r_sm <- r_sn # successsful treatment
				r_sp <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
				r_sb <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
				rf_sm = rf_sn # failed treatment
			   rf_sp =  r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a) 
			   rf_sb = r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a)
				a_m = 0   # acquired resistance
				a_p = acq_rate
				prot_t_m = 1 # post-treatment protection against MQ rsistant strains
				prot_t_p = 0
			} else if (region ==2) { #east # still using ASMQ in eastern Cambodia
				r_sm <- r_s_0*app*s_m + r_s_0*(1-app)*s_a  # successful treatment
				r_sp <- r_sn
				r_sb <- r_s_0*app*s_m + r_s_0*(1-app)*s_a
				rf_sm = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)  # failed treatment
		    	rf_sp = rf_sn 
			    rf_sb = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)  
				# acquired resistance
				 a_m = 0 # acquired MQ resistance introduced in model in Eastern Cambodia in 2010
				a_p = 0
				prot_t_m = 0
				prot_t_p = 1
			}	
			
		}  else if (year < 2017) { # Using DP in whole country 
			r_sm <- r_sn # successful treatment
			r_sp <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
			r_sb <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
			rf_sm = rf_sn # failed treatment
			rf_sp = r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a) 
			rf_sb = r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a) 
			a_m = 0  #acquired resistance
			a_p = acq_rate
			prot_t_m = 1 #post-treatment protection
			prot_t_p = 0
			
		} else if (year < 2018) {  # Switch back to ASMQ occurs first in Western Cambodia
			
			if (region ==1) { # west # using MQ in western Cambodia
				r_sm <- r_s_0*app*s_m + r_s_0*(1-app)*s_a # successful treatment
				r_sp <- r_sn
				r_sb <- r_s_0*app*s_m + r_s_0*(1-app)*s_a
				rf_sm = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)  #  failed treatment
		    	rf_sp = rf_sn 
			    rf_sb = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a)  
				a_m = acq_rate # acquired resistance
				a_p = 0
				prot_t_m = 0
				prot_t_p = 1
			} else if (region ==2) { #east # still using DP in eastern Cambodia
				r_sm <- r_sn  #successful treatment
				r_sp <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
				r_sb <- r_s_0*app*s_p + r_s_0*(1-app)*s_a
				rf_sm = rf_sn # NEW - failed treatment
			   rf_sp = r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a) 
			   rf_sb = r_s_0*app*(1-s_p) + r_s_0*(1-app)*(1-s_a) 
				a_m = 0   # acquired resistance
				a_p = acq_rate
				prot_t_m = 1
				prot_t_p = 0
			}	
			
		 } else { # using ASMQ in whole country
			
			r_sm <- r_s_0*app*s_m + r_s_0*(1-app)*s_a # successful treatment
			r_sp <- r_sn
			r_sb <- r_s_0*app*s_m + r_s_0*(1-app)*s_a
			
			rf_sm = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a) #  failed treatment
			rf_sp = rf_sn 
			rf_sb = r_s_0*app*(1-s_m) + r_s_0*(1-app)*(1-s_a) 
			
			a_m = acq_rate    # acquired resistance
			a_p = 0
			prot_t_m = 0
			prot_t_p = 1
		}
			
		
		# population exit rates (declining population at risk)
		# input, mu_0, is in terms of days similar to other model inputs
		# 2010-2013: decline from about 7.5% to 5% in both E and W = 0.67
		# 2004 - 2007: decline from about 20% to 15% overall (not given separately) = 0.75
		# 2007-2010: question appears to have changed from work in forest to spend night in forest
		start_decl = 2004 
		lambda = log(0.7)/((2013 - 2010)*366) # a compromise
		if (year < start_decl) {
			mu = mu_0
		} else if (year < 2018) {
			mu = mu_0 - lambda # note lambda is negative
		} else { # fix from 2018 onwards
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
		beta_m = mosq_a_t*mosq_c
		beta_h = mosq_m*mosq_a_t*mosq_b
		
		# trying to force trough near april, peak near september through 
	    amp <- (1 - a*sin(2*pi*(times)/366)) 	 
			
		# arbitrarily, applying fitness costs and seasonality only one direction (mosquito -> human)	
				
		beta_mn <- beta_m # human infecting mosquito doesn't change over time, let's say
		beta_mm <- beta_mn
		beta_mp <- beta_mn
		beta_mb <- beta_mn
		
		beta_hn <- beta_h*amp # accounting for seasonality
		beta_hm <- beta_hn*f_m # and fitness costs
		beta_hp <- beta_hn*f_p
		beta_hb <- beta_hn*f_b
		
		I_hn <- I_hsn + b*I_han # infectiousness from symptomatic + asymptomatic
		I_hm <- I_hsm + b*I_ham
		I_hp <- I_hsp + b*I_hap
		I_hb <- I_hsb + b*I_hab
		
		####################### MOSQUITO EQUATIONS
		
		# Holding population sizes fixed and letting beta vary for seasonality
		# (Considering only % infected, not raw numbers)
		
		# Susceptible mosquitoes 
		dS_m <- mu_m*(Nm - S_m) +
			- beta_mn*I_hn/Nh*S_m - beta_mm*I_hm/Nh*S_m +
			- beta_mp*I_hp/Nh*S_m - beta_mb*I_hb/Nh*S_m
		
		# Exposed mosquitoes
		dE_mn <- beta_mn*I_hn/Nh*S_m - (mu_m + p_m)*E_mn # no resistance
		
		dE_mm <- beta_mm*I_hm/Nh*S_m - (mu_m + p_m)*E_mm # genotypic MQ resistance
		
		dE_mp <- beta_mp*I_hp/Nh*S_m - (mu_m + p_m)*E_mp # genotypic PPQ resistance
		
		dE_mb <- beta_mb*I_hb/Nh*S_m - (mu_m + p_m)*E_mb # genotypic MQ and PPQ resistance ("both")
		
		# Infectious mosquitoes
		dI_mn <- p_m*E_mn - mu_m*I_mn
		
		dI_mm <- p_m*E_mm - mu_m*I_mm
		
		dI_mp <- p_m*E_mp - mu_m*I_mp
		
		dI_mb <- p_m*E_mb - mu_m*I_mb
		
		
		###################### HUMAN EQUATIONS
		
		# Susceptible humans
		dS_h <- mu_0*Nh +
			- beta_hn*I_mn*S_h - beta_hm*I_mm*S_h - beta_hp*I_mp*S_h - beta_hb*I_mb*S_h +
			+ w*R - mu*S_h
			
		# Exposed humans
			
		dE_hn <- beta_hn*I_mn*S_h - (mu + p_h)*E_hn
		 
		dE_hm <- beta_hm*I_mm*S_h - (mu + p_h)*E_hm
		
		dE_hp <- beta_hp*I_mp*S_h - (mu + p_h)*E_hp
		
		dE_hb <- beta_hb*I_mb*S_h - (mu + p_h)*E_hb
		
		# Infectious, symptomatic humans
		
		dI_hsn <-  p_h*E_hn - mu*I_hsn +
					- r_sn*I_hsn + - 0.1*rf_sn*I_hsn +  
					- a_m*I_hsn - a_p*I_hsn 
		
		dI_hsm <-  p_h*E_hm - mu*I_hsm +
					- r_sm*I_hsm + - 0.1*rf_sm*I_hsm +
					- a_p*I_hsm + a_m*I_hsn 
		
		dI_hsp <-  p_h*E_hp - mu*I_hsp +
					- r_sp*I_hsp + - 0.1*rf_sp*I_hsp +
					- a_m*I_hsp + a_p*I_hsn			
		
		dI_hsb <-  p_h*E_hb - mu*I_hsb +
					- r_sb*I_hsb + - 0.1*rf_sb*I_hsb +
					a_m*I_hsp + a_p*I_hsm 
		

		# protected by treatment drug

		dRt <- app*(r_sn*I_hsn + r_sm*I_hsm + r_sp*I_hsp + r_sb*I_hsb) +
			- (1-s_m)*(1-prot_t_m)*beta_hm*I_mm*Rt - (1-s_p)*(1-prot_t_p)*beta_hp*I_mp*Rt +
			- ((1-s_m)*(1-prot_t_m) + (1-s_p)*(1-prot_t_p))*beta_hb*I_mb*Rt +
			- w_t*Rt - mu*Rt 
			
		# partially immune, no longer protected by treatment drug
		
		dR <- (1-app)*(r_sn*I_hsn + r_sm*I_hsm + r_sp*I_hsp + r_sb*I_hsb) +
			w_t*Rt + r_a*(I_han + I_ham + I_hap + I_hab) +
			- (beta_hn*I_mn + beta_hm*I_mm)*R + 
			- (beta_hp*I_mp + beta_hb*I_mb)*R +
			- w*R - mu*R
		
		# infectious, asymptomatic 
		
		rc = 0.05 # recombination coefficient
		
		dI_han <- beta_hn*I_mn*R - r_a*I_han - mu*I_han +  0.1*rf_sn*I_hsn +
					+ 0.5*beta_hn*I_mn*(I_ham + I_hap) + (0.5 - rc)*beta_hn*I_mn*I_hab + 
					- 0.5*(beta_hm*I_mm + beta_hp*I_mp)*I_han - (0.5 + rc)*beta_hb*I_mb*I_han + 
					+ rc*beta_hm*I_mm*I_hap + rc*beta_hp*I_mp*I_ham 
		
		dI_ham <- beta_hm*I_mm*R + (1-prot_t_m)*beta_hm*I_mm*Rt + 0.1*rf_sm*I_hsm +
			- r_a*I_ham - mu*I_ham +
			+ 0.5*beta_hm*I_mm*(I_han + I_hab) + (0.5 - rc)*beta_hm*I_mm*I_hap + 
			- 0.5*(beta_hn*I_mn + beta_hb*I_mb)*I_ham - (0.5 + rc)*beta_hp*I_mp*I_ham  + 
			+ rc*beta_hn*I_mn*I_hab + rc*beta_hb*I_mb*I_han 
		
		dI_hap <- beta_hp*I_mp*R + (1-prot_t_p)*beta_hp*I_mp*Rt + 0.1*rf_sp*I_hsp +
			- r_a*I_hap - mu*I_hap  +
			+ 0.5*beta_hp*I_mp*(I_han + I_hab) + (0.5 - rc)*beta_hp*I_mp*I_ham + 
			- 0.5*(beta_hn*I_mn + beta_hb*I_mb)*I_hap - (0.5 + rc)*beta_hm*I_mm*I_hap + 
			+ rc*beta_hn*I_mn*I_hab + rc*beta_hb*I_mb*I_han 
		
		dI_hab <- beta_hb*I_mb*R + 0.1*rf_sb*I_hsb +
		((1-s_m)*(1-prot_t_m) + (1-s_p)*(1-prot_t_p))*beta_hb*I_mb*Rt +
		    - r_a*I_hab - mu*I_hab  +
			+ 0.5*beta_hb*I_mb*(I_ham + I_hap) + (0.5 - rc)*beta_hb*I_mb*I_han + 
			- 0.5*(beta_hm*I_mm + beta_hp*I_mp)*I_hab - (0.5 + rc)*beta_hn*I_mn*I_hab + 
			+ rc*beta_hm*I_mm*I_hap + rc*beta_hp*I_mp*I_ham 

	# Parameters for likelihoods (new infections)
	
		dNew_n <- beta_hn*I_mn*S_h
		dNew_m <- beta_hm*I_mm*S_h
		dNew_p <- beta_hp*I_mp*S_h
		dNew_b <- beta_hb*I_mb*S_h
		
	
	return(list(c(dS_m, dE_mn, dE_mm, dE_mp, dE_mb,
			 dI_mn, dI_mm, dI_mp, dI_mb,
			 dS_h, dE_hn, dE_hm, dE_hp, dE_hb,
			 dI_hsn, dI_hsm, dI_hsp, dI_hsb,
			 dRt, dR, dI_han, dI_ham, dI_hap, dI_hab,
			 dNew_n, dNew_m,dNew_p,dNew_b)))}) 
	}
	
