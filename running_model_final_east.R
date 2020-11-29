## This function is used to run the model after the parameter sets have already been selected
## It includes the intervention scenarios beginning in 2020
## It outputs CSV files used for plotting each intervention
## Also outputs directly plots showing model fit

running_model_final_east <- function(scenario, which_proph, fitness_sens_ana, fixed_inputs, variable_inputs_all) {
	
# which intervention scenario is being consideered
fixed_inputs["which_proph"] <- 0 # none
if (which_proph == "PPQ") {
	fixed_inputs["which_proph"] <- 1 # chemoprophylaxis
}
if (which_proph == "TRT") {
	fixed_inputs["which_proph"] <- 2 #triple ACT
}
	
# NOTE: slightly different terminology for the version of the model used for this code, than for the version used to fit
yinit <- c(S_q = 0.99, E_qn = 0, E_qm = 0, E_qp = 0, E_qmp = 0, 
			 I_qn = 0.01, I_qm = 0, I_qp = 0,  I_qmp = 0, 
			 S_h = NA, E_hn = 0, E_hm = 0, E_hp = 0,  E_hmp = 0, # To fill in based on fitting	 
			 I_hsn = 100, I_hsm = 0, I_hsp = 0,  I_hsmp = 0, 
			 Rt = 0, R = 0, I_han = 0, I_ham = 0, I_hap = 0, I_hamp = 0, 
			 New_n = 0, New_m = 0, New_p = 0, New_mp = 0,
			  S_pr = 0, E_hp_pr = 0, E_hmp_pr = 0,
			  I_hsp_pr = 0, I_hsmp_pr = 0,
			 I_hap_pr = 0, I_hamp_pr = 0,
			  R_pr = 0, Rt_pr =0)

### Setting time periods
start_time = 0 # start date (days). letting day 0 correspond to jan 1 1995
end_time_init = (2000-1995)*366
end_time_init2 = (2010-1995)*366 # 
proph_time = (2020-1995)*366 # end date (days)
end_time = (2030-1995)*366 # end date (days)
times <- seq(start_time, end_time, by = 30.5) #for simplicity, as if 12 equal months of 30.5 days 
times_init <- seq(start_time,end_time_init, by = 30.5)
times_init2 <- seq(end_time_init, end_time_init2, by = 30.5)
times_noproph <- seq(end_time_init2, proph_time, by = 30.5)
times_proph <- seq(proph_time, end_time, by = 30.5)

# Initializing data frames with 1 column per run
prev_pcr_df <- data.frame(matrix(0, ncol=dim(variable_inputs_all)[1]+1, nrow = length(times)))
rep_case_df <- data.frame(matrix(0, ncol=dim(variable_inputs_all)[1]+1, nrow = length(times)))
prev_pfmdr_df <- data.frame(matrix(0, ncol=dim(variable_inputs_all)[1]+1, nrow = length(times)))
prev_pfpm_df <- data.frame(matrix(0, ncol=dim(variable_inputs_all)[1]+1, nrow = length(times)))
prev_pfb_df <- data.frame(matrix(0, ncol=dim(variable_inputs_all)[1]+1, nrow = length(times)))
mres_prop_df <- data.frame(matrix(0, ncol=dim(variable_inputs_all)[1]+1, nrow = length(times)))
pres_prop_df <- data.frame(matrix(0, ncol=dim(variable_inputs_all)[1]+1, nrow = length(times)))
bres_prop_df <- data.frame(matrix(0, ncol=dim(variable_inputs_all)[1]+1, nrow = length(times)))
prev_mosq_df <- data.frame(matrix(0, ncol=dim(variable_inputs_all)[1]+1, nrow = length(times)))
prop_trans_asym_df <- data.frame(matrix(0, ncol=dim(variable_inputs_all)[1]+1, nrow = length(times)))
colnames(prev_pcr_df)[1] <- "year"
colnames(rep_case_df)[1] <- "year"
colnames(prev_pfmdr_df)[1] <- "year"
colnames(prev_pfpm_df)[1] <- "year"
colnames(prev_pfb_df)[1] <- "year"
colnames(mres_prop_df)[1] <- "year"
colnames(pres_prop_df)[1] <- "year"
colnames(bres_prop_df)[1] <- "year"
colnames(prev_mosq_df)[1] <- "year"
colnames(prop_trans_asym_df)[1] <- "year"
prev_pcr_df[,1] <- times/366+1995
rep_case_df[,1] <- times/366+1995
prev_pfmdr_df[,1] <- times/366+1995
prev_pfpm_df[,1] <- times/366+1995
prev_pfb_df[,1] <- times/366+1995
mres_prop_df[,1] <- times/366+1995
pres_prop_df[,1] <- times/366+1995
bres_prop_df[,1] <- times/366+1995
prev_mosq_df[,1] <- times/366+1995
prop_trans_asym_df[,1] <- times/366+1995

# loop through all variables
for (i in 1:dim(variable_inputs_all)[1]) {
	print(i)
	
	# defining parameter set and initial conditions for run
	variable_inputs <- variable_inputs_all[i,]
	if (fitness_sens_ana == F) {
		variable_inputs["f_mp"] <- variable_inputs["f_m"]*variable_inputs["f_p"] # standard assumption
	} else {
		variable_inputs["f_mp"] <- min(variable_inputs["f_p"], variable_inputs["f_m"]) - 0.01 # assumption for sensitivity analyses
	}
	pars <- data.frame(as.list(c(variable_inputs, fixed_inputs)))
	pars_vec <- pars # called by 2 diff names
    yinit["S_h"] <- pars$N_init-100

	##### model run - first time period
	results_init <- as.data.frame(ode(y=yinit, times= times_init, 
			func=malaria_model, pars, method='lsoda'))
	yinit2 <- subset(results_init,select=-time)[dim(results_init)[1],]
	ynames <- names(yinit2)
	yinit2 <- as.numeric(yinit2)
	names(yinit2) <- ynames
	
		##### model run - second time period
	results_init2 <- as.data.frame(ode(y=yinit2, times=times_init2, 
			func=malaria_model, pars, method='lsoda'))		
		yinit3 <- subset(results_init2,select=-time)[dim(results_init2)[1],]

	# Introducing MQ resistance
	# to mosquitoes
	yinit3$E_qm <- pars_vec$mdr_2000*yinit3$E_qn
	yinit3$E_qn <- (1-pars_vec$mdr_2000)*yinit3$E_qn
	yinit3$I_qm <- pars_vec$mdr_2000*yinit3$I_qn
	yinit3$I_qn <- (1-pars_vec$mdr_2000)*yinit3$I_qn
	# to humans
	yinit3$E_hm <- pars_vec$mdr_2000*yinit3$E_hn
	yinit3$E_hn <- (1-pars_vec$mdr_2000)*yinit3$E_hn
	yinit3$I_hsm <- pars_vec$mdr_2000*yinit3$I_hsn
	yinit3$I_hsn <- (1-pars_vec$mdr_2000)*yinit3$I_hsn
	yinit3$I_ham <- pars_vec$mdr_2000*yinit3$I_han
	yinit3$I_han <- (1-pars_vec$mdr_2000)*yinit3$I_han
	
	# introducing PPQ resistance
	ppq_init <- pars$ppq_init
	# to mosquitoes
	yinit3$E_qp <- with(yinit3, ppq_init*E_qn*E_qn/(E_qn+E_qm))
	yinit3$E_qmp <- with(yinit3, ppq_init*E_qn*E_qm/(E_qn+E_qm))
	yinit3$E_qn <- with(yinit3, (1-ppq_init)*E_qn)
	yinit3$I_qp <- with(yinit3, ppq_init*I_qn*I_qn/(I_qn+I_qm))
	yinit3$I_qmp <- with(yinit3, ppq_init*I_qn*I_qm/(I_qn+I_qm))
	yinit3$I_qn <- with(yinit3, (1-ppq_init)*I_qn)
	# to humans
	yinit3$E_hp <- with(yinit3, ppq_init*E_hn*E_hn/(E_hn+E_hm))
	yinit3$E_hmp <- with(yinit3, ppq_init*E_hn*E_hm/(E_hn+E_hm))
	yinit3$E_hn <- with(yinit3, (1-ppq_init)*E_hn)
	yinit3$I_hsp <- with(yinit3, ppq_init*I_hsn*I_hsn/(I_hsn+I_hsm))
	yinit3$I_hsmp <- with(yinit3, ppq_init*I_hsn*I_hsm/(I_hsn+I_hsm))
	yinit3$I_hsn <- with(yinit3, (1-ppq_init)*I_hsn)
	yinit3$I_hap <- with(yinit3, ppq_init*I_han*I_han/(I_han+I_ham))
	yinit3$I_hamp <- with(yinit3, ppq_init*I_han*I_ham/(I_han+I_ham))
	yinit3$I_han <- with(yinit3, (1-ppq_init)*I_han)
	
	ynames <- names(yinit3)
	yinit3 <- as.numeric(yinit3)
	names(yinit3) <- ynames	
			
	## Running for 3rd time period
	results1 <- as.data.frame(ode(y=yinit3, times= times_noproph, 
			func=malaria_model, pars, method='lsoda'))
	yinit4 <- results1[length(times_noproph),2:dim(results1)[2]]
	tmp <- names(yinit4)
	yinit4 <- as.numeric(yinit4)
	names(yinit4) <- tmp
	
	## Running through the end
	results2 <- as.data.frame(ode(y=yinit4, times= times_proph, 
			func=malaria_model, pars, method='lsoda'))
	results_init2 <- results_init2[2:dim(results_init2)[1],] 
	results1 <- results1[2:dim(results1)[1],] # repeated twice, removing once
	results2 <- results2[2:dim(results2)[1],] # 2019 repeated twice, removing once
	
	## Formatting output. 	
	
	results <- rbind(results_init, results_init2, results1, results2)
	results$year <- results$time/366+1995
	
	# size of human population at risk
	results$Nh <- with(results,{S_h + E_hn + E_hm + E_hp +  E_hmp + 
			 I_hsn + I_hsm + I_hsp +  I_hsmp +
			 Rt + R + I_han + I_ham + I_hap + I_hamp + 
			 S_pr +  E_hp_pr + E_hmp_pr  +
			 I_hsp_pr + I_hsmp_pr  + I_hap_pr + I_hamp_pr  + R_pr + Rt_pr})
		
		# mosquito population size		 
	results$Nm <- with(results,{S_q + E_qn + E_qm + E_qp + 
			E_qmp +  I_qn + I_qm + I_qp +  I_qmp })
					 
	# symptomatic and asymptomatic cases		 
	results$I_s <- with(results,{I_hsn + I_hsm + I_hsp + 
			 I_hsmp + I_hsp_pr + I_hsmp_pr })	
	results$I_a <- with(results,{ I_han + I_ham + I_hap +
		  I_hamp +  I_hap_pr + I_hamp_pr })	
	results$prop_trans_asym <- with(results,{ b*I_a/(I_s + b*I_a) })	 #transmission from asymptomatic
	results$prev_pcr <- with(results,{(I_s+I_a)/560000})  # PCR prevalence - whole population
	results$prev_mosq <- with(results,{(I_qn + I_qm + I_qp + I_qmp)/Nm}) # prevalence in mosquitoes

	# incident and reported incident human cases
	firstdiff <- function(x) {
	   shifted <- c(0,x[1:(length(x)-1)])
	   x-shifted
	   }
	results$New <- with(results,{New_n + New_m + New_p + New_mp})	
	results$rep <- pars$p*firstdiff(results$New) # reported cases
	
	# resistance - among new casses only
		results$I_q_inc <- with(results,{I_qn + pars$f_m*I_qm + pars$f_p*I_qp + pars$f_mp*I_qmp })
	results$prev_pfmdr <- with(results,{(pars$f_m*I_qm + pars$f_mp*I_qmp)/(I_q_inc)}) 
	 results$prev_pfpm <-  with(results,{(pars$f_p*I_qp + pars$f_mp*I_qmp)/(I_q_inc)}) 
	 results$prev_pfb <-  with(results,{(pars$f_mp*I_qmp)/(I_q_inc)}) 
	results$prev_mres <- (1-pars$s_m)*results$prev_pfmdr 
	results$prev_pres <- (1-pars$s_p)*results$prev_pfpm
	results$prev_bres <- (1-pars$s_both)*results$prev_pfb
			 
	# making data frames for plotting
	prev_pcr_df[,i+1] <- results$prev_pcr
	rep_case_df[,i+1] <- results$rep
	prev_pfmdr_df[,i+1] <- results$prev_pfmdr
	prev_pfpm_df[,i+1] <- results$prev_pfpm
	prev_pfb_df[,i+1] <- results$prev_pfb
	mres_prop_df[,i+1] <- results$prev_mres
	pres_prop_df[,i+1] <- results$prev_pres
	bres_prop_df[,i+1] <- results$prev_bres
	prev_mosq_df[,i+1] <- results$prev_mosq
	prop_trans_asym_df[,i+1] <- results$prop_trans_asym
	
}

if(scenario!="short_baseline")  {

# Creating output CSVs: incidence, prevalence, pfmdr1, pfpm2, mqres, ppqres
# prevalence
prev_pcr_out <- data.frame(matrix(0, ncol=5, nrow = length(times)))
colnames(prev_pcr_out) <- c('Year', 'Lower', 'Median', 'Mean', 'Upper')
prev_pcr_out$Year <- prev_pcr_df$year
prev_pcr_out$Lower <- apply(prev_pcr_df[,2:ncol(prev_pcr_df)], 1,  function(x) quantile(x, probs=0.025))
prev_pcr_out$Median <- apply(prev_pcr_df[,2:ncol(prev_pcr_df)], 1,  function(x) quantile(x, probs=0.5))
prev_pcr_out$Mean <- apply(prev_pcr_df[,2:ncol(prev_pcr_df)], 1,  mean)
prev_pcr_out$Upper <- apply(prev_pcr_df[,2:ncol(prev_pcr_df)], 1,  function(x) quantile(x, probs=0.975))
write.csv(prev_pcr_out, file= paste(scenario,which_proph,"east_prev_pcr.csv",sep="_"))

# incidence
rep_case_out <- data.frame(matrix(0, ncol=5, nrow = length(times)))
colnames(rep_case_out) <- c('Year', 'Lower', 'Median', 'Mean', 'Upper')
rep_case_out$Year <- rep_case_df$year
rep_case_out$Lower <- apply(rep_case_df[,2:ncol(rep_case_df)], 1,  function(x) quantile(x, probs=0.025))
rep_case_out$Median <- apply(rep_case_df[,2:ncol(rep_case_df)], 1,  function(x) quantile(x, probs=0.5))
rep_case_out$Mean <- apply(rep_case_df[,2:ncol(rep_case_df)], 1,  mean)
rep_case_out$Upper <- apply(rep_case_df[,2:ncol(rep_case_df)], 1,  function(x) quantile(x, probs=0.975))
write.csv(rep_case_out, file= paste(scenario,which_proph,"east_rep_case.csv",sep="_"))

# pfmdr1
prev_pfmdr_out <- data.frame(matrix(0, ncol=5, nrow = length(times)))
colnames(prev_pfmdr_out) <- c('Year', 'Lower', 'Median', 'Mean', 'Upper')
prev_pfmdr_out$Year <- prev_pfmdr_df$year
prev_pfmdr_out$Lower <- apply(prev_pfmdr_df[,2:ncol(prev_pfmdr_df)], 1,  function(x) quantile(x, probs=0.025))
prev_pfmdr_out$Median <- apply(prev_pfmdr_df[,2:ncol(prev_pfmdr_df)], 1,  function(x) quantile(x, probs=0.5))
prev_pfmdr_out$Mean <- apply(prev_pfmdr_df[,2:ncol(prev_pfmdr_df)], 1,  mean)
prev_pfmdr_out$Upper <- apply(prev_pfmdr_df[,2:ncol(prev_pfmdr_df)], 1,  function(x) quantile(x, probs=0.975))
write.csv(prev_pfmdr_out, file= paste(scenario,which_proph,"east_prev_pfmdr.csv",sep="_"))

# pfpm2
prev_pfpm_out <- data.frame(matrix(0, ncol=5, nrow = length(times)))
colnames(prev_pfpm_out) <- c('Year', 'Lower', 'Median', 'Mean', 'Upper')
prev_pfpm_out$Year <- prev_pfpm_df$year
prev_pfpm_out$Lower <- apply(prev_pfpm_df[,2:ncol(prev_pfpm_df)], 1,  function(x) quantile(x, probs=0.025))
prev_pfpm_out$Median <- apply(prev_pfpm_df[,2:ncol(prev_pfpm_df)], 1,  function(x) quantile(x, probs=0.5))
prev_pfpm_out$Mean <- apply(prev_pfpm_df[,2:ncol(prev_pfpm_df)], 1,  mean)
prev_pfpm_out$Upper <- apply(prev_pfpm_df[,2:ncol(prev_pfpm_df)], 1,  function(x) quantile(x, probs=0.975))
write.csv(prev_pfpm_out, file= paste(scenario,which_proph,"east_prev_pfpm.csv",sep="_"))

# pfmdr + pfpm2
prev_pfb_out <- data.frame(matrix(0, ncol=5, nrow = length(times)))
colnames(prev_pfb_out) <- c('Year', 'Lower', 'Median', 'Mean', 'Upper')
prev_pfb_out$Year <- prev_pfb_df$year
prev_pfb_out$Lower <- apply(prev_pfb_df[,2:ncol(prev_pfb_df)], 1,  function(x) quantile(x, probs=0.025))
prev_pfb_out$Median <- apply(prev_pfb_df[,2:ncol(prev_pfb_df)], 1,  function(x) quantile(x, probs=0.5))
prev_pfb_out$Mean <- apply(prev_pfb_df[,2:ncol(prev_pfb_df)], 1,  mean)
prev_pfb_out$Upper <- apply(prev_pfb_df[,2:ncol(prev_pfb_df)], 1,  function(x) quantile(x, probs=0.975))
write.csv(prev_pfb_out, file= paste(scenario,which_proph,"east_prev_pfb.csv",sep="_"))


# MQ res
mres_prop_out <- data.frame(matrix(0, ncol=5, nrow = length(times)))
colnames(mres_prop_out) <- c('Year', 'Lower', 'Median', 'Mean', 'Upper')
mres_prop_out$Year <- mres_prop_df$year
mres_prop_out$Lower <- apply(mres_prop_df[,2:ncol(mres_prop_df)], 1,  function(x) quantile(x, probs=0.025))
mres_prop_out$Median <- apply(mres_prop_df[,2:ncol(mres_prop_df)], 1,  function(x) quantile(x, probs=0.5))
mres_prop_out$Mean <- apply(mres_prop_df[,2:ncol(mres_prop_df)], 1,  mean)
mres_prop_out$Upper <- apply(mres_prop_df[,2:ncol(mres_prop_df)], 1,  function(x) quantile(x, probs=0.975))
write.csv(mres_prop_out, file= paste(scenario,which_proph,"east_mres_prop.csv",sep="_"))


# PPQ res
pres_prop_out <- data.frame(matrix(0, ncol=5, nrow = length(times)))
colnames(pres_prop_out) <- c('Year', 'Lower', 'Median', 'Mean', 'Upper')
pres_prop_out$Year <- pres_prop_df$year
pres_prop_out$Lower <- apply(pres_prop_df[,2:ncol(pres_prop_df)], 1,  function(x) quantile(x, probs=0.025))
pres_prop_out$Median <- apply(pres_prop_df[,2:ncol(pres_prop_df)], 1,  function(x) quantile(x, probs=0.5))
pres_prop_out$Mean <- apply(pres_prop_df[,2:ncol(pres_prop_df)], 1,  mean)
pres_prop_out$Upper <- apply(pres_prop_df[,2:ncol(pres_prop_df)], 1,  function(x) quantile(x, probs=0.975))
write.csv(pres_prop_out, file= paste(scenario,which_proph,"east_pres_prop.csv",sep="_"))

# MQ-PPQ Res
bres_prop_out <- data.frame(matrix(0, ncol=5, nrow = length(times)))
colnames(bres_prop_out) <- c('Year', 'Lower', 'Median', 'Mean', 'Upper')
bres_prop_out$Year <- bres_prop_df$year
bres_prop_out$Lower <- apply(bres_prop_df[,2:ncol(bres_prop_df)], 1,  function(x) quantile(x, probs=0.025))
bres_prop_out$Median <- apply(bres_prop_df[,2:ncol(bres_prop_df)], 1,  function(x) quantile(x, probs=0.5))
bres_prop_out$Mean <- apply(bres_prop_df[,2:ncol(bres_prop_df)], 1,  mean)
bres_prop_out$Upper <- apply(bres_prop_df[,2:ncol(bres_prop_df)], 1,  function(x) quantile(x, probs=0.975))
write.csv(bres_prop_out, file= paste(scenario,which_proph,"east_bres_prop.csv",sep="_"))

}

# ###########################################################################################
############################## PLOTTING LINES AND POINTS  #########################
###########################################################################################
# This is just to compare data and model pre-intervention
# run only for the short_baseline scenario

if(which_proph=="NONE" & scenario=="short_baseline") {
	
  pdf(file= paste(which_proph,"_east.pdf",sep=""))
  
  # to indicate which drug was used when
  drug_mat=as.data.frame(matrix(nrow=3,ncol=3))
  colnames(drug_mat) <- c('Start','End','Drug')
  drug_mat$Start <- as.numeric(c(2000, 2010, 2018))
drug_mat$End <- as.numeric(c(2010, 2018, Inf))
drug_mat$Drug <- c('ASMQ','DP','ASMQ')

######## Monthly case counts

rep_results <- melt(rep_case_df,id.vars="year")

dates <- seq(2004.25,2019.75,0.5)
cases <- c(227,1451,202,1398,439,2193,414,1149,299,1369,633,1857,311,2065,304,1436,278,790,180,503, 220,906, 209,861, 137,536, 245,1006, 187,770, 102,419)
cases_upper <- cases + 1.96*sqrt(cases)
cases_lower <- cases - 1.96*sqrt(cases)

results_tmp <- bind_rows(rep_results,data.frame(year=dates,cases = cases, cases_upper=cases_upper, cases_lower = cases_lower))

p <- ggplot(data = results_tmp) + geom_path(aes(x=year,y= value,group=variable)) + 
	geom_point(aes(x=year,y=cases),col='orangered', size=3) + 
	geom_rect(data=drug_mat, aes(xmin = Start, xmax = End,   ymin = -Inf, ymax = Inf, fill = Drug), alpha = 0.2 ) + 
  geom_text(aes(x = Start, y = 4500, label = Drug), data = drug_mat,  size = 3, vjust = 0, hjust = 0, nudge_x = 0) +
   ylim(0,5000) + ylab('Monthly cases') + theme_bw() +
   scale_fill_manual(values = c("gray80", "gray100"))  + scale_x_continuous(expand = c(0, 0), limits=c(2000,2020))
print(p)

# ######## Prevalence

pcr_results <- melt(prev_pcr_df,id.vars="year")

dates_prev <- c(2004.75,2007.75,2010.75,2013.75) # surveys took place in the fall
	
# East
gen_slide_N <- c(2885,3330,4505,4482)
gen_slide_pf <-c(128,114,25,2)
fg_slide_prop <- gen_slide_pf/gen_slide_N # this is about transforming from pcr to slide %
tmp <- 0.954 + 0.868*log(fg_slide_prop/(1-fg_slide_prop))
fg_pcr_prop <- exp(tmp)/(1+exp(tmp))
fg_pcr_num <- round(fg_pcr_prop*gen_slide_N)

pcr_results <- bind_rows(pcr_results,data.frame(year= dates_prev, fg_pcr_prop = fg_pcr_prop))

p <- ggplot(data = pcr_results) + geom_path(aes(x=year,y= value,group=variable)) + geom_point(aes(x=year,y= fg_pcr_prop),col='orangered', size=3) +  ylab('PCR Prevalence') + 
	geom_rect(data=drug_mat, aes(xmin = Start, xmax = End,   ymin = -Inf, ymax = Inf, fill = Drug), alpha = 0.2 ) + 
  geom_text(aes(x = Start, y = 0.23, label = Drug), data = drug_mat,  size = 3, vjust = 0, hjust = 0, nudge_x = 0) +
   ylim(0,0.25) +  theme_bw() +
   scale_fill_manual(values = c("gray80", "gray100"))  + scale_x_continuous(expand = c(0, 0), limits=c(2000,2020))
print(p)

# # ############ Plot pfmdr1 copy numbers - east

pfmdr_results <- melt(prev_pfmdr_df,id.vars="year")

dates <- c(2010.5, 2011.5, 2012.5, 2013.5, 2010.5, 2011.5, 2012.5, 2013.5, 2005, 
	 2017.5, 2015.5, 2016.5, 2016.5, 2017.5, 2018.5, 2017.5, 2017.5, 2018.5, 
	2016.5, 2017.5, 2015.5, 2016.5, 2018.5,
	2019.5, 2019.5) 

N_tested <- c(54, 100, 52, 34, 49, 72, 23, 13, 471, 247, 52, 44, 55, 59, 53, 68, 60, 44, 31, 149, 48, 14, 24,
          8,52 ) 
N_multicopies <-c(0, 9, 0, 1, 0, 1, 0, 0, 30, 1,  1, 15, 17, 0, 4, 1, 2, 0, 0, 2, 0,  0, 1,
          0, 26 ) 


pfmdr_multicopy_prop <- N_multicopies/N_tested
results_tmp <- bind_rows(pfmdr_results,data.frame(year= dates, pfmdr_multicopy_prop = pfmdr_multicopy_prop))

p <- ggplot(data = results_tmp) + geom_path(aes(x=year,y= value,group=variable)) + geom_point(aes(x=year,y= pfmdr_multicopy_prop),col='orangered', size=3) + ylim(0,1) + ylab('pfmdr1 multi-copies') + 
	geom_rect(data=drug_mat, aes(xmin = Start, xmax = End,   ymin = -Inf, ymax = Inf, fill = Drug), alpha = 0.2 ) + 
  geom_text(aes(x = Start, y = 0.85, label = Drug), data = drug_mat,  size = 3, vjust = 0, hjust = 0, nudge_x = 0) +
   scale_fill_manual(values = c("gray80", "gray100"))  + scale_x_continuous(expand = c(0, 0), limits=c(2000,2020)) +
   theme_bw()
print(p)

# # ############ Plot pfpm2 copy numbers - east

pfpm_results <- melt(prev_pfpm_df,id.vars="year")

dates <- c(2003, 2005, 2007, 2011, 2011, 2013, 2013, 2013, 2015, 2015, 2015, 2010.5, 
	2011.5, 2012.5, 2013.5, 2017.5, 2015.5, 2016.5, 2016.5, 2017.5, 
	2018.5, 2017.5, 2016.5, 2017.5, 2018.5, 2017.5, 2015.5, 2016.5, 2018.5,
	2019.5,2019.5) 

N_tested <- c(53, 55, 67, 30, 51, 31, 38, 22, 88, 55, 46, 49, 72, 23, 13, 247, 52, 45, 53, 
	59, 53, 68, 31, 60, 44, 149, 48, 15, 24,
	8, 52) 

N_multicopies <-c(1, 3, 1, 0, 1, 1, 3, 1, 40, 14, 29, 0, 1, 0, 2, 153, 26, 7, 30, 39, 8, 43, 4, 34, 7, 101, 36, 9, 6,
	1, 2) 

pfpm_multicopy_prop <- N_multicopies/N_tested
results_tmp <- bind_rows(pfpm_results,data.frame(year= dates, pfpm_multicopy_prop = pfpm_multicopy_prop))

p <- ggplot(data = results_tmp) + geom_path(aes(x=year,y= value,group=variable)) + geom_point(aes(x=year,y= pfpm_multicopy_prop),col='orangered', size=3)   + ylim(0,1) + ylab('pfpm2 multi-copies') + 
	geom_rect(data=drug_mat, aes(xmin = Start, xmax = End,   ymin = -Inf, ymax = Inf, fill = Drug), alpha = 0.2 ) + 
  geom_text(aes(x = Start, y = 0.85, label = Drug), data = drug_mat,  size = 3, vjust = 0, hjust = 0, nudge_x = 0) +
   scale_fill_manual(values = c("gray80", "gray100"))  + scale_x_continuous(expand = c(0, 0), limits=c(2000,2020)) +
   theme_bw()
print(p)

# # ############ Plot MQ TES - east

mres_results <- melt(mres_prop_df,id.vars="year")

dates <- c(2001.5, 2002.5, 2003.5, 2004.5, 2017.5, 2014.5, 2015.5, 2016.5, 
	2018.5, 2018.5, 2018.5, 2018.5, 2018.5,
	2019.5, 2019.5) 
	
N_tested <- c(50, 71, 63, 80, 59, 91, 54, 19, 29, 53, 6, 44, 24,
	8, 50) 

N_LTF <-c(0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0,
	0, 0) 

mres_prop <- N_LTF/N_tested
results_tmp <- bind_rows(mres_results,data.frame(year= dates, mres_prop = mres_prop))

p <- ggplot(data = results_tmp) + geom_path(aes(x=year,y= value,group=variable)) + geom_point(aes(x=year,y= mres_prop),col='orangered', size=3) + ylim(0,1) + ylab('Phenotypic MQ resistance')  + 
	geom_rect(data=drug_mat, aes(xmin = Start, xmax = End,   ymin = -Inf, ymax = Inf, fill = Drug), alpha = 0.2 ) + 
  geom_text(aes(x = Start, y = 0.85, label = Drug), data = drug_mat,  size = 3, vjust = 0, hjust = 0, nudge_x = 0) +
   scale_fill_manual(values = c("gray80", "gray100"))  + scale_x_continuous(expand = c(0, 0), limits=c(2000,2020)) +
   theme_bw()

print(p)

 # ############ Plot PPQ TES - east

pres_results <- melt(pres_prop_df,id.vars="year")

dates <- c(2012.5, 2011.5, 2013.5, 2014.5, 2014.5, 2015.5, 2016.5, 2010.5, 2013.5, 
		2016.5, 2013.5, 2009.5)
	
N_tested <- c(60, 55, 22, 40, 32, 60, 45, 59, 60, 58, 60, 55)
	
N_LTF <-c(0, 2, 0, 4, 11, 5, 12, 0, 1, 8, 1, 0)

pres_prop <- N_LTF/N_tested
results_tmp <- bind_rows(pres_results,data.frame(year= dates, pres_prop = pres_prop))

p <- ggplot(data = results_tmp) + geom_path(aes(x=year,y= value,group=variable)) + geom_point(aes(x=year,y= pres_prop),col='orangered', size=3) + ylim(0,1) + ylab('Phenotypic PPQ resistance') + 
	geom_rect(data=drug_mat, aes(xmin = Start, xmax = End,   ymin = -Inf, ymax = Inf, fill = Drug), alpha = 0.2 ) + 
  geom_text(aes(x = Start, y = 0.85, label = Drug), data = drug_mat,  size = 3, vjust = 0, hjust = 0, nudge_x = 0) +
   scale_fill_manual(values = c("gray80", "gray100"))  + scale_x_continuous(expand = c(0, 0), limits=c(2000,2020)) +
   theme_bw()

print(p)

dev.off()	
}
}