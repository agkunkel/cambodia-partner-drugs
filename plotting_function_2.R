########### This function formats the plots seen in the main text of the paper

plotting_function <- function(scenario) {
source('multiplot.R')	
	
inc_none_east <- read.csv(paste(scenario,'_NONE_east_rep_case.csv',sep=""))
inc_none_east$Strategy <- 'Status Quo'
inc_none_east$Region <- 'East'
inc_none_east$Indicator <- 'Number of cases'
inc_ppq_east <- read.csv(paste(scenario,'_PPQ_east_rep_case.csv',sep=""))
inc_ppq_east$Strategy <- 'Piperaquine Prophylaxis'
inc_ppq_east $Region <- 'East'
inc_ppq_east $Indicator <- 'Number of cases'
inc_trt_east <- read.csv(paste(scenario,'_TRT_east_rep_case.csv',sep=""))
inc_trt_east$Strategy <- 'Triple ACT'
inc_trt_east $Region <- 'East'
inc_trt_east $Indicator <- 'Number of cases'
inc_east <- bind_rows(inc_none_east, inc_ppq_east, inc_trt_east)
	

inc_none_west <- read.csv(paste(scenario,'_NONE_west_rep_case.csv',sep=""))
inc_none_west$Strategy <- 'Status Quo'
inc_none_west$Region <- 'West'
inc_none_west$Indicator <- 'Number of cases'
inc_ppq_west <- read.csv(paste(scenario,'_PPQ_west_rep_case.csv',sep=""))
inc_ppq_west$Strategy <- 'Piperaquine Prophylaxis'
inc_ppq_west$Region <- 'West'
inc_ppq_west$Indicator <- 'Number of cases'
inc_trt_west <- read.csv(paste(scenario,'_TRT_west_rep_case.csv',sep=""))
inc_trt_west$Strategy <- 'Triple ACT'
inc_trt_west$Region <- 'West'
inc_trt_west$Indicator <- 'Number of cases'

inc_west <- bind_rows(inc_none_west, inc_ppq_west, inc_trt_west)
	
	
	# pfmdr1
pfmdr_none_east <- read.csv(paste(scenario,'_NONE_east_prev_pfmdr.csv',sep=""))
pfmdr_none_east$Strategy <- 'Status Quo'
pfmdr_none_east $Region <- 'East'
pfmdr_none_east $Indicator <- 'pfmdr1'
pfmdr_ppq_east <- read.csv(paste(scenario,'_PPQ_east_prev_pfmdr.csv',sep=""))
pfmdr_ppq_east$Strategy <- 'Piperaquine Prophylaxis'
pfmdr_ppq_east $Region <- 'East'
pfmdr_ppq_east $Indicator <- 'pfmdr1'
pfmdr_trt_east <- read.csv(paste(scenario,'_TRT_east_prev_pfmdr.csv',sep=""))
pfmdr_trt_east$Strategy <- 'Triple ACT'
pfmdr_trt_east $Region <- 'East'
pfmdr_trt_east $Indicator <- 'pfmdr1'

pfmdr_east <- bind_rows(pfmdr_none_east, pfmdr_ppq_east, pfmdr_trt_east)


pfmdr_none_west <- read.csv(paste(scenario,'_NONE_west_prev_pfmdr.csv',sep=""))
pfmdr_none_west$Strategy <- 'Status Quo'
pfmdr_none_west $Region <- 'West'
pfmdr_none_west $Indicator <- 'pfmdr1'
pfmdr_ppq_west <- read.csv(paste(scenario,'_PPQ_west_prev_pfmdr.csv',sep=""))
pfmdr_ppq_west$Strategy <- 'Piperaquine Prophylaxis'
pfmdr_ppq_west $Region <- 'West'
pfmdr_ppq_west $Indicator <- 'pfmdr1'
pfmdr_trt_west <- read.csv(paste(scenario,'_TRT_west_prev_pfmdr.csv',sep=""))
pfmdr_trt_west$Strategy <- 'Triple ACT'
pfmdr_trt_west $Region <- 'West'
pfmdr_trt_west $Indicator <- 'pfmdr1'

pfmdr_west <- bind_rows(pfmdr_none_west, pfmdr_ppq_west, pfmdr_trt_west)

	
	
	# pfpm2
pfpm_none_east <- read.csv(paste(scenario,'_NONE_east_prev_pfpm.csv',sep=""))
pfpm_none_east$Strategy <- 'Status Quo'
pfpm_none_east $Region <- 'East'
pfpm_none_east $Indicator <- 'pfpm2'
pfpm_ppq_east <- read.csv(paste(scenario,'_PPQ_east_prev_pfpm.csv',sep=""))
pfpm_ppq_east$Strategy <- 'Piperaquine Prophylaxis'
pfpm_ppq_east $Region <- 'East'
pfpm_ppq_east $Indicator <- 'pfpm2'
pfpm_trt_east <- read.csv(paste(scenario,'_TRT_east_prev_pfpm.csv',sep=""))
pfpm_trt_east$Strategy <- 'Triple ACT'
pfpm_trt_east $Region <- 'East'
pfpm_trt_east $Indicator <- 'pfpm2'

pfpm_east <- bind_rows(pfpm_none_east, pfpm_ppq_east, pfpm_trt_east)


pfpm_none_west <- read.csv(paste(scenario,'_NONE_west_prev_pfpm.csv',sep=""))
pfpm_none_west$Strategy <- 'Status Quo'
pfpm_none_west $Region <- 'West'
pfpm_none_west $Indicator <- 'pfpm2'
pfpm_ppq_west <- read.csv(paste(scenario,'_PPQ_west_prev_pfpm.csv',sep=""))
pfpm_ppq_west$Strategy <- 'Piperaquine Prophylaxis'
pfpm_ppq_west $Region <- 'West'
pfpm_ppq_west $Indicator <- 'pfpm2'
pfpm_trt_west <- read.csv(paste(scenario,'_TRT_west_prev_pfpm.csv',sep=""))
pfpm_trt_west$Strategy <- 'Triple ACT'
pfpm_trt_west $Region <- 'West'
pfpm_trt_west $Indicator <- 'pfpm2'

pfpm_west <- bind_rows(pfpm_none_west, pfpm_ppq_west, pfpm_trt_west)

	
	
	# pfmdr1 + pfpm2
pfb_none_east <- read.csv(paste(scenario,'_NONE_east_prev_pfb.csv',sep=""))
pfb_none_east$Strategy <- 'Status Quo'
pfb_none_east $Region <- 'East'
pfb_none_east $Indicator <- 'pfpm2 + pfmdr1'
pfb_ppq_east <- read.csv(paste(scenario,'_PPQ_east_prev_pfb.csv',sep=""))
pfb_ppq_east$Strategy <- 'Piperaquine Prophylaxis'
pfb_ppq_east $Region <- 'East'
pfb_ppq_east $Indicator <- 'pfpm2 + pfmdr1'
pfb_trt_east <- read.csv(paste(scenario,'_TRT_east_prev_pfb.csv',sep=""))
pfb_trt_east$Strategy <- 'Triple ACT'
pfb_trt_east $Region <- 'East'
pfb_trt_east $Indicator <- 'pfpm2 + pfmdr1'

pfb_east <- bind_rows(pfb_none_east, pfb_ppq_east, pfb_trt_east)


pfb_none_west <- read.csv(paste(scenario,'_NONE_west_prev_pfb.csv',sep=""))
pfb_none_west$Strategy <- 'Status Quo'
pfb_none_west $Region <- 'West'
pfb_none_west $Indicator <- 'pfpm2 + pfmdr1'
pfb_ppq_west <- read.csv(paste(scenario,'_PPQ_west_prev_pfb.csv',sep=""))
pfb_ppq_west$Strategy <- 'Piperaquine Prophylaxis'
pfb_ppq_west $Region <- 'West'
pfb_ppq_west $Indicator <- 'pfpm2 + pfmdr1'
pfb_trt_west <- read.csv(paste(scenario,'_TRT_west_prev_pfb.csv',sep=""))
pfb_trt_west$Strategy <- 'Triple ACT'
pfb_trt_west $Region <- 'West'
pfb_trt_west $Indicator <- 'pfpm2 + pfmdr1'

pfb_west <- bind_rows(pfb_none_west, pfb_ppq_west, pfb_trt_west)

p_all <- bind_rows(inc_east, inc_west, pfmdr_east, pfmdr_west, pfpm_east, pfpm_west, pfb_east, pfb_west)
p_all <- subset(p_all, Year>2015)
p_all$Strategy <- factor(p_all$Strategy, levels=c("Triple ACT","Piperaquine Prophylaxis","Status Quo"))

blank_data <- data.frame(Indicator = c("Number of cases", "Number of cases", "Number of cases", "Number of cases", 'pfmdr1', 'pfmdr1', 'pfmdr1', 'pfmdr1', 'pfpm2', 'pfpm2', 'pfpm2', 'pfpm2', 'pfpm2 + pfmdr1', 'pfpm2 + pfmdr1', 'pfpm2 + pfmdr1', 'pfpm2 + pfmdr1'),  Region = c('East','East','West','West','East','East','West','West','East','East','West','West','East','East','West','West'), x = 2020, 
y = c(0, 1200, 0, 1200, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1), lab=c('','A','','B','','C','','D','','E','','F','','G','','H'))
a <- ggplot(data=p_all) + xlim(2018,2025) + geom_ribbon(aes(x=Year, ymin=Lower, ymax=Upper, fill=Strategy), alpha=0.5)  + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_grid(Indicator~Region, scales="free") + geom_line(aes(x=Year, y=Median, col=Strategy)) + geom_vline(xintercept=2020, size=0.3) + geom_vline(xintercept=2021,linetype="dashed", size=0.3)  + geom_text(data=blank_data, aes(label=lab, y=y),x=2024)+ theme_bw() 
ggsave(plot=a, filename=paste(scenario,'.pdf',sep=""))

	}