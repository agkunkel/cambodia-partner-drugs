# cambodia-partner-drugs

This repository includes the code used to create the analyses in "Novel anti-malarial drug strategies to prevent artemisinin partner drug resistance: a model-based analysis". 

The entire process, from fitting model parameters to genearting output figures, can be run from the file "driver.R" which calls other files.

The files "malaria_IMIS_east.R" and "malaria_IMIS_west.R" contain the main functions involved in fitting the model for Eastern and Western Cambodia, respectively. It also relies on functions IMIS_mod_par.R, likelihoods.R, and model_fitting.R (which contains the actual differential equation model).

The files east_pars_0_0.csv and west_pars_0_0.csv contain parameter sets drawn from the posterior distributions of the model parameters after running the above mentioned codes.

The remaining R functions are involved in running the model for the final parameter sets and plotting the output.

The main folder shows baseline results; the folder new_s_m_s_p has new versions of s_m and s_p forces s_m=s_p for sensitivity analyses.

More details about the model and parameters can be found in the paper and appendix.

