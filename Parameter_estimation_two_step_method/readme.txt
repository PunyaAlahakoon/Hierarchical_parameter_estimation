This folder includes the two-step method to estimate parameters of an SIRS model. 
ONLY BETA IS CONSIDERED TO BE UNKNOWN, RATE OF RECOVERY=1, WANING IMMUNITY RATE=0.06. 

Step1_a === use this to independently, on parellel, to sample from posteriors
             using the ABC-SMC (Toni et al. version). 
		required functions needed to run this:
			Gillespie4, sample_gen and abc_ind (INCLUDED IN THE FOLDER)

****NOTE: you can use other software/ packages to do this step. 

Step1_b === use this to run the MCMC algorithm using the information in step1_a


Step 2== use hyper parameters to estimate sub-pop specific parameters 


**** MATLAB parallel computing toolbox is required (uses parfor function)