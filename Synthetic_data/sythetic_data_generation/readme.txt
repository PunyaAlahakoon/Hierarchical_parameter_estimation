THIS FLODER CONTAINS THE CODES THAT GENERATED SYNTHETIC DATA FOR 
THE 2 POPULATIONS EACH HAVINH 15 SUB-POPS EACH. 



popA_data= synthetic data generation for Population A 
popB_data=  synthetic data generation for Population B

how to generate data: 

1) Scenario 1: identical model parameters with stochastic effects  
beta = 2; 
gamma =1;
mu=0.06;
k=15; %number of pops 

2) Scenario 2: variable model parameters with stochastic effects
pd = makedist('Normal','mu',2,'sigma',0.5);
tpd=truncate(pd,1,10);

mean(tpd) = 2.0276
sigma=std(tpd) = 0.4708  
t_betas=random(tpd,1,k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std(t_betas)=0.4821 %in the data that were generated
mean(t_betas)= 2.1805 %in the data that were generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NOTE: sub-pop 16 $\beta$ was also generated similar to that of Scenario 2 above.