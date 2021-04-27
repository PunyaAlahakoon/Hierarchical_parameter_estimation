%Step (1)--b
%use this to estimate the hyper-parameters (ppsi_beta, sigma_beta) for an SIRS model 
%kept epsilon and gamma fixed 
%add a unif prior for hyper-parameter  
%functions need to run this:
        %weight3 --CALCULATES THE ESTIMATED LIKELIHOOD---MAKE SURE TO CHANGE
        %THE CONDITIONAL PRIOR IF REQUIRED 

%load parameter sets of all the populations as matrices
me=beta_smc; %from Step (1)--a  

%N=length(gamma);
B=size(me,1); %length of particles in Step (1)--a
dim=size(me,2); %number of sub-pops 

N=10000; %number of MCMC iterations 


%store hyper-parameters 
hb=zeros(1,N); %mean hyper parameters for beta ----psi_beta
hb_sig=zeros(1,N);%std hyper-parameteer  ---- sigma_beta
%store probabilities of accepting 
p=zeros(1,N); %not a necessary step

%MCMC initialisation 
hb(1)= 3;
hb_sig(1)=0.1;
%calculate the weight: 
[w,poste]=weight3(me,hb(1),hb_sig(1),dim,B);
poste2=poste;
p(1)=w*unifpdf(hb(1),0.001,10)*unifpdf(hb_sig(1),0,2.5);

%MCMC 
for i=2:N
   %propose values from beta and gamma matrices as vectors 
   hm_r=abs(normrnd(hb(i-1),0.1));
   sig_r=abs(normrnd(hb_sig(i-1),0.1));
    
   [w,poste]=weight3(me,hm_r,sig_r,dim,B);
   prior=unifpdf(hm_r,0.001,10)*unifpdf(hb_sig(1),0,2.5);
   
   p_s=w*prior;
   prop=p_s/p(i-1);

   alpha=min(1,prop);
   u=rand(1);
   
   if u<alpha
       hb(i)=hm_r;
       hb_sig(i)=sig_r;
       p(i)=p_s;
       poste2=poste;
   else
       hb(i)=hb(i-1);
        hb_sig(i)=hb_sig(i-1);
       p(i)=p(i-1);
   end
end


figure(1)
plot(1:N,hb);
    yline(10,'--','color','red','LineWidth',1);
    title('Trace plot for the hyper-parameter (mean) under the MCMC method') 
    
figure(2)
histogram(hb);
    xline(median(hb),'--','color','blue','LineWidth',1);
    xline(10,'--','color','red','LineWidth',1);
 title('Posterior of the hyper-parameter (mean) under the MCMC method')

 
save('hb.mat','hb');
save('hb_sig.mat','hb_sig');