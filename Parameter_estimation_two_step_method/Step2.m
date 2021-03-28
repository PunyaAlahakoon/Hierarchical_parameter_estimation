%Step 2
%use this to estimate sub-population specific parameters for an SIRS model
%given hyper-parameters 
%functions need to run this:
        %abc_hie2
        
% conststruct a parellel version for several populations 
tic 
%This is for an SIRS model  for two parameters  
%Number of particles/ parameters sets to run 
B=5000; 
data=load('sc1_data.mat', 'data');
data=data.data;
%data=data(1:30,:);
y_all=data;

dim=size(data,2); %number of pos
T=size(data,1);

pd_b = @(f,h) makedist('Normal','mu',f,'sigma',h);
tpd=@(f,h) truncate(pd_b(f,h),1,10);

    %load hyper-parameters 
    hb=load('hb.mat', 'hb');
    hb=hb.hb;
    hb=hb(5001:10000); 
    hb_sig=load('hb_sig.mat', 'hb_sig');
    hb_sig=hb_sig.hb_sig;
    hb_sig=hb_sig(5001:10000);
 

%m=1;
%the model and initial conditions: SIRS   
%initial conditions 
s0 =999;
i0 = 1;
r0=0;
ini_state=[s0 i0 r0]; %initial population sizes in each compartment

stoi= [-1 1 0;0 -1 1;1 0 -1]; %stoichimetry matrix 
time = 0; %start time to consider 
stp1= @(n) n(2)==0; %stopping criteria  


%Known parameters for all the clusters (assume they are equal in all clusters)
epsilon=0.06;
gamma=1;


eta=1; %number of samples to generate for each parameter set 

%store the final smc-proposed parameters and other needed values in a matrix 
beta_smc=zeros(B,dim);
w_smc=zeros(B,dim); %weights 
E_smc=zeros(G+1,dim);%tolerance values 
AG_smc=zeros(G,dim);%number of steps generated to get B parameters 
s_x= zeros(B,dim);%save the distance criteria

for k=1:dim
    %consider the kth population  
    y=y_all(:,k); %convert the array to a vector
    %stopping criteria, 
    T=length(y);
    stp2=T;
    t_seq=1:T;
    %Store the sets of tolerance values 
    E=[150 70];
    e=E(1);
    %E(1)=e;
    %set a counter for number of iterations to run for each gen 
    AG=zeros(1,G);

      
    %store values for the current generation
    betas0=zeros(1,B); %store the cluster based beta values from posterior 
    w0=zeros(1,B); %store weights 
    ag=0;%set the counter 
    ag0=zeros(1,B);%set the counter 
    rho_m=zeros(1,B);%store the distance values 

     parfor a=1:B %particle number     
     [betas0(a),w0(a),rho_m(a),ag0(a)]=abc_hie2(hb(a),hb_sig(a),y,e,ini_state,stoi,time,stp1,stp2,eta,gamma,epsilon);
     end 
        s_xx= rho_m';%save the distance criteria
        betas=betas0; %store the cluster based beta values from posterior 
        AG(g)=sum(ag0);
        
%store the final values of the population k
beta_smc(:,k)=betas;
w_smc(:,k)=w; %weights 
E_smc(:,k)=E;%tolerance values 
AG_smc(:,k)=AG;%number of steps generated to get B parameters 
s_x(:,k)=s_xx;
end
save('s_x_smc.mat','s_x');
save('beta_smc.mat','beta_smc');
save('E_smc','E_smc');
save('AG_smc.mat','AG_smc');
save('w_smc.mat','w_smc');

toc


