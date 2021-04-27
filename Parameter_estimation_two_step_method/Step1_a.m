%Step (1)--a
%THIS IS THE CODE FOR STEP (1)--a IN THE PAPER--- YOU CAN REPLACE THIS STEP
%WITH OTHER SOFTWARE PACKAGES 
%use this to estimate betas of an SIRS model  using ABC-SMC by Toni et al 
% conststruct a parellel version for several populations 
%functions need to run this:
    %Gillespie4-- The Gillespie algorithm
    %sample_gen --- reads data at discrete times
    %abc_ind --- ABC step--- MAKE SURE TO CHANGE THE PRIOR AS REQUIRED! 
    
tic 
%load data
data=load('sc1_data.mat', 'data'); % load the dataset you need
y_all=data.data; % you can insert data with different number of observed points
                    %in each sub-population

dim=size(y_all,2);%number of sub-populations

%the model and initial conditions: SIRS    
%initial conditions 
s0 =999; % #of  susceptibles 
i0 = 1; % # of infectious 
r0=0;% # of recovered 
ini_state=[s0 i0 r0]; %initial sizes in each compartment

stoi= [-1 1 0;0 -1 1;1 0 -1]; %stoichimetry matrix 
time = 0; %start time to consider 
stp1= @(n) n(2)==0; % 1st stopping criteria for Gillespie (i.e., stop the algorithm when the #of infectious=0)
%Known parameters for all the sub-populations:
epsilon=0.06; % waning immunity rate 
gamma=1; % rate of recovery 

%Number of particles/ parameters sets to sample using ABC-SMC
B=5000; 

%ABC-SMC initilisation 
%number of generations to run 
G=7; 
eta=1; %number of samples to generate for each parameter set 

%store the final smc-sampled parameters and other needed values in a matrix 
beta_smc=zeros(B,dim);
w_smc=zeros(B,dim); % store weights
E_smc=zeros(G+1,dim);%tolerance values-- not neccssary unless you decide to
                       %use dynamic methods to calculate the tolarance 
AG_smc=zeros(G,dim);%store the # of particles generated to get B parameters(if needed) 
s_x= zeros(B,dim);%save the distance criteria (if needed)

for k=1:dim
    %consider the kth population  
    y=y_all(:,k); %convert the array to a vector
    %stopping criteria, 
    T=length(y);
    stp2=T; % 2nd stopping criteria for Gillespie (i.e., the time-period 
            %the sample path is generated if the first criteria is not met
            %during this perid)
    t_seq=1:T;
    %Store the sets of tolerance values 
    E=[350 300 270 250 200 170 150 90]; %these are the values used for Pop A
    %starting tolerance levels for each pop: starting values 
    e=E(1);
    %set a counter for number of iterations to run for each gen 
    AG=zeros(1,G);

    %store parameter sets
    betas=zeros(1,B); %store betas
    w=zeros(1,B); %store weights 
    
    g=1;
    
while (g<1+G) %number of generation 
    %store values for the current generation
    betas0=zeros(1,B); %store the sub-pop based beta values from posterior 
    w0=zeros(1,B); %store weights 
    ag=0;%set the counter 
    ag0=zeros(1,B);%set the counter 
    rho_m=zeros(1,B);%store the distance values 

     parfor a=1:B %particle number     
     [betas0(a),w0(a),rho_m(a),ag0(a)]=abc_ind(g,B,betas,w,y,e,ini_state,stoi,time,stp1,stp2,eta,gamma,epsilon);
     end 
        s_xx= rho_m';%save the distance criteria
        e=E(g+1);
        %normalize the weights 
        w0=normalize(w0,'norm',1);
        %replecae the previous gen values from the new gen values
        betas=betas0; %store the cluster based beta values from posterior 
        w=w0; %store weights 
        AG(g)=sum(ag0);
        g=g+1; %update the generation number 
end
%store the final values of the population k
beta_smc(:,k)=betas;
w_smc(:,k)=w; %weights 
E_smc(:,k)=E;%tolerance values 
AG_smc(:,k)=AG;%number of steps generated to get B parameters 
s_x(:,k)=s_xx;
end
save('s_x_smc.mat','s_x');
save('ind_beta_smc.mat','beta_smc');
save('E_smc','E_smc');
save('AG_smc.mat','AG_smc');
save('w_smc.mat','w_smc');
toc


