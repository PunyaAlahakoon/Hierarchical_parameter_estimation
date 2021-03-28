%Scenario 2: sub-populations with variable model parameters: 
% nededed functions:
%Gillespie4 and sample_gen

% generate random betas with a hierarchy 
pd = makedist('Normal','mu',2,'sigma',0.5);
tpd=truncate(pd,1,10);

k=15;
t_betas=random(tpd,1,k);
histogram(t_betas);
save('true_sc2_betas.mat','t_betas');

%use this to generate synthetic data from a hierachical structure
k=15; %number of clusters 
%load h_betas of hierachical models: h_beta=2 
betas=load('true_sc2_betas.mat', 't_betas');
betas=betas.t_betas;
%mu value for extinction prob=0.5
mu=0.06;
gamma=1;
%the model:

%%%%%%population parameters%%%%%
%setting parametrs for SIR with demography

s0 =999;
i0 = 1;
r0=0;
ini_state=[s0 i0 r0]; %initial population sizes in each compartment

stoi= [-1 1 0;0 -1 1;1 0 -1]; %stoichimetry matrix 
time = 0; %start time to consider 
stp1= @(n) n(2)==0; %stopping criteria a 
stp2=30;
t_seq=1:stp2;

m=1; %number of gillespie paths created in Gillespie 4 (keep this at 1 for this algorithm)
data=zeros(stp2,k); %store data of observed value

for j=1:k
   subplot(3,5,j);
    Ri = {@(n) betas(j)*n(1)*n(2)/(sum(n(1)+n(2)+n(3))-1);...
            @(n) gamma*n(2);@(n) (mu)*n(3)}; %reactions for par1
    %X=zeros(stp2,m);% store number of infecteds (observed)
    %[Times,paths]=Gillespe4(ini_state,time,stoi,Ri,stp1,stp2,m);

    data(:,j)=sample_gen(t_seq,ini_state,time,stoi,Ri,stp1,stp2,m);
    plot(t_seq,data(:,j));
    ylim([0 300]);
   
end

%save('sc2_data.mat','data');

