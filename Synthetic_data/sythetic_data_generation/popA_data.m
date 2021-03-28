%generate sub-populations with identical model parameters

% nededed functions:
%Gillespie4 and sample_gen

beta=2;
%mu value for extinction prob=0.5
mu=0.06; %waning immunity rate 
gamma=1; %rate of recovery 
%the model:
k=15; 
%%%%%%population parameters%%%%%
%setting parametrs for SIRS

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
       Ri = {@(n) beta*n(1)*n(2)/(sum(n(1)+n(2)+n(3))-1);...
            @(n) gamma*n(2);@(n) (mu)*n(3)}; %reactions for par1
    %X=zeros(stp2,m);% store number of infecteds (observed)
    data(:,j)=sample_gen(t_seq,ini_state,time,stoi,Ri,stp1,stp2,m);
    plot(t_seq,data(:,j));
    %plot(T{i},Gi{i})
    ylim([0 300]);
    hold on 
end 

%save('sc1_data.mat','data');
