%inistate = initial population sizes in each compartment 
%time = start time to consider 
%stoi= stoichimetry matrix 
%R = reactions
%stp1= stopping criteria 1 (absorption state)
%stp2= time criteria  time to end the algorithm (T<=tx)
%samples = number of Gillespie paths to create 

function [Times,paths]= Gillespe4(ini_state,time,stoi,R,stp1,stp2,samples)
n = size(stoi,1); % Number of transitions
%m = size(stoi,2); % Number of states
B=samples;
paths={1,B}; %store states of all the Gillespie paths
Times={1,B};%store times of all the Gillespie paths  

parfor j=1:B
states= []; 
T=[];%store times
%initilization 
  states= [states;ini_state]; 
  T=[T time];
 
     i=1;
     while ((T<=stp2) &(~stp1(states(i,:))))
      
        %calculate rates 
        current_R = zeros(1,n);
        for k = 1:n
        func = R{k};
        current_R(k) = func(states(i,:));
        end
    
    % Total rate
    rtot = sum(current_R);
    
    %create random nums
    r=rand(1,2);
    
    %create time of the next event 
    tau=(1/rtot)*log(1/r(1));
    Tnext=T(i)+tau;
    T=[T Tnext];
    %choose transition 
    g=rtot*r(2);
    trans=find(cumsum(current_R)>=g, 1 );
    e=min(trans);
    
     state_next= states(i,:)+ stoi(e,:); 
     states=[states;state_next];
     i=i+1;
     
     end 
     paths{1,j}=states;
     Times{1,j}=T;
end
     
         

 
 

