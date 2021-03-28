%use this to create samples from Gillespie to read at each time point you
%want 
%input:
%time as a sequence 
%all the other inputs for the Gillespie4 algorithm 
%output:
%xs= gives the numbers in the second compartment. 

function[xs]=sample_gen(t_seq,ini_state,time,stoi,R,stp1,stp2,samples)

[Times,paths]=Gillespe4(ini_state,time,stoi,R,stp1,stp2,samples);
m=length(t_seq);
xs=zeros(m,samples);

parfor i=1:samples
e=paths{1,i}(:,2);
T=Times{1,i};

statesi=zeros(1,m);
for j=1:m
    l=max(find(T<=t_seq(j))); %index of T that is closest to tj
    statesi(j)=e(l);
end

xs(:,i)=statesi;
end 
end