%ABC step ====DO NOT CHANGE ANYTHING
%INPUT PARAMETERS:
%g = current generation number 
%B = # of particles to obtain from posterior 
%betas = samples particles in the previous generation 
%w= weights in the previous generation 
%y= data
%e = tolerance value 

function [betas0,w0,rho0,ag0]= abc_ind(g,B,betas,w,y,e,ini_state,stoi,time,stp1,stp2,eta,gamma,epsilon)
%model definition 
aa=0; %to make sure that you get an accepted particle from the output 
ag=0; %set the counter
%initilise
betas0=0;w0=0;
rho0=0;
while(aa<2)
    ag=ag+1;
 if g==1 
        be=unifrnd(1,10); % sample from prior    
    else    

        %find the index to the parameter set to use 
        ind=randsample(1:B,1,true,w);
        be0=betas(ind);
        %perturbate:
        ss=std(betas);
        be=abs(normrnd(be0,ss));
   
 end  
         p1=unifpdf(be,1,10);
         
    if p1>0
                
                t_seq=1:stp2;
                Rj = {@(n) be*n(1)*n(2)/(sum(n(1)+n(2)+n(3))-1);...
                @(n) gamma*n(2);@(n) (epsilon)*n(3)};
                %create eta sample sets using Rj:
                X=sample_gen(t_seq,ini_state,time,stoi,Rj,stp1,stp2,eta);
                %distance metric
                l=abs(X-y);
                cx=sqrt(sum((l.^2)));
                S=cx;
            if S<=e
                betas0=be;
                %store s
                rho0=S;
                %store the weight w
                if g==1
                    w0=1;
                else
                %denominator:
                     pd=normpdf(be,betas,ss);
                     den=w.*pd;
                    w0=p1/sum(den);   
                end
                %update a=a+1
                aa=aa+1;
                ag0=ag;
            end   
    end
end
end

