%DO NOT CHANGE ANY NOTATION. THIS FUNCTION IS USED IN STEP 2

function [betas0,rho0,ag0]= abc_hie2(hb,hb_sig,y,e,ini_state,stoi,time,stp1,stp2,eta,gamma,epsilon)
%model definition 
aa=0; %to make sure that you get an accepted particle from the output 
ag=0; %set the counter
pd_b = @(f,h) makedist('Normal','mu',f,'sigma',h);
tpd=@(f,h) truncate(pd_b(f,h),0.001,10);
tpdb=tpd(hb,hb_sig);
 %p1=unifpdf(hb,1,10)*unifpdf(hb_sig,0,2.5);
%initilise
betas0=0;
rho0=0;
while(aa<2)
    ag=ag+1;
 
        be=random(tpdb);
        p2=pdf(tpdb,be);
         
    if p2>0
                                
                t_seq=1:stp2;
                Rj = {@(n) be*n(1)*n(2)/(sum(n(1)+n(2)+n(3))-1);...
                @(n) gamma*n(2);@(n) (epsilon)*n(3)};
                %create eta sample sets using Rj
                X=sample_gen(t_seq,ini_state,time,stoi,Rj,stp1,stp2,eta);
                %distance metric
                l=abs(X-y);
                cx=sqrt(sum((l.^2)));
                S=cx;
            if S<=e
                betas0=be;
                %store s
                rho0=S;
                
                %update a=a+1
                aa=aa+1;
                ag0=ag;
            end   
    end
end
end

