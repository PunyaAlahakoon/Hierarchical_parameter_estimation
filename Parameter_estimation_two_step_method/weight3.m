%DO NOT CHANGE ANY NOTATION. THIS FN IS USED IN STEP1_B
%me = samled betas from step 1_a
%b= proposed psi_beta,
%hbs= proposed sigma_beta,
%dim= # of sub-pops,
%B = # of samled betas from step 1_a

function[w,poste]=weight3(me,b,hbs,dim,B)
    %distributions of theta|psi:
     w=[];
    %poste=zeros(#OF PARTICLES,# OF SUB-POPS);
    poste=zeros(B,dim);
parfor k=1:dim
    w0=@(x,n) normpdf(me(:,n)',x,hbs);
    wj=w0(b,k);  
    poste(:,k)=wj;
     w=[w sum(wj(wj>0))];
end
w=prod(w);
end