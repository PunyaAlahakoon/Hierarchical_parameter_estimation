%use this algorithm to estimate HPD (highest posterior density) for any
%sub-pops
% intervals using the Chen-Shao algorithm
%Inputs: MCMC sample of the posteriors of the model 

%Input hyper-parameters as a matrix (hy)
%Input pop-specific parameters as an array 


function[hbd]=HPD_estimation2(betas)
    hbd={1,2};

n=size(betas,1); %number of sampled values of sub-pop-specifc params 
k=size(betas,2); %number of pops 

    %sort all the matrices 

    s_betas=sort(betas);
 
    %2)for pops
    hpd_betas=zeros(k,2);
    nc=(0.95*n);
    n_ci=n- nc; %number of credible intervals to calc
    for j=1:k
    ci_betaj=zeros(n_ci,2); %store the credible intervals 
    dist_beta=zeros(1,n_ci); %calculate the distance between the cis 
    for i=1:n_ci
    ci_betaj(i,:)=[s_betas(i,j) s_betas(i+nc,j)];
    dist_beta(i)=abs(s_betas(i+nc,j)- s_betas(i,j));
    end
%find the credible interval with the shortest widths
ind_b=find(dist_beta==min(dist_beta));
if length(ind_b)>1
    ind_b=ind_b(1);
end
hpd_betas(j,:)=ci_betaj(ind_b,:);
end
hbd{1,1}=hpd_betas;
end