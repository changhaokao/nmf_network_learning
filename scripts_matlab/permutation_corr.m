function [rho_data, pval_perm] = permutation_corr(xData,yData,nPerm)

% number of permutation
if nargin<3
    nPerm = 10000;
end

nSub = numel(xData);

rho_perm = NaN(nPerm,1);
for idx_p = 1:nPerm
    
    x_perm = xData(randperm(nSub));
    [rho] = corr(x_perm,yData);
    rho_perm(idx_p,1) = rho;
    
end
rho_data = corr(xData,yData);
if rho_data>=0
    pval_perm = mean(rho_perm>rho_data);
else
    pval_perm = mean(rho_perm<rho_data);
end
