% Function: Generate samples for the optimization
function genSamples(n_samples,E_mean,E_std,rho_mean,rho_std)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   n_samples - number of samples
%   E_mean - Mean value of the Young's modulus
%   E_std - Standard deviation value of the Young's modulus
%   rho_mean - Mean value of the density
%   rho_std - Standard deviation value of the density

% Outputs:
%
%   U_p - Robust objective function (without costs)
%   U_mean - Expected value of the performance
%   U_var - Variance of the performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Val_lim=[E_mean E_std;
         rho_mean rho_std];

n_th=size(Val_lim,1); % Number of parameters 

% Prior PDF definition of model parameters:
% (1) - Young's modulus
% (2) - Density
th=zeros(n_samples,n_th);
% Sampling the prior space for each parameter - NORMAL PDFs
for i=1:n_th
   th(:,i)=normrnd(Val_lim(i,1),Val_lim(i,2),[n_samples,1]);
end

save('./dat/th.mat','th')