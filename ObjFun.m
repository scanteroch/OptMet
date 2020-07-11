% Objective function to assess the design of a U beam with resonators
function [U_p, U_mean, U_var]=ObjFun(m, n, th, A, fmin, fmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   m  - The percentage mass
%   n  - Number of resonators
%   th - Samples acquired from the prior p(\th)
%   A  - Trade-off variable between expectation and variance in [0,1]
%   fmin - Minimum frequency in the design interval for vibration
%          attenuation
%   fmax - Maximum frequency in the design interval for vibration
%          attenuation
%
% Outputs:
%
%   U_p - Robust objective function (without costs)
%   U_mean - Expected value of the performance
%   U_var - Variance of the performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Monte Carlo approximation of the expectation over \theta
	n_samples=length(th(:,end));
    U_mean_aux=zeros(n_samples,1);
    for i=1:n_samples
        U_mean_aux(i)=FRFresp(n,m,th(i,:),fmin,fmax);
    end
    U_mean=1/n_samples*sum(U_mean_aux);
    U_var=1/n_samples*sum((U_mean_aux-U_mean).^2);
    
    % Objective function (without costs):
    U_p=(1-A)*U_mean+A*U_var; % U_mean<0, U_var>0

end