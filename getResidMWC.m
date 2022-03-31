function [T] = getResidMWC(alpha,ku,eL,E0,N)
%GETRESIDMWC computes the mean TF residence time analytically for the
%standard (equilibrium) MWC model.
%   [T] = getResidMWC(alpha,ku,eL,E0,N) returns the mean TF residence time
%   T given the coopertivity parameter alpha, the TF unbinding rate ku,
%   the Mediator rate ratio eL = kb_M/ku_M, the targeted expression level
%   E0 and the number of TF binding sites N.

% The analytical solution can be found in 
% Grah et al. 2020, DOI:10.1073/pnas.2006731117

x = (eL*(1/E0-1))^(1/N);
eE = (x-1)./(1-x.*alpha);
W1 = (1+eE).^(N-1);
W2 = eL.*(1+eE.*alpha).^(N-1);
T = (1./ku).*(W1+W2.*alpha)./(W1+W2);
end

