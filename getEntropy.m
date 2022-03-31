function [S,K] = getEntropy(M,P)
%GETENTROPY computes the steady state entropy production S and the currents
%K of the model M.
%   [S,K] = getEntropy(M,P) returns the steady state entropy production S
%   and the currents K of the model, computed according to Schnakenberg
%   1976. The function takes as input the state rate matrix of the reaction 
%   network M (the Laplacian matrix) and the steady state occupancies P of
%   the model.

% The derivation of the entropy production can be found in 
% Schnakenberg 1976, DOI:https://doi.org/10.1103/RevModPhys.48.571

% Size of the system
N = size(M,1);

% Compute the currents at steady state
K = M .* repmat(P',N,1) - M' .* repmat(P,1,N);

% Compute entropy production (medium) at steady state
% the system entropy vanishes at steady state
S = K.*log(M./M');
S = 0.5*nansum(S(:));

end

