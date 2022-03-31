function [Phi] = getNoise(M,Ip,tp)
%GETNoise computes the steady state noise propagation of the model M.
%   [Phi] = getNoise(M,Ip,tp) returns the steady state value of the 
%   filtering function Phi, which is normalized from 0 to 1.
%   The function takes as input the state rate matrix of the reaction 
%   network M (the Laplacian matrix), a logical vector Ip defining which
%   states are considered active (leads to expression) and a filtering
%   time scale tp, which typically corresponds to the mean lifetime of
%   the expressed molecules (mRNAs or proteins).
%
%   Copyright (c) 2022, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree.

% The derivation of the propagated noise Phi can be found in 
% Grah et al. 2020, DOI:10.1073/pnas.2006731117 

% Reduced matrix (M is degenerate)
Mr = M;
N = size(Mr,1)-1;
k0 = Mr(1:N,end);
Mr = Mr(1:N,1:N) - k0*ones(1,N);

% Compute steady state
P = -Mr\k0;

% Build v0 & v1 vector for the noise filtering function involved in
% noise propagation
Np = sum(Ip);
mn = sum(P(Ip));
sn = mn*(1-mn);

v0 = zeros(1,N);
v0(Ip) = 1;

v1 = zeros(N,1);
ii = find(Ip);
for i=1:Np
    mi = P(ii(i));
    vi = -mi*P;
    vi(ii(i)) = mi*(1-mi);
    v1 = v1 + vi;
end

% Compute the filtering function (noise propagation)
Phi = v0*((eye(N)-tp*Mr)\v1)/sn;

end

