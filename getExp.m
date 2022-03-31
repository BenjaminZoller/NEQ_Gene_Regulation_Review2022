function [E,P] = getExp(M,Ip)
%GETEXP computes the steady state expression level E and the occupancies of
%the model M.
%   [E,P] = getExp(M,Ip) returns the steady state expression level E 
%   normalized from 0 to 1, and the steady state occupancies P of the
%   model. The function takes as input the state rate matrix of the reaction 
%   network M (the Laplacian matrix) and a logical vector Ip defining which
%   states are considered active
%
%   Copyright (c) 2022, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree.

% The steady state probabilities are obtained by solving the master
% equation with lhs (time derivative) set to zero.
% Further explanation can be found in 
% Grah et al. 2020, DOI:10.1073/pnas.2006731117 

% Reduced matrix (M is degenerate)
Mr = M;
N = size(Mr,1)-1;
k0 = Mr(1:N,end);
Mr = Mr(1:N,1:N) - k0*ones(1,N);

% Compute the steady state occupancies
P = -Mr\k0;
P = [P;1-sum(P)];

% Compute the expression level
E = sum(P(Ip));

end

