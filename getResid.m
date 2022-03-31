function [m1T,s1T] = getResid(M,Ip,P)
%GETRESID computes the mean and std residence time for a set of states
%of the model M.
%   [m1T,s1T] = getResid(M,Ip,P) returns the mean m1T and the standard
%   deviation s1T of the residence time within specified states.
%   The function takes as input the state rate matrix of the reaction 
%   network M (the Laplacian matrix), a logical vector Ip defining for which
%   states the residence time must be computed, and the steady state
%   occupancies P of the model (which can be computed using the getExp
%   function). Providing the steady state occupancies P is optional, but
%   it is highly recommended if there are more than 1 active state (sum(Ip)>1).
%   Otherwise the residence time might be inccorect, as the system will be
%   assumed to have just settled in the first Ip state.
%
%   Copyright (c) 2022, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree.

% Residence time pdf is given by continuous phase-type distribution.
% Further explanation can be found in 
% Grah et al. 2020, DOI:10.1073/pnas.2006731117 

% Initial condition corresponds to the probability that the system has
% just settle in states Ip
if nargin < 3
    N = size(M,1);
    a = zeros(N,1);
    a = a(Ip);
    a(1) = 1;
else
    a = M(Ip,~Ip)*P(~Ip);
    a = a/sum(a);
end

% Reduced matrix, non Ip states become absorbing states
M = M(Ip,Ip);
I = ones(1,sum(Ip));
T = -I/M;
% Compute mean, variance and std of the residence time
m1T = T*a;
s2T = abs(2*(I/M^2)*a - m1T^2);
s1T = sqrt(s2T);

end

