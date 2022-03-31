function [P,Dph,S,K,Resid] = computePheno(M,Ip,W)
%COMPUTEPHENO computes dynamical regulatory phenotypes of the model M.
%   [P,Dph,S,K,Resid] = computePheno(M,Ip,W) returns various regulatory
%   phenotypes, such as P the steady state probability, Dph dynamical
%   phenotypes including noise, power spectrum, auto-correlation, transient
%   relaxation and entropy production during transient. The function also
%   returns the steady state entropy production S and current K, and the
%   residence time properties Resid of the active and inactive states. The
%   function takes as input the state rate matrix of the reaction 
%   network M (the Laplacian matrix), a logical vector Ip defining which
%   states are considered active and a time/frequency vector W over which to
%   compute the dynamical phenotypes.
%
%   Copyright (c) 2022, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

% The steady state probabilities and the propagated noise Phi are obtained
% by solving the master equation with lhs (time derivative) set to zero.
% The residence time pdf is given by continuous phase-type distribution.
% Further explanation can be found in 
% Grah et al. 2020, DOI: 10.1073/pnas.2006731117 

% Regarding the power spectrum and the auto-correlation computation,
% detailed explanations can be found in
% Lestas et al. 2010, DOI: 10.1109/TAC.2007.911347

% The derivation of the entropy production can be found in 
% Schnakenberg 1976, DOI: https://doi.org/10.1103/RevModPhys.48.571

% Size of the system
N = size(M,1);

% W acts as either a time or frequency vector
Nw = length(W);

% Reduced matrix (M is degenerate)
Mr = M;
Nr = size(Mr,1)-1;
k0 = Mr(1:Nr,end);
Mr = Mr(1:Nr,1:Nr) - k0*ones(1,Nr);

% Compute steady state
Pr = -Mr\k0;
P = [Pr;1-sum(Pr)];

% Build v0 & v1 vector for the noise filtering function related to
% noise propagation
Np = sum(Ip);
mn = sum(Pr(Ip));
sn = mn*(1-mn);

v0 = zeros(1,Nr);
v0(Ip) = 1;

v1 = zeros(Nr,1);
ii = find(Ip);
for i=1:Np
    mi = Pr(ii(i));
    vi = -mi*Pr;
    vi(ii(i)) = mi*(1-mi);
    v1 = v1 + vi;
end

% Build steady state covariance V for power spectrum and auto-correlation
% function. V is the multinomial covariance.
% off diagonal elements
V = -P*P';
% diagonal elements
V(1:N+1:end) = P.*(1-P);
I = double(Ip);

% Binomial variance of the overall active state
% can be used to normalized the power spectrum and the auto-correlation
s2 = I'*V*I; 

% Initial condition for transient relaxation to the active state(s)
P0 = P;
P0(Ip) = 0;
P0 = P0/sum(P0);

% Variables to store dynamical phenotypes
Dph = nan(6,Nw);

% Complex variable for power spectrum computation
myc = complex(0,1);

for i=1:Nw
    % Compute the filtering function (noise propagation)
    Dph(1,i) = v0*((eye(Nr)-W(i)*Mr)\v1) /sn;
    
    % Compute the power spectrum
    Si = (myc*W(i)*eye(N)-M)\V + (-myc*eye(N)*W(i)-M)\V;
    Dph(2,i) = I'*Si*I;
    
    % Compute the propagator for auto-correlation and relaxation
    Ei = expm(M*W(i));
    
    % Compute the auto-correlation
    Si = Ei*V;
    Dph(3,i) = I'*Si*I;
    
    % Compute the transient relaxation
    Pi = Ei*P0;
    Dph(4,i) = I'*Pi;
    
    % Compute the entropy production during transient
    Ki = M .* repmat(Pi',N,1) - M' .* repmat(Pi,1,N);
    Ssys = Ki.*log(repmat(Pi',N,1)./repmat(Pi,1,N));
    Dph(5,i) = 0.5*nansum(Ssys(:));
    Smed = Ki.*log(M./M');
    Dph(6,i) = 0.5*nansum(Smed(:));
end

% The power specta density is normalized such that the integral of pws from
% zero to infinity is equal to the binomial variance s2.
Dph(2,:) = Dph(2,:)/pi;

% The auto-correlation is normalized such that value at lag zero is equal
% to the binomial variance s2. Alternatively, on could divide the ac by s2
% to normalized it to 1 at lag zero (line below).
%Dph(3,i) = Dph(3,:)/s2;

% Transient relaxation of active state(s) is normalized such that it
% reaches 1 at steady state
Dph(4,:) = Dph(4,:)/sum(P(Ip));

% Compute the currents at steady state
K = M .* repmat(P',N,1) - M' .* repmat(P,1,N);

% Compute the entropy production (medium) at steady state
% the system entropy vanishes at steady state
S = K.*log(M./M');
S = 0.5*nansum(S(:));

% Compute residence time
% for active states
[Resid.m1Ta,Resid.s1Ta,Resid.tTa,Resid.FTa] = ComputeResid(M,P,Ip,Nw);
% for inactive states
[Resid.m1Ti,Resid.s1Ti,Resid.tTi,Resid.FTi] = ComputeResid(M,P,~Ip,Nw);

end

function [m1,s1,t,pdf] = ComputeResid(M,P,Is,Nw)
% Residence time pdf is given by continuous phase-type distribution.
% Initial condition corresponds to the probability that the system has
% just settle in states Is
a = M(Is,~Is)*P(~Is);
a = a/sum(a);

% Reduced matrix, non Is states become absorbing states
M = M(Is,Is);
I = ones(1,sum(Is));
T = -I/M;
% Compute mean, variance and std of the residence time
m1 = T*a;
s2 = abs(2*(I/M^2)*a - m1^2);
s1 = sqrt(s2);

% Time vector over which to compute the pdf
t = [0,logspace(min(-log10(m1),-3),log10(1e1*m1),Nw-1)];
pdf = nan(1,Nw);

% Compute the pdf as a continuous phase-type distribution
m0 = -I*M;
for i=1:Nw
    expTF = expm(M*t(i));
    pdf(i) = m0*expTF*a;
end
end