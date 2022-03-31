function [S,T] = genRastTrajectories(M,Ip,kr,tau,dt,tend,k0)
%GENRASTTRAJECTORIES generates stochastic realisations of the model M
%following a raster approach based on the propagator of the master equation.
%   [S,T] = genRastTrajectories(M,Ip,kr,tau,dt,tend,k0) returns a single
%   stochastic trajectory as a matrix S, where the first column corresponds
%   to the occupied state, the second column to the occurence of any active
%   states, and the third column to the accumulation produced/expressed
%   molecules. The function also returns a time vector T corresponding to
%   each row of S. The function takes as input the state rate matrix of
%   the reaction network M (the Laplacian matrix), a logical vector Ip
%   defining which states are considered active, a production rate kr for
%   generating molecules in the active states, a mean liftime for the
%   molecules, a sampling interval dt over which to represent the trajectory,
%   a total duration tend for the trajectory, and an optional initial
%   state k0 (if not provided, k0 will be drawn from the steady state
%   probability).

% Determine the number of states needed to describe the temporal evolution
% of the network.
Nx = size(M,1);

% Compute the propagator PP of the master 
PP = expm(M*dt);
% PP columns are supposed to be normalized to 1
% But for numerical reasons the normalization can be slightly off
% Here we force the normalization to 1 to avoid sampling issues below
PP = PP .* repmat(1./sum(PP),Nx,1);
V = double(Ip);

% In case the initial condition is not provided,
% initialize state a t0 according to the steady state probability
% X is the state vector at time t
X = zeros(Nx,1);
r = rand(1);
if nargin < 7
    % compute steady state
    [~,P] = getExp(M,Ip);
    C = cumsum(P);
    
    k = 1;
    while r > C(k)
        k = k+1;
    end
    X(k) = 1;
else
    X(k0) = 1;
end

% Sampling the stochastic trajectories according to a raster approach
% based on the propagator of the master equation. This approach can be much
% faster when dealing with models whose underlying reaction rates are
% vastly different (several order of magnitudes). However, depending on dt,
% the raster approach will skip over fast occuring reactions.
T = 0:dt:tend;
Nt = numel(T);
S = zeros(Nt,3);
S(1,1) = k;
S(1,2) = V(k);
for i=2:Nt
    P = PP*X;
    C = cumsum(P);
    r = rand(1);
    k = 1;
    while C(k) < r
        k = k+1;
    end
    S(i,1) = k;
    S(i,2) = V(k);
    X = zeros(Nx,1);
    X(k) = 1;
end

% Let's do convolution
tt = 0:dt:5*tau;
kernel = kr*dt*exp(-tt/tau);
P = conv(S(:,2),kernel');
S(:,3) = P(1:length(T));

end
