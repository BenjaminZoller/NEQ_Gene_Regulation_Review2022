function [S,T] = genStoTrajectories(M,Ip,kr,dt,tend,k0)
%GENSTOTRAJECTORIES generates stochastic realisations of the model M 
%according to the Gillespie algorithm.
%   [S,T] = genStoTrajectories(M,Ip,kr,dt,tend,k0) returns a single
%   stochastic trajectory as a matrix S, where the first column corresponds
%   to the occupied state, the second column to the occurence of production
%   events, and the third column to the winding number. The function also
%   returns a time vector T corresponding to each row of S. The
%   function takes as input the state rate matrix of the reaction 
%   network M (the Laplacian matrix), a logical vector Ip defining which
%   states are considered active, a production rate kr for generating
%   molecules in the active states, a sampling interval dt over which to
%   represent the trajectory, a total duration tend for the trajectory, and
%   an optional initial state k0 (if not provided, k0 will be drawn from
%   the steady state probability).

% Compute the number of reactions Nr involved in the network
% Also add one reaction per active state for molecule production
Im = M > 0;
Nr = sum(Im(:))+sum(Ip);

% Determine the number of states needed to describe the temporal evolution
% of the network. There is one additional state for the production events
Nxm = size(M,1);
Nx = Nxm+1;

% Build the state change vector V, the reaction rate vector R and a state
% change vector D specific to the winding number
V = zeros(Nr,Nx);
R = zeros(Nr,1);
D = zeros(Nr,1);
ii = 0;
for i=1:Nxm
    Im = M(i,:) > 0;
    for j=find(Im)
        ii = ii+1;
        V(ii,i) = 1;
        V(ii,j) = -1;
        R(ii) = M(i,j);
        if j>i
            if i==1 && j==Nxm
                D(ii) = 1;
            else
                D(ii) = -1;
            end
        else
            if i==Nxm && j==1
                D(ii) = -1;
            else
                D(ii) = 1;
            end
        end
    end
end
V((ii+1):end,end) = 1;
R((ii+1):end) = kr;
V0 = V<0;
j = find(Ip);
for i=1:sum(Ip)
    V0(ii+i,j(i)) = 1;
end

% In case the initial condition is not provided,
% initialize state a t0 according to the steady state probability
% X is the state vector at time t
X = zeros(1,Nx);
r = rand(1);
if nargin < 6
    % compute steady state
    [~,P] = getExp(M,Ip);
    A = cumsum(P);
    
    k = 1;
    while r > A(k)
        k = k+1;
    end
    X(k) = 1;
else
    X(k0) = 1;
end

% Sampling the stochastic trajectories according to the Gillespie algorithm
T = 0:dt:tend;
Nt = numel(T);
S = zeros(Nt,3);
W = 0;
j = 1;
tc = 0;
while tc < tend
    Ix = logical(X);
    Ix(end) = false;
    Ir = V0(:,Ix);
    P = Ir.*R;    
    A = cumsum(P);
    r = rand(2,1);
    t = -log(r(1))/A(end);
    A = A/A(end);
    
    k = 1;
    while r(2) > A(k)
        k = k+1;
    end
    
    tc = tc + t;
    
    while tc > T(j)
        S(j,1) = find(X(1:(end-1)));
        S(j,2) = X(end);
        S(j,3) = W;
        X(end) = 0;
        j = j + 1;
        if j > Nt
            break
        end
    end
    
    X = X + V(k,:);
    W = W + D(k);
end  

end



