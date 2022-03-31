function [M,configBound,configTF,configLink,Ed] = ...
    makeRateMatrixNeqMWC(k_bind,k_unbind,k_bindM,k_unbindM,k_link,k_unlink,alpha,n)
%MAKERATEMATRIXNEQMWC makes the state rate matrix (the Laplacian
%matrix of the master equation) for the non-equilibrium MWC model in Fig. 3.
%   [M,configBound,configTF,configLink,Ed] = makeRateMatrixNeqMWC(k_bind,
%   k_unbind,k_bindM,k_unbindM,k_link,k_unlink,alpha,n) returns the state
%   rate matrix of the reaction network M and a matrix configBound that
%   encodes for each state the number of bound TF, the number of links and
%   whether Mediator is bound or not. The function also returns two
%   matrices configTF and configLink that determine for each state which TF
%   binding sites is occupied and which TF is linked to Mediator. Lastly,
%   the function also returns a matrix Ed whose size is identical to M,
%   which labels the different reactions (the edges of the network).
%   The function takes as input arguments the different rates of the model,
%   namely the TF binding rate k_bind, the TF unbinding rate k_unbind,
%   the Mediator binding rate k_bindM, the Mediator unbinding rates
%   k_unbindM, the linking rate k_link, the unlinking rate k_unlink, the
%   cooperativity parameter alpha and the number of TF binding sites n.
%
%   Copyright (c) 2022, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree.

% number of states
m0 = 2^n;
m1 = 0;
for i=0:n
   m1 = m1 + nchoosek(n,i)*2^i;
end
m = m0+m1;
M = zeros(m,m);

% list of states
% [sm, s1, s2, ..., sn, m1, m2, ..., mn]

X0 = zeros(m0,1+2*n);
for i=1:n
    x = [ones(2^(n-i),1);zeros(2^(n-i),1)];
    X0(:,i+1) = repmat(x,2^(i-1),1);
end

X1 = zeros(m1,1+2*n);
X1(:,1) = 1;
ii = 0;
for i=1:m0
    v = X0(i,2:(n+1));
    k = find(v);
    s = sum(v);
    l = 2^s;
    w = zeros(l,n);
    for j=1:s
        x = [ones(2^(s-j),1);zeros(2^(s-j),1)];
        w(:,k(j)) = repmat(x,2^(j-1),1);
    end
    X1(ii+(1:l),2:end) = [repmat(v,l,1),w];
    ii = ii + l;
end

X = [X1;X0];

Vb = sum(X(:,1+(1:n)),2);
Vl = sum(X(:,n+1+(1:n)),2);

D = X(:,2:(n+1))+X(:,n+1+(1:n));
Ed = zeros(size(M));

for i=1:m
    xi = X(i,:);
    d1 = X(:,1:(n+1)) - repmat(xi(1:(n+1)),m,1);
    d2 = X(:,n+1+(1:n)) - repmat(xi(n+1+(1:n)),m,1);
    d3 = D - repmat(D(i,:),m,1);
    d3(any(d3>1,2),:) = 0;
    d3(d3~=0) = 1;
    v1 = sum(abs(d1),2);
    v2 = sum(abs(d2),2);
    v3 = sum(d3,2);
    
    rm = find(v1==1&abs(d1(:,1))==1&~any(d2>0,2));
    r1 = find(v1==1&d1(:,1)==0&v3==1);
    r2 = find(v1==0&v2==1);
    
    % mediator binding and unbinding    
    if d1(rm,1) == -1
        M(rm,i) = k_unbindM/alpha^Vl(i);
        Ed(rm,i) = -0.5;
    else
        M(rm,i) = k_bindM;
        Ed(rm,i) = 0.5;
    end
        
    % TF binding and unbinding
    for j=r1'
        k = find(d1(j,:));
        if d1(j,k) == -1         
            M(j,i) = k_unbind/alpha^xi(n+k);
            Ed(j,i) = -1;
        else
            M(j,i) = k_bind;
            Ed(j,i) = 1;
        end
    end
    
    % links binding and unbinding
    for j=r2'
        k = find(d2(j,:));
        if d2(j,k) == -1 
            M(j,i) = k_unlink;
            Ed(j,i) = -1.5;
        else
            M(j,i) = k_link;
            Ed(j,i) = 1.5;
        end
    end    
end

Z = sum(M,1);
M = M - diag(Z);

% configBound = [number of TF bound, number of links, presence of Mediator]
configBound = [Vb,Vl,X(:,1)]; 
configTF = X(:,1+(1:n));
configLink = X(:,n+1+(1:n));

end

