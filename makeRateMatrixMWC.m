function [M,configBound,configTF,configLink] = ...
    makeRateMatrixMWC(k_bind,k_unbind,k_bindM,k_unbindM,alpha,n)
%MAKERATEMATRIXMWC makes the state rate matrix (the Laplacian
%matrix of the master equation) for the equilibrium MWC model in Fig. 3.
%   [M,configBound,configTF,configLink] = makeRateMatrixNeqMWC(k_bind,
%   k_unbind,k_bindM,k_unbindM,alpha,n) returns the state
%   rate matrix of the reaction network M and a matrix configBound that
%   encodes for each state the number of bound TF, the number of links and
%   whether Mediator is bound or not. The function also returns two
%   matrices configTF and configLink that determine for each state which TF
%   binding sites is occupied and which TF is linked to Mediator.
%   The function takes as input arguments the different rates of the model,
%   namely the TF binding rate k_bind, the TF unbinding rate k_unbind,
%   the Mediator binding rate k_bindM, the Mediator unbinding rates
%   k_unbindM, the cooperativity parameter alpha and the number of TF
%   binding sites n.

m = 2^(1+n);
M = zeros(m,m);

X = zeros(m,1+n);
for i=1:(1+n)
    x = [ones(2^(n-i+1),1);zeros(2^(n-i+1),1)];
    X(:,i) = repmat(x,2^(i-1),1);
end

Vb = sum(X(:,1+(1:n)),2);
Vl = zeros(m,1);
Vl(1:(m/2)) = Vb(1:(m/2));

for i=1:m
    xi = X(i,:);
    mi = X(i,1);
    d1 = X(:,1:(n+1)) - repmat(xi(1:(n+1)),m,1);
    v1 = sum(abs(d1),2);
    
    rm = find(v1==1&abs(d1(:,1))==1);
    r1 = find(v1==1&d1(:,1)==0);
    
    % mediator binding and unbinding    
    if d1(rm,1) == -1
        M(rm,i) = k_unbindM/alpha^Vl(i);
    else
        M(rm,i) = k_bindM;
    end
        
    % TF binding and unbinding
    for j=r1'
        k = find(d1(j,:));
        if d1(j,k) == -1         
            M(j,i) = k_unbind/alpha^mi;
        else
            M(j,i) = k_bind;
        end
    end  
end

Z = sum(M,1);
M = M - diag(Z);

configBound = [Vb,Vl,X(:,1)];
configTF = X(:,1+(1:n));
configLink = [X(1:(m/2),1+(1:n));zeros(m/2,n)];

end