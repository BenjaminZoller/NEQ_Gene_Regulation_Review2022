function [M,X,V,E] = makeRateMatrixCooperativity(kb,ku,co)
%MAKERATEMATRIXCOOPERATIVITY makes the state rate matrix (the Laplacian
%matrix of the master equation) for the cooperativity model in Fig. 2 right.
%   [M,X,V,E] = makeRateMatrixCooperativity(kb,ku,co) returns the state
%   rate matrix of the reaction network M, a matrix X that encodes for each
%   state which TF binding sites are occupied, a vector V that determines
%   how many TF are bound for each state and a matrix E whose size is
%   identical to M, which labels the different reactions (the edges of the
%   network). The function takes as input the TF binding rate kn, the TF
%   unbinding rate ku, and a cooperativity vector whose length is equal to
%   the number of TF binding sites.

%given the TF binding rate kb, the TF unbinding
%   rate ku, and the proof-reading rate kq.

%co: cooperativity vector
N = length(co);

m = 2^N; % number of states
M = zeros(m,m);

X = zeros(m,N);
for i=1:N
    x = [ones(2^(N-i),1);zeros(2^(N-i),1)];
    X(:,i) = repmat(x,2^(i-1),1);
end
V = sum(X,2); % #bound TFs
E = zeros(m,m);

for i=1:m
    xi = X(i,:);
    d1 = X - repmat(xi,m,1);
    v1 = sum(abs(d1),2); 
    r1 = find(v1==1);
    
    % TF binding and unbinding
    for j=r1'
        k = find(d1(j,:));
        if d1(j,k) == -1
            %unbind
            M(j,i) = ku/(co(k)^(V(i)-1));
            E(j,i) = -k;
        else
            %bind
            M(j,i) = kb;
            E(j,i) = 1;
        end
    end  
end

Z = sum(M,1);
M = M - diag(Z);

end

