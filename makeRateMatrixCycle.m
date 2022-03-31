function [M] = makeRateMatrixCycle(wf,wb)
%MAKERATEMATRIXCYCLE makes the state rate matrix (the Laplacian matrix of
%the master equation) for the cycle model in Fig. 1.
%   [M] = makeRateMatrixCycle(wf,wb) returns the state rate matrix of 
%   the reaction network M, given the forward rate wf and backward rate wb.

N = length(wf);
if N > 2
    M = diag(-(wf+wb),0)...
        + diag(wf(1:(N-1)),-1)...
        + diag(wb(2:N),1);
    
    M(1,end) = wf(end);
    M(end,1) = wb(1);
else
    M = [[-wf,wb];[wf,-wb]];
end

end