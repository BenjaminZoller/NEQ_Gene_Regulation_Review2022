function [M] = makeRateMatrixBinding(kb,ku,kq)
%MAKERATEMATRIXBINDING makes the state rate matrix (the Laplacian matrix of
%the master equation) for the proof-reading model in Fig. 2 left.
%   [M] = makeRateMatrixBinding(kb,ku,kq) returns the state rate matrix of 
%   the reaction network M, given the TF binding rate kb, the TF unbinding
%   rate ku, and the proof-reading rate kq.

M = [[-ku,kq,0];[0,-(kq+ku),kb];[ku,ku,-kb]];

end

