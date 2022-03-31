function [M] = makeRateMatrixBinding(kb,ku,kq)
%MAKERATEMATRIXBINDING makes the state rate matrix (the Laplacian matrix of
%the master equation) for the proof-reading model in Fig. 2 left.
%   [M] = makeRateMatrixBinding(kb,ku,kq) returns the state rate matrix of 
%   the reaction network M, given the TF binding rate kb, the TF unbinding
%   rate ku, and the proof-reading rate kq.
%
%   Copyright (c) 2022, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree.

M = [[-ku,kq,0];[0,-(kq+ku),kb];[ku,ku,-kb]];

end

