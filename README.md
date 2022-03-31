# NEQ_Gene_Regulation_Review2022
Simulation code for toy-models

The Matlab sripts and functions are provided to reproduce all the
figures and explore the regulatory phenotypes of the four toy-models
illustrated in the review.

The code has been tested on MATLAB R2020a.

One can find in the repository the following piece of codes:
1. Functions to build the state rate matrix (Laplacian matrix of the
    master equation) for each model as a function of input parameters.
    - makeRateMatrixCycle.m – Fig. 1
    - makeRateMatrixBinding.m – Fig. 2 left
    - makeRateMatrixCooperativity.m – Fig. 2 right
    - makeRateMatrixNeqMWC.m – Fig. 3
    - makeRateMatrixMWC.m – Fig. 3
2. Generic functions to compute the various regulatory phenotypes and 
    other features (waiting time distributions, fluxes, entropy production,
    etc.) from the state rate matrix.
    - computePheno.m
    - getExp.m
    - getNoise.m
    - getResid.m
    - getEntropy.m
3. Two functions to generate stochastic realizations of the models.
    - genStoTrajectories.m
    - genRastTrajectories.m
4. Scripts to generate the different figures.
    - MainToyModelCycle.m – Fig. 1
    - MainToyModelBinding.m – Fig. 2 left
    - MainToyModelCooperativity.m – Fig. 2 right
    - MainToyModelNeqMWC.m – Fig. 3

%   Copyright (c) 2022, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree.
