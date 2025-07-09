clear; close all; clc;
% This script implements Multi-fidelity Monte Carlo (MFMC) method to estimate the
% value of E[XY] in order to evalute the coefficients of a polynomial regression of 
% max degree d. The results compare the performance of Multi-fidelity methods 
% compared to using only high fidelity methods to estimate f = exp(x)

%% Problem Set-up
addpath('..\mfmcExamples\functions')

% Models
models(1).f = @(x) exp(x); % high-fidelity
models(2).f = @(x) -0.9*exp(0.5.*x); % low fidelity
a = 0; b = 5; % bounds

% Model Costs
w = [1; 0.001]; % [Hi-Fidelity, Low Fidelity]

% Computational Budgets
p = [10; 100; 1000];

% Number of Estimates
R = 500; 

%% Getting Necessary Statistics for Analysis
% Gaussian Quadrature Calculations for scalar statistical values of models
N_gauss = 10;
stats = getStatsGauss(models, N_gauss, a, b);
mu_true = stats.mus(1); % Expected value of truth function

%% Polynomial Set-up, d is the max polynomial order; (d+1) features in X
d = 4;

%% Calculating Exact Analytical Values for Cxx and Cxy
% Gaussian Quadrature
exactLR = getExactLR(models(1).f, d, "exp", a, b);

%% Getting Parameter Values, \alpha^{mean}
[mfmc, mc] = getParameters(stats, w, p); 

%% Numerical Experiment: Sampling
% Initializing variables
% beta will be a d+1 x R matrix with pages associated with each cost
% For each cost there is R replicates of the linear regression
betaMFMC = zeros(d+1, R, length(p));
betaMC = zeros(d+1, R, length(p));
CxyMFMC = zeros(d+1, R, length(p));
CxyMC = zeros(d+1, R, length(p));

% MC/MFMC Linear Regression
for i = 1:length(p)
    [betaMC(:, :, i), CxyMC(:, :, i)] = doMCLR(models(1).f, exactLR.Cxx, d, mc.m(i), R, "exp");
    [betaMFMC(:, :, i), CxyMFMC(:, :, i)] = doMFMCLR(models, mfmc, exactLR.Cxx, d, R, i, "exp");
end

%% Generates Plots for MFMC vs MC Results in Linear Regression
% Figure 1: Generates 3 subplots
%   Subplot 1: This plots the first element of \hat{c}_{XY} or the
%   estimated value of XY with computational cost to analyze convergence.
%   Subplot 2: This plots the first element of \hat{beta} with
%   computational cost
%   Subplot 3: This plots the value of \hat{f}(5) of the polynomial
%   regression with computational cost
% Figure 2: Generates 3 subplots with multi-fidelity versus high-fidelity
% regressions at different computational budgets.
%   Subplot 1: p = 10
%   Subplot 2: p = 100
%   Subplot 3: p = 1000
% Figure 3: Plotting MFMC and MC convergence to best fit polynomial
plotMFLinearExp(models, betaMFMC, betaMC, CxyMFMC, CxyMC, exactLR, p, a, b, d, R)




