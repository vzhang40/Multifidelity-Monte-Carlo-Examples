clear; close all; clc;
% This script implements Multi-fidelity Monte Carlo (MFMC) method to estimate the
% value of E[XY] in order to evalute the coefficients of a polynomial regression of 
% max degree d. The results compare the performance of Multi-fidelity methods 
% compared to using only high fidelity methods to estimate a function for
% the maximum temperature in the chamber of a CDR model with inputs and
% parameters (A in samples.mat):
% Ti: inlet temperature
% To: left wall temperature
% phi: fuel:oxidizer ratio of premixed inflow
% Ape: pre-exponential factor of Arrhenius Equation
% E: activation energy
 
%% Problem Set-up
load('samples.mat')
addpath('..\mfmcExamples\functions')

models.inputs = A;
models.f = yA;
models.testInputs = B;
models.testf = yB(:, 1);

clearvars -except models

% Cost vector 
w = [1.94; 6.20e-3]; % [Hi-Fidelity, Low Fidelity]

% Computational Budgets
p = [10; 100; 1000];

% number of replicates (how many estimators)
R = 500;

%% Getting necessary statistics for analysis
% Statistics are generated with full access to the data
stats = getStats(models.f);
mu_true = stats.mus(1); % Expected value of truth function

%% Polynomial Set-up
% max polynomial order
d = 2;

%% Calculating Exact Analytical Values for Cxx and Cxy
% Gaussian Quadrature
exactLR = getExactLR(models.f(:, 1), d, "cdr", models.inputs);

%% Getting Parameter Values, \alpha^{mean}
[mfmc, mc] = getParameters(stats, w, p); 

%% Numerical Experiment: Sampling
% Initializing variables
% beta will be a (# of features) x R matrix with pages associated with each cost
% For each cost there is R replicates of the linear regression
nFeat = size(exactLR.Cxx, 1);
betaMFMC = zeros(nFeat, R, length(p));
betaMC = zeros(nFeat, R, length(p));
CxyMFMC = zeros(nFeat, R, length(p));
CxyMC = zeros(nFeat, R, length(p));

%% MC/MFMC Linear Regression
for i = 1:length(p)
    [betaMC(:, :, i), CxyMC(:, :, i)] = doMCLR(models.f(:, 1), exactLR.Cxx, d, mc.m(i), R, "cdr", models.inputs);
    [betaMFMC(:, :, i), CxyMFMC(:, :, i)] = doMFMCLR(models.f, mfmc, exactLR.Cxx, d, R, i, "cdr", models.inputs);
end

%% Generates Plots for MFMC vs MC Results in Linear Regression
% Figure 1: This plots the first element of \hat{c}_{XY} or the
%   estimated value of XY with computational cost to analyze convergence.
% Figure 2: Figure 3: Plotting MFMC and MC convergence to best fit polynomial
plotMFLinearCDR(betaMFMC, betaMC, CxyMFMC, CxyMC, exactLR, p, d, R, models.testInputs) 



