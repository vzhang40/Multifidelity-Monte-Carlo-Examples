clear; close all; clc
% This script implements Multi-fidelity Monte Carlo (MFMC) method to estimate the
% mean value of exp(x). The results compare the performance of
% Multi-fidelity methods compared to using only high fidelity methods to
% estimate.

%% Problem Set-up
addpath('..\mfmcExamples\functions') 

% Models
models(1).f = @(x) exp(x); % high-fidelity
models(2).f = @(x) 0.9*exp(0.5.*x); % low fidelity

a = 0; b = 5; % bounds

% Model Costs
w = [1; 0.001]; % [high-fidelity, low-fidelity]

% Computational Budgets
p = [10; 100; 1000];

% Number of Replicates
R = 500; 

%% Getting Necessary Statistics for Analysis
% Gaussian Quadrature Calculations
N_gauss = 10; % number of points sampled for stats, accurate up to 19-order polynomial.
stats = getStatsGauss(models, N_gauss, a, b);
mu_true = stats.mus(1); % Expected value of truth function

%% Numerical Experiment: Sampling
% Initializing variables
s_MFMC = zeros(R, length(p)); % MFMC estimates
s_MC = zeros(R, length(p)); % MC estimates

% Getting parameters to calculate MFMC and MC estimator
[mfmc, mc] = getParameters(stats, w, p);

for n = 1:length(p)
    % Implementing MFMC
    s_MFMC(:, n) = doMFMC(models, mfmc, R, n, "exp");
    
    % Implementing MC
    s_MC(:, n) = doMC(models(1).f, mc.m(n), R, "exp");
end

%% Experimental and Analytical Statistics
[mfmc, mc] = getExperimentAnalysisStats(mfmc, mc, s_MFMC, s_MC, stats, p, w, mu_true);

%% Generates Plots and Table Values for MFMC vs MC Results
% Table compares Analytical versus Experimental Statistics
% Figure 1: Plotting MFMC and MC convergence to true mean: exp1.png
% Figure 2: Plotting Analytical Mean Square Error with Observed values:
% exp2.png
showResultsMFMC(mu_true, p, mfmc, mc, "exp")





