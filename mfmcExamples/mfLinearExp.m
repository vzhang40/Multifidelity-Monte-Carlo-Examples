clear; close all; clc;

%% Problem Set-up
addpath('..\mfmcExamples\functions')

% Models
models(1).f = @(x) exp(x); % high-fidelity
models(2).f = @(x) -0.9*exp(0.5.*x); % low fidelity
a = 0; b = 5; % bounds

% Model Costs
w = [1; 0.001]; % [high-fidelity, low-fidelity]

% Computational Budgets
p = [10; 100; 1000];

% number of replicates (how many estimators)
R = 500; 

%% Getting scalar stats values
N_gauss = 10;
stats = getStatsGauss(models, N_gauss, a, b);
mu_true = stats.mus(1);

%% Polynomial Set-up, d is the polynomial order
d = 4;

%% Calculating Exact Analytical Values for Cxx and Cxy
% Gaussian Quadrature
exactLR = getExactLR(models(1).f, d, a, b);

%% Getting Parameter Values
[mfmc, mc] = getParameters(stats, w, p); 

%% MC/MFMC Linear Regression
% beta will be a d+1 x R matrix with pages associated with each cost
% For each cost there is R replicates of the linear regression
betaMFMC = zeros(d+1, R, length(p));
betaMC = zeros(d+1, R, length(p));
CxyMFMC = zeros(d+1, R, length(p));
CxyMC = zeros(d+1, R, length(p));

for i = 1:length(p)
    [betaMC(:, :, i), CxyMC(:, :, i)] = doMCLR(models(1).f, exactLR.Cxx, d, mc.m(i), R);
    [betaMFMC(:, :, i), CxyMFMC(:, :, i)] = doMFMCLR(models, mfmc, exactLR.Cxx, d, R, i);
end

%% Plotting:
plotMFLinearExp(models, betaMFMC, betaMC, CxyMFMC, CxyMC, exactLR, p, a, b, d, R)




