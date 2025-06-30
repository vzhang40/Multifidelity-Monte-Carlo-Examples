clear; close all; clc
% This script implements all the analytical example in Elizabeth's paper,
% comparing MC and MFMC to estimate the mean value of exp(x)
% things I didn't implement:
% --> ordering the functions by correlation: theres only 2

%% Problem Set-up
addpath('..\mfmcExamples\functions')

% Models
models(1).f = @(x) exp(x); % high-fidelity
models(2).f = @(x) 0.9*exp(0.5.*x); % low fidelity
% models(3).f = @(x) 0.5*exp(0.4.*x); % lower fidelity/???
% models(4).f = @(x) -0.1*exp(0.1.*x);
a = 0; b = 5; % bounds

% Model Costs
w = [1; 0.001]; % [high-fidelity, low-fidelity]

% Plotting models
plotModels(models, a, b, 1)

% Computational Budgets
p = [10; 40; 100; 1000];

% number of replicates (how many estimators)
R = 50; 

%% Getting necessary statistics for analysis
% Gaussian Quadrature Calculations, works if you add more functions
N_gauss = 100; % number of points sampled for stats
stats = getStatsGauss(models, N_gauss, a, b);
mu_true = stats.mus(1);

%% Sampling
s_MFMC = zeros(R, length(p));
s_MC = zeros(R, length(p));
[mfmc, mc] = getParameters(stats, w, p);

for n = 1:length(p)
    % Getting parameters to calculate MFMC and MC estimator
    % Implementing MFMC
    s_MFMC(:, n) = doMFMC(models, mfmc, R, n, "exp");
    
    % Implementing MC
    s_MC(:, n) = doMC(models(1).f, mc.m(n), R, "exp");
end

%% Experimental and Analytical Statistics
% Experimental stats
mfmc.means = mean(s_MFMC)';
mfmc.sigma = std(s_MFMC)';
mc.means = mean(s_MC)';
mc.sigma = std(s_MC)';
mfmc.mse = mean((mu_true - s_MFMC).^2)';
mc.mse = mean((mu_true - s_MC).^2)';

% Analytical stats
mfmc.mseAnal = stats.sigma(1)^2*(1-stats.rho(2)^2)*p./(mfmc.m(:,1).^2*w(1));
mc.mseAnal = stats.sigma(1)^2./mc.m;
mc.sigmaAnal = stats.sigma(1)./sqrt(mc.m);
mfmc.sigmaAnal = sqrt(stats.sigma(1)^2./mfmc.m(:, 1) - sum((1./mfmc.m(:, 1:end-1) - 1./mfmc.m(:, 2:end))*stats.rho(2:end).^2.*stats.sigma(1)^2, 2));

% Comparing Numerical Results
showTableVals(p, mfmc, mc)

%% Plotting Comptuational Cost versus Mean Estimate
plotEstimates(p, mu_true, mfmc, mc, 2)

%% Plotting Computational Cost versus MSE
plotMSE(p, mfmc, mc, 3)
