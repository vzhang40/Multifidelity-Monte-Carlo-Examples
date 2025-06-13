clear; close all; clc
% This script implements all the analytical example in Elizabeth's paper,
% comparing MC and MFMC to estimate the mean value of exp(x)
% things I didn't implement:
% --> ordering the functions by correlation: theres only 2

%% Problem Set-up

% Models
f1 = @(x) exp(x); % high-fidelity
f2 = @(x) 0.9*exp(0.5.*x); % low fidelity
f = @(x) [f1(x), f2(x)]; % [high-fidelity, surrogates]

% Model Costs
w = [1; 0.001]; % [high-fidelity, low-fidelity]

% Computational Budgets
p = [10; 40; 100; 1000];

% number of replicates (how many estimators)
R = 50; 

%% Getting necessary statistics for analysis

% Monte Carlo Estimates for Statistics
% N_MC = 1e6; % number of points sampled for stats
% [mu_true, sigma, rho] = getStatsMC(f, N_MC);

% Gaussian Quadrature Calculations
N_gauss = 100; % number of points sampled for stats
a = 0; b = 5; % bounds defined by r.v Z
[mu_true, sigma, rho] = getStatsGauss(f1, f, N_gauss, a, b);

%% Sampling

% initializing vectors
mfmc_means = zeros(length(p), 1);
mfmc_sigma = zeros(length(p), 1);
mc_means = zeros(length(p), 1);
mc_sigma = zeros(length(p), 1);
mfmc_sigma_anal = zeros(length(p), 1);
mc_sigma_anal = zeros(length(p), 1);

for n = 1:length(p)
    % Getting parameters to calculate MFMC and MC estimator
    [alpha, m_MFMC, m_MC] = getParameters(sigma, rho, w, p(n));
    
    % Implementing MFMC
    fs = {f1, f2};
    s_MFMC = doMFMC(fs, alpha, m_MFMC, R);
    
    % Implementing MC
    s_MC = doMC(f1, m_MC, R);
    
    % Experimental estimator statistics
    mfmc_means(n) = mean(s_MFMC);
    mfmc_sigma(n) = std(s_MFMC);
    mc_means(n) = mean(s_MC);
    mc_sigma(n) = std(s_MC);

     % Analytical standard deviations of estimates
    [mc_sigma_anal(n), mfmc_sigma_anal(n)] = calcAnalSTD(sigma(1), rho, m_MFMC, m_MC);
end
%% Plotting
figure(1); clf(1)
xscale('log') 
hold on
plot(p, mu_true.*ones(length(p),1), "k-")
plot(p, mfmc_means, 'r--')
plot(p, mc_means, 'b--')
patch([p; flip(p)], [mfmc_means-mfmc_sigma; flip(mfmc_means+mfmc_sigma)], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none')
patch([p; flip(p)], [mc_means-mc_sigma; flip(mc_means+mc_sigma)], 'b', 'FaceAlpha', 0.1, 'EdgeColor','none')
legend("True $\mu$", "$\mu_{MFMC}$", "$\mu_{MC}$", "$\pm \sigma_{MFMC}$", "$\pm \sigma_{MC}$", "Interpreter", "latex")
xlabel("Computational Budget", "Interpreter", "latex")
ylabel("Estimate Mean Model Output", "Interpreter", "latex")
title("Analytical MFMC Example Results", "Interpreter", "latex")

displayStuff = table(p, mfmc_means, mfmc_sigma, mfmc_sigma_anal, mc_means, mc_sigma, mc_sigma_anal)