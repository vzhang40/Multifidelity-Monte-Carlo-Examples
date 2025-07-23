clear; close all; clc
% This script implements Multi-fidelity Monte Carlo (MFMC) method to estimate the
% mean value of the maximum temperature in a Convection-Diffusion-Reaction (CDR) Experiment. 
% The results compare the performance of Multi-fidelity methods compared to using only high 
% fidelity models to estimate.

%% Problem Set-up
addpath('..\mfmcExamples\functions')
load('samples.mat') % Loads CDR data into the workspace

% Models
models.f = yA; % Column 1: Hi-Fidelity, Column 2: Surrogate/Low Fidelity

clearvars -except models

% Cost vector 
w = [1.94; 6.20e-3]; % [Hi-Fidelity, Low Fidelity]

% Computational Budgets
p = [10; 100; 1000];

% Number of Replicates
R = 500;

%% Getting necessary statistics for analysis
% Statistics are generated with full access to the data
stats = getStats(models.f);
mu_true = stats.mus(1); % Expected value of truth function

%% Numerical Experiment: Sampling
% Initializing variables
s_MFMC = zeros(R, length(p));
s_MC = zeros(R, length(p));

% Getting parameters to calculate MFMC and MC estimator
[mfmc, mc] = getParameters(stats, w, p);

for n = 1:length(p)
    % Implementing MFMC
    s_MFMC(:, n) = doMFMC(models, mfmc, R, n, "cdr");
    
    % Implementing MC
    s_MC(:, n) = doMC(models.f(:, 1), mc.m(n), R, "cdr");
end

%% Experimental and Analytical Statistics
[mfmc, mc] = getExperimentAnalysisStats(mfmc, mc, s_MFMC, s_MC, stats, p, w, mu_true);

%% Generates Plots and Table Values for MFMC vs MC Results
% Table compares Analytical versus Experimental Statistics
% Figure 1: Plotting MFMC and MC convergence to true mean: cdr1.png
% Figure 2: Plotting Analytical Mean Square Error with Observed values:
% cdr2.png
showResultsMFMC(mu_true, p, mfmc, mc, "cdr")
