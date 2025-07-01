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
models.testing = yB(:, 1);

% Cost vector 
w = [1.94; 6.20e-3]; % [Hi-Fidelity, Low Fidelity]

% Computational Budgets
p = [10; 100; 1000];

% number of replicates (how many estimators)
R = 50;

%% Getting necessary statistics for analysis
% Statistics are generated with full access to the data
stats = getStats(models.f);
mu_true = stats.mus(1); % Expected value of truth function

%% Polynomial Set-up, d is the max polynomial order
d = 2;

