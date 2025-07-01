clear; close all; clc;

%% Problem Set-up
load('samples.mat')
addpath('..\mfmcExamples\functions')

models.inputs = A;
models.f = yA;
models.testing = yB(:, 1);

% Cost vector [hi-fidelity, low fidelity]
w = [1.94; 6.20e-3];

% Computational Budgets
p = [10; 100; 1000];

% number of replicates (how many estimators)
R = 50;
