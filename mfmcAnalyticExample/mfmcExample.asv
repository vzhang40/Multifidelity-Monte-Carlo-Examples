clear; close all; clc
% This script implements all the analytical example in Elizabeth's paper,
% comparing MC and MFMC to estimate the mean value of exp(x)
% things I didn't implement:
% --> ordering the functions by correlation: theres only 2

%% Problem Set-up
% Models
models(1).f = @(x) exp(x); % high-fidelity
models(2).f = @(x) -0.9*exp(0.5.*x); % low fidelity
% models(3).f = @(x) 0.5*exp(0.4.*x); % lower fidelity/???
% models(4).f = @(x) -0.1*exp(0.1.*x);
a = 0; b = 5; % bounds

% Model Costs
w = [1; 0.001]; % [high-fidelity, low-fidelity]

% Plotting models
figure(1); clf(1);
plotModels(models, a, b)

% Computational Budgets
p = [10; 40; 100; 1000];

% number of replicates (how many estimators)
R = 50; 

%% Getting necessary statistics for analysis

% Monte Carlo Estimates for Statistics, does not work if you add more
% functions 6/16
% N_MC = 1e6; % number of points sampled for stats
% stats = getStatsMC(models, N_MC);
% mu_true = stats.mus(1);

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
    s_MFMC(:, n) = doMFMC(models, mfmc, R, n);
    
    % Implementing MC
    s_MC(:, n) = doMC(models(1).f, mc.m(n), R);
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
figure(2); clf(2)
plotEstimates(p, mu_true, mfmc, mc)

%% Plotting Computational Cost versus MSE
figure(3); clf(3);
plotMSE(p, mfmc, mc)

%% ~~~~~~~~~~~~~~~~~  Plotting functions  ~~~~~~~~~~~~~~~~~~~~~
function plotModels(models, a, b)
    x = linspace(a, b, 100);
    legendtext = cell(size({models.f}));
    legendtext{1} = "Truth $f^{(1)}$";
    plot(x, models(1).f(x))
    hold on
    for i = 2:length({models.f})
        plot(x, models(i).f(x))
        legendtext{i} = "$f^{(" + i + ")}$";
    end
    legend(legendtext, "Interpreter", "latex")
    xlabel("$x$", "Interpreter", "latex")
    ylabel("$f^{(i)}(x)$", "Interpreter", "latex")
    title("Truth and Surrogate Models")
end

function showTableVals(p, mfmc, mc)
    values = table;
    values.Computational_Budget = p;
    values.MFMC_Means = mfmc.means;
    values.MFMC_Std_Dev = mfmc.sigma;
    values.MFMC_Std_Dev_Analytical = mfmc.sigmaAnal;
    values.MFMC_MSE = mfmc.mse;
    values.MFMC_MSE_Analytical = mfmc.mseAnal;
    values.MC_Means = mc.means;
    values.MC_Std_Dev = mc.sigma;
    values.MC_Std_Dev_Analytical = mc.sigmaAnal;
    values.MC_MSE = mc.mse;
    values.MC_MSE_Analytical = mc.mseAnal;
    disp(values)
end

function plotEstimates(p, mu_true, mfmc, mc)
    xscale('log') 
    hold on
    plot(p, mu_true.*ones(length(p),1), "k-")
    plot(p, mfmc.means, 'Color', [0.8500 0.3250 0.0980], "LineStyle","-")
    plot(p, mc.means, 'Color', [0 0.4470 0.7410], "LineStyle","-")
    patch([p; flip(p)], [mfmc.means-mfmc.sigmaAnal; flip(mfmc.means+mfmc.sigmaAnal)], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none')
    patch([p; flip(p)], [mc.means-mc.sigmaAnal; flip(mc.means+mc.sigmaAnal)], 'b', 'FaceAlpha', 0.1, 'EdgeColor','none')
    legend("True $\mu$", "$\mu_{MFMC}$", "$\mu_{MC}$", "$\pm \sigma_{MFMC}$", "$\pm \sigma_{MC}$", "Interpreter", "latex")
    xlabel("Computational Budget", "Interpreter", "latex")
    ylabel("Estimate Mean Model Output", "Interpreter", "latex")
    title("Analytical MFMC Example Results", "Interpreter", "latex")
end

function plotMSE(p, mfmc, mc)
    loglog(p, mfmc.mseAnal, 'Color', [0.8500 0.3250 0.0980], "LineStyle","-")
    hold on
    loglog(p, mc.mseAnal, 'Color', [0 0.4470 0.7410], "LineStyle","-")
    loglog(p, mfmc.mse, 'Color', [0.8500 0.3250 0.0980], "Marker","o", "LineStyle","none")
    loglog(p, mc.mse, 'Color', [0 0.4470 0.7410], "Marker","o", "LineStyle","none")
    legend("Analytical MFMC", "Analytical MC", "Observed MFMC", "Observed MC", "Interpreter", "latex")
    xlabel("Computational Cost $p$", "Interpreter", "latex")
    ylabel("Mean Square Error", "Interpreter", "latex")
    title("Mean Square Error in Sampling Methods", "Interpreter", "latex")
end