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
xTest = linspace(a, b, 100)';
VDM = xTest.^(0:d)'; % VDM matrix 

%% One instance of a MC linear regression for each cost
figure(1); clf(1);
xlim([0,5])
rbg = orderedcolors("gem");
plot(xTest, models(1).f(xTest), "Color", [0, 0, 0, 1], "DisplayName", "True Value")
hold on
for i = 1:length(p)
    fhat = VDM'*betaMC(:, 1, i);
    plot(xTest, fhat, "Color", rbg(i, :), "DisplayName", "p = " + p(i));
end
legend("Location", "best", "Interpreter", "latex")
xlabel("x")
ylabel("y")
title("One Instance of MCLR for each cost")

%% One instance of a MFMC linear regression for each cost
figure(2); clf(2);
xlim([0,5])
rbg = orderedcolors("gem");
plot(xTest, models(1).f(xTest), "Color", [0, 0, 0, 1], "DisplayName", "True Value")
hold on
for i = 1:length(p)
    fhat = VDM'*betaMFMC(:, 1, i);
    plot(xTest, fhat, "Color", rbg(i, :), "DisplayName", "p = " + p(i));
end
legend("Location", "best", "Interpreter", "latex")
xlabel("x")
ylabel("y")
title("One Instance of MFMCLR for each cost")

%% All of the MC replicates for each cost 
figure(3); clf(3);
plot(xTest, models(1).f(xTest), "Color", [0, 0, 0, 1], "DisplayName", "True Value")
xlim([0,5])
hold on
for j = 1:R
    for i = 1:length(p)
        fhat = VDM'*betaMC(:, j, i);
        plot(xTest, fhat, "Color", [rbg(i, :) 0.5], "LineWidth", 0.01)
    end
end
plot(xTest, models(1).f(xTest), "Color", [0, 0, 0, 1])
legend("True Value", "p = 10", "p = 100", "p = 1000")
xlabel("x")
ylabel("y")
title("All MC LR replicates for each cost")

%% All of the MFMC replicates for each cost 
figure(4); clf(4);
plot(xTest, models(1).f(xTest), "Color", [0, 0, 0, 1], "DisplayName", "True Value")
xlim([0,5])
hold on
for j = 1:R
    for i = 1:length(p)
        fhat = VDM'*betaMFMC(:, j, i);
        plot(xTest, fhat, "Color", [rbg(i, :) 0.5], "LineWidth", 0.01)
    end
 
end
plot(xTest, models(1).f(xTest), "Color", [0, 0, 0, 1])
legend("True Value", "p = 10", "p = 100", "p = 1000")
xlabel("x")
ylabel("y")
title("All MFMC LR replicates for each cost")

%% First Element of Cxy versus p
figure(5); clf(5)
subplot(3, 1, 1)
semilogx(p, exactLR.Cxy(1).*ones(length(p),1), "Color", [0, 0, 0, 1])
hold on
Cxy1MC = reshape(CxyMC(1,:,:), [R, length(p)]);
Cxy1MF = reshape(CxyMFMC(1,:,:), [R, length(p)]);
semilogx(p, Cxy1MC, "Color", [rbg(1, :) 0.1], "LineWidth", 0.01)
semilogx(p, Cxy1MF, "Color", [rbg(2, :) 0.1], "LineWidth", 0.01)
semilogx(p, exactLR.Cxy(1).*ones(length(p),1), "Color", [0, 0, 0, 1])
semilogx(p, mean(Cxy1MC, 1), "Color", [rbg(1, :) 1], "Linestyle", "--", "Marker", "o")
semilogx(p, mean(Cxy1MF, 1), "Color", [rbg(2, :) 1], "Linestyle", "--", "Marker", "x")
xlabel("x")
ylabel("y")
title("First Element of $\hat{c}_{XY}$; Average value of $f^{(1)}$", "Interpreter", "latex")

%% First Element of beta versus p
subplot(3, 1, 2)
semilogx(p, exactLR.beta(1).*ones(length(p),1), "Color", [0, 0, 0, 1])
hold on
beta1MC = reshape(betaMC(1,:,:),  [R, length(p)]);
beta1MF = reshape(betaMFMC(1,:,:), [R, length(p)]);
semilogx(p, beta1MC, "Color", [rbg(1, :) 0.1], "LineWidth", 0.01)
semilogx(p, beta1MF, "Color", [rbg(2, :) 0.1], "LineWidth", 0.01)
semilogx(p, exactLR.beta(1).*ones(length(p),1), "Color", [0, 0, 0, 1])
semilogx(p, mean(beta1MC, 1)', "Color", [rbg(1, :) 1], "Linestyle", "--", "Marker", "o")
semilogx(p, mean(beta1MF, 1)', "Color", [rbg(2, :) 1], "Linestyle", "--", "Marker", "x")
xlabel("x")
ylabel("y")
title("First Element of $\hat{\beta}$", "Interpreter", "latex")

%% fhat(5) versus p
subplot(3, 1, 3)
semilogx(p, exactLR.poly(5).*ones(length(p),1), "Color", [0, 0, 0, 1])
hold on
mc5 = zeros(R, length(p));
mf5 = zeros(R, length(p));
for i = 1:length(p)
    mc5(:, i) = (5).^(0:d)*betaMC(:, : , i);
    mf5(:, i) = (5).^(0:d)*betaMFMC(:, : , i);
end
semilogx(p, mc5, "Color", [rbg(1, :) 0.1], "LineWidth", 0.01)
semilogx(p, mf5, "Color", [rbg(2, :) 0.1], "LineWidth", 0.01)
semilogx(p, exactLR.poly(5).*ones(length(p),1), "Color", [0, 0, 0, 1], "DisplayName", "True Value")
semilogx(p, mean(mc5, 1)', "Color", [rbg(1, :) 1], "Linestyle", "--", "Marker", "o", "DisplayName", "Monte Carlo")
semilogx(p, mean(mf5, 1)', "Color", [rbg(2, :) 1], "Linestyle", "--", "Marker", "x", "DisplayName", "Multi-fidelity Monte Carlo")
xlabel("x")
ylabel("y")
title("Predicting $f^{(1)}(5)$", "Interpreter", "latex")

% pl{1} = plot(nan, "Color", [0, 0, 0, 1]);
% pl{2} = plot(nan, "Color", [rbg(1, :) 0.1]);
% pl{3} = plot(nan, "Color", [rbg(2, :) 0.1]);
% pl{4} = plot(nan, "Color", [rbg(1, :) 1], "Linestyle", "--", "Marker", "o");
% pl{5} = plot(nan, "Color", [rbg(2, :) 1], "Linestyle", "--", "Marker", "x");

% legend([pl{:}], {'True Value','Monte Carlo','Multi-fidelity MC', 'Monte Carlo Mean', 'Multi-fidelity MC Mean'}, 'location', 'southoutside')