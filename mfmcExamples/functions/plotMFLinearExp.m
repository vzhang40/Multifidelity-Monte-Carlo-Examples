% This function plots all the results from the MFMC vs MC numerical
% experiment for the Exponential Linear Regression Example.
function plotMFLinearExp(models, betaMFMC, betaMC, CxyMFMC, CxyMC, exactLR, p, a, b, d, R)
    
    rbg = orderedcolors("gem");
    figure(1); clf(1);
    x = linspace(a, b, 100)';
    VDM = x.^(0:d)'; % VDM matrix

    %% First Element of Cxy versus p
    figure(1); clf(1)
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
    pl{1} = semilogx(p, exactLR.poly(5).*ones(length(p),1), "Color", [0, 0, 0, 1]);
    hold on
    mc5 = zeros(R, length(p));
    mf5 = zeros(R, length(p));
    for i = 1:length(p)
        mc5(:, i) = (5).^(0:d)*betaMC(:, : , i);
        mf5(:, i) = (5).^(0:d)*betaMFMC(:, : , i);
    end
    pl{2} = semilogx(p, mc5(1, :), "Color", [rbg(1, :) 0.1], "LineWidth", 0.01);
    semilogx(p, mc5(2:end, :), "Color", [rbg(1, :) 0.1], "LineWidth", 0.01);
    pl{3} = semilogx(p, mf5(1,:), "Color", [rbg(2, :) 0.1], "LineWidth", 0.01);
    semilogx(p, mf5(2:end, :), "Color", [rbg(2, :) 0.1], "LineWidth", 0.01);
    semilogx(p, exactLR.poly(5).*ones(length(p),1), "Color", [0, 0, 0, 1], "DisplayName", "True Value");
    pl{4} = semilogx(p, mean(mc5, 1)', "Color", [rbg(1, :) 1], "Linestyle", "--", "Marker", "o", "DisplayName", "Monte Carlo");
    pl{5} = semilogx(p, mean(mf5, 1)', "Color", [rbg(2, :) 1], "Linestyle", "--", "Marker", "x", "DisplayName", "Multi-fidelity Monte Carlo");
    xlabel("x")
    ylabel("y")
    title("Predicting $f^{(1)}(5)$", "Interpreter", "latex")
    
    legend([pl{1}, pl{4}, pl{5}], 'True Value','Monte Carlo', 'Multi-fidelity MC', 'Location', 'best');

    if ~isfolder("Plots")
         mkdir("Plots");
    end

    if ~exist('.../Plots/expLR1.png', 'file')
        f1 = fullfile("Plots", "expLR1.png"); 
        saveas(gcf, f1)
    end

    %% Regression Plots for Each Computational Budget
    figure(2); clf(2);
    xlim([0,5])
    for i = 1:length(p)
        subplot(length(p), 1, i)
        plot(x, models(1).f(x), "Color", [0, 0, 0, 1], "DisplayName", "True Value")
        hold on
        for j = 1:R
            fhat = VDM'*betaMFMC(:, j, i);
            plot(x, fhat, "Color", [rbg(2, :) 0.5], "LineWidth", 0.01)
            fhat = VDM'*betaMC(:, j, i);
            plot(x, fhat, "Color", [rbg(1, :) 0.5], "LineWidth", 0.01)
        end
        xlabel("x", "Interpreter", "latex")
        ylabel("y", "Interpreter", "latex")
        title("Regression MF versus hi-fidelity; p = " + p(i), "Interpreter", "latex")
        plot(x, models(1).f(x), "Color", [0, 0, 0, 1], "DisplayName", "True Value")
        if i == length(p)
            legend("True Value", "Multifidelity Regression", "Hi-Fidelity Linear Regression", "Interpreter", "latex")
        end
    end

    if ~exist('.../Plots/expLR2.png', 'file')
        f2 = fullfile("Plots", "expLR2.png"); 
        saveas(gcf, f2)
    end

    %% 1000 Testing Points
    xTest = linspace(a, b, 1000)';
    XTest = repmat((xTest.^(0:d)')', [1, 1, length(p)]);
    bestf = repmat(exactLR.poly(xTest), [1, R, length(p)]);
    
    fMF = pagemtimes(XTest, betaMFMC); 
    fMC = pagemtimes(XTest, betaMC); 
    
    errorsMF = abs((fMF - bestf)./bestf);
    errorsMC = abs((fMC - bestf)./bestf);
    
    errorMF = squeeze(mean(mean(errorsMF))); 
    errorMC =squeeze(mean(mean(errorsMC))); 
    
    stdMF = squeeze(std(mean(errorsMF, 1)));
    stdMC = squeeze(std(mean(errorsMC, 1)));
    
    %% Plotting convergence of MF Linear Regression versus HF Linear Regression with computational budget
    figure(3); clf(3);
    xscale('log') 
    hold on
    plot(p, zeros(length(p),1), "k-")
    plot(p, errorMF, 'Color', [0.8500 0.3250 0.0980], "LineStyle","-")
    plot(p, errorMC, 'Color', [0 0.4470 0.7410], "LineStyle","-")
    patch([p; flip(p)], [errorMF-stdMF; flip(errorMF+stdMF)], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none')
    patch([p; flip(p)], [errorMC-stdMC; flip(errorMC+stdMC)], 'b', 'FaceAlpha', 0.1, 'EdgeColor','none')
    legend("$\hat{f}(z, \beta^*)$", "$\hat{f}(z, \beta^{MF})$", "$\hat{f}(z, \beta^{HF})$", "$\pm \sigma_{MF}$", "$\pm \sigma_{HF}$", "Interpreter", "latex")
    xlabel("Computational Budget", "Interpreter", "latex")
    ylabel("Mean Relative Error", "Interpreter", "latex")
    title("Analytical MFMC Linear Regression Example Results", "Interpreter", "latex")

    if ~exist('.../Plots/expLR3.png', 'file')
        f3 = fullfile("Plots", "expLR3.png"); 
        saveas(gcf, f3)
    end
end