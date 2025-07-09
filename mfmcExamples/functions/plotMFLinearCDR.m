% This function plots all the results from the MFMC vs MC numerical
% experiment for the CDR Linear Regression Example.
function plotMFLinearCDR(betaMFMC, betaMC, CxyMFMC, CxyMC, exactLR, p, d, R, B) 
%% Cxy(1) versus p 
    rbg = orderedcolors("gem");
    figure(1); clf(1)
    pl{1} = semilogx(p, exactLR.Cxy(1).*ones(length(p),1), "Color", [0, 0, 0, 1]);
    hold on
    Cxy1MC = reshape(CxyMC(1,:,:), [R, length(p)]);
    Cxy1MF = reshape(CxyMFMC(1,:,:), [R, length(p)]);
    pl{2} = semilogx(p, Cxy1MC(1, :), "Color", [rbg(1, :) 0.1], "LineWidth", 0.01);
    semilogx(p, Cxy1MC(2:end, :), "Color", [rbg(1, :) 0.1], "LineWidth", 0.01);
    pl{3} = semilogx(p, Cxy1MF(1, :), "Color", [rbg(2, :) 0.1], "LineWidth", 0.01);
    semilogx(p, Cxy1MF(2:end, :), "Color", [rbg(2, :) 0.1], "LineWidth", 0.01);
    semilogx(p, exactLR.Cxy(1).*ones(length(p),1), "Color", [0, 0, 0, 1])
    pl{4} = semilogx(p, mean(Cxy1MC, 1), "Color", [rbg(1, :) 1], "Linestyle", "--", "Marker", "o");
    pl{5} = semilogx(p, mean(Cxy1MF, 1), "Color", [rbg(2, :) 1], "Linestyle", "--", "Marker", "x");
    xlabel("x")
    ylabel("y")
    title("First Element of $\hat{c}_{XY}$; Average value of $f^{(1)}$", "Interpreter", "latex")
    legend([pl{1}, pl{2}, pl{3}, pl{4}, pl{5}], 'True Value','Monte Carlo', 'Multi-fidelity MC', 'Monte Carlo Mean', 'Multi-fidelity MC Mean');
    
    %% Generates Plots and Table Values for MFMC vs MC Results in Linear Regression
    % Testing points
    [XTest, ~] = getX(B(:, :), d);
    bestf = repmat(XTest'*exactLR.beta, [1, R, length(p)]);
    XTest = repmat(XTest', [1, 1, length(p)]);
    
    fMF = pagemtimes(XTest, betaMFMC); 
    fMC = pagemtimes(XTest, betaMC); 
    
    errorsMF = abs((fMF - bestf)./bestf);
    errorsMC = abs((fMC - bestf)./bestf);
    
    errorMF = squeeze(mean(mean(errorsMF))); 
    errorMC =squeeze(mean(mean(errorsMC))); 
    
    stdMF = squeeze(std(mean(errorsMF, 1)));
    stdMC = squeeze(std(mean(errorsMC, 1)));

    if ~exist('.../Plots/cdr1.png', 'file')
        f1 = fullfile("Plots", "cdr1.png"); 
        saveas(gcf, f1)
    end
    
    %% Plotting convergence of MF Linear Regression versus HF Linear Regression with computational budget
    figure(2); clf(2);
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
    title("CDR MFMC Linear Regression Example Results", "Interpreter", "latex")

    if ~exist('.../Plots/cdr2.png', 'file')
        f2 = fullfile("Plots", "cdr2.png"); 
        saveas(gcf, f2)
    end
end