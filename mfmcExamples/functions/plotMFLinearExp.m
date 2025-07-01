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
    
    legend([pl{1}, pl{2}, pl{3}, pl{4}, pl{5}], 'True Value','Monte Carlo', 'Multi-fidelity MC', 'Monte Carlo Mean', 'Multi-fidelity MC Mean', 'Location', 'best');

    %% 
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
        xlabel("x")
        ylabel("y")
        title("Regression MF versus hi-fidelity; p = " + p(i))
        plot(x, models(1).f(x), "Color", [0, 0, 0, 1], "DisplayName", "True Value")
        if i == length(p)
            legend("True Value", "Multifidelity Regression", "Hi-Fidelity Linear Regression")
        end
    end
end