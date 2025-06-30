function plotEstimates(p, mu_true, mfmc, mc, fig)
    figure(fig); clf(fig);
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