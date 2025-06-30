function plotMSE(p, mfmc, mc, fig)
    figure(fig); clf(fig);
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