% This function plots all the results from the MFMC vs MC numerical
% experiment for the Exponential and CDR Example
function showResultsMFMC(mu_true, p, mfmc, mc, flag)
    %% Showing Table of values
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

    %% Plotting Estimates
    if flag == "exp"
        name = "Exponential";
    elseif flag == "cdr"
        name = "CDR";
    end
    figure(1); clf(1);
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
    title(name +  " MFMC Example Results", "Interpreter", "latex")

    if ~isfolder("Plots")
         mkdir("Plots");
    end
    
    if flag == "exp"
        if ~exist('.../Plots/exp1.png', 'file')
        f1 = fullfile("Plots", "exp1.png"); 
        saveas(gcf, f1)
        end
    elseif flag == "cdr"
        if ~exist('.../Plots/cdr1.png', 'file')
        f1 = fullfile("Plots", "cdr1.png"); 
        saveas(gcf, f1)
        end
    end

    %% Plotting Mean Square Error
    figure(2); clf(2);
    loglog(p, mfmc.mseAnal, 'Color', [0.8500 0.3250 0.0980], "LineStyle","-")
    hold on
    loglog(p, mc.mseAnal, 'Color', [0 0.4470 0.7410], "LineStyle","-")
    loglog(p, mfmc.mse, 'Color', [0.8500 0.3250 0.0980], "Marker","o", "LineStyle","none")
    loglog(p, mc.mse, 'Color', [0 0.4470 0.7410], "Marker","o", "LineStyle","none")
    legend("Analytical MFMC", "Analytical MC", "Observed MFMC", "Observed MC", "Interpreter", "latex")
    xlabel("Computational Cost $p$", "Interpreter", "latex")
    ylabel("Mean Square Error", "Interpreter", "latex")
    title("Mean Square Error in Sampling Methods: " + name, "Interpreter", "latex")

    if flag == "exp"
        if ~exist('.../Plots/exp2.png', 'file')
        f1 = fullfile("Plots", "exp2.png"); 
        saveas(gcf, f1)
        end
    elseif flag == "cdr"
        if ~exist('.../Plots/cdr2.png', 'file')
        f1 = fullfile("Plots", "cdr2.png"); 
        saveas(gcf, f1)
        end
    end
    
end