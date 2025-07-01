% This function updates the mfmc and mc structures to be updated with
% experimental statistical results and the analytical statistical results
% using the equations in the MFMC paper.
function [mfmc, mc] = getExperimentAnalysisStats(mfmc, mc, s_MFMC, s_MC, stats, p, w, mu_true)
    % IN:
    %   mfmc: structure with mfmc parameters to be updated
    %   mc: structure with mc parameters to be updated
    %   s_MFMC: matrix with R replicates of MFMC results for each computational budget
    %   s_MC: matrix with R replicates of MC results for each computational budget
    %   stats: structre with true statistical values fields from the models:
    %       sigma: standard deviation of each model
    %       rho: correlation coefficient of each model with the truth model
    %   p: computational budget vector
    %   w: cost vector
    %   mu_true: true mean value for truth model
    % 
    % OUT:
    %   mfmc: structure with new fields associated with the mfmc method:
    %       means: experimental means for each computational budget
    %       sigma: experimental standard deviations for each computational
    %       budget
    %       mse: experimental mean square error for each computational
    %       budget
    %       sigmaAnal: analytical standard deviations for each
    %       computational budget
    %       mseAnal: analytical mean square errors for each
    %       computational budget
    %   mc: structure with new fields associated with the mc method:
    %       same as mfmc above

    % Experimental statistics from numerical experiment results
    mfmc.means = mean(s_MFMC)';
    mfmc.sigma = std(s_MFMC)';
    mc.means = mean(s_MC)';
    mc.sigma = std(s_MC)';
    mfmc.mse = mean((mu_true - s_MFMC).^2)';
    mc.mse = mean((mu_true - s_MC).^2)';
    
    % Analytical statistics, equations are from MFMC paper
    mfmc.mseAnal = stats.sigma(1)^2*(1-stats.rho(2)^2)*p./(mfmc.m(:,1).^2*w(1));
    mc.mseAnal = stats.sigma(1)^2./mc.m;
    mc.sigmaAnal = stats.sigma(1)./sqrt(mc.m);
    mfmc.sigmaAnal = sqrt(stats.sigma(1)^2./mfmc.m(:, 1) - sum((1./mfmc.m(:, 1:end-1) - 1./mfmc.m(:, 2:end))*stats.rho(2:end).^2.*stats.sigma(1)^2, 2));
end