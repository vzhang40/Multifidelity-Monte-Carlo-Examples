% This function uses the statistical values to calculate the weights and
% sample numbers needed for MFMC
function [mfmc, mc] = getParameters(stats, w, p)
    % IN:
    %   stats.sigma: standard deviations of models
    %   stats.rho: correlation coefficients with high fidelity model
    %   w: costs of models
    %   p: computational budget
    % OUT: 
    %   mfmc: field with weights (alpha) and sample number (m)
    %   mc: field with sample number (m)
    rho = [stats.rho; 0]; % setting rho_{1, k+1} = 0
    r = sqrt(w(1).*(rho(1:end-1).^2 - rho(2:end).^2)./(w.*(1-rho(2).^2)));
    m1 = p./(w'*r);
    mfmc.m = floor(m1.*r');
    mfmc.alpha = (rho(1:end-1)./stats.sigma).*stats.sigma(1);
    mc.m = p/w(1);
end