% This function uses the statistical values to calculate the weights and
% sample numbers needed for MFMC
function [mfmc, mc] = getParameters(stats, w, p)
    % IN:
    %   stats: structure with fields:
    %       sigma: standard deviations of models
    %       rho: correlation coefficients with high fidelity model
    %   w: costs of models
    %   p: computational budget
    % OUT: 
    %   mfmc: structure with fields:
    %       alpha: weights associated with surrogate models
    %       m: sample number associated with models/cost
    %   mc: structure with field:
    %       m: sample number associated with cost
    rho = [stats.rho; 0]; % setting rho_{1, k+1} = 0
    r = sqrt(w(1).*(rho(1:end-1).^2 - rho(2:end).^2)./(w.*(1-rho(2).^2)));
    m1 = p./(w'*r);
    mfmc.m = floor(m1.*r');
    mfmc.alpha = (rho(1:end-1)./stats.sigma).*stats.sigma(1);
    mc.m = floor(p./w(1));
end