% This function uses the statistical values to calculate the weights and
% sample numbers needed for MFMC
function [alpha, m_MFMC, m_MC] = getParameters(sigma, rho, w, p)
    % IN:
    %   sigma: standard deviations of models
    %   rho: correlation coefficients with high fidelity model
    %   w: costs of models
    %   p: computational budget
    % OUT: 
    %   alpha: weight for models
    rho = [rho; 0]; % setting rho_{1, k+1} = 0
    r = sqrt(w(1).*(rho(1:end-1).^2 - rho(2:end).^2)./(w.*(1-rho(2).^2)));
    m1 = p./(w'*r);
    m_MFMC = floor(m1.*r);
    alpha = (rho(1:end-1)./sigma).*sigma(1);
    m_MC = p/w(1);
end