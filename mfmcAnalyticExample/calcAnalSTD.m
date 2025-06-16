% This function calculates the analytical standard deviations of the MFMC
% estimator and the MC estimator.
function [sigma_MC, sigma_MFMC] = calcAnalSTD(sigma1, rho, m_MFMC, m_MC)
    % IN:
    %   sigma1: standard deviation of the high-fidelity model
    %   rho: vector of correlation coefficients
    %   m_MFMC: sample sizes for MFMC
    %   m_MC: number of samples for MC
    % OUT:
    %   sigma_MC: the analytical standard deviation for the MC estimator
    %   sigma_MFMC: the analytical standard deviation for the MFMC
    %               estimator
    sigma_MC = sigma1./sqrt(m_MC);
    sigma_MFMC = sqrt(sigma1^2./m_MFMC(1) - sum((1./m_MFMC(1:end-1) - 1./m_MFMC(2:end)).*rho(2:end).^2.*sigma1^2));
end
