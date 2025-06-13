% This function uses gauss quadrature to generate statisical values
function [mu_true, sigma, rho] = getStatsGauss(f1, f, N, a, b)
    % IN:
    %   f - functions ordered by correlation to hi-fidelity model
    %   N - number of points used
    %   a - lower bound of domain
    %   b - upper bound of domain
    % OUT: 
    %   mu_true - the average value of the truth function
    %   sigma - the standard deviations of the models
    %   rho - the correlation coefficients of the models compared to truth
    mus = (1./(b-a)).*gauss(f, N, a, b); 
    mu_true = mus(1);
    sigma = sqrt((1./(b-a)).*gauss(@(x) (f(x) - mus).^2, N, a, b))';
    rho = (1./(sigma(1).*sigma)).*(1./(b-a)).*gauss(@(x) (f1(x) - mu_true).*(f(x) - mus), N, a, b)';
end
