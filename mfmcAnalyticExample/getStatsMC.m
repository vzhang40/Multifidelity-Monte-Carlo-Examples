% This function uses MC method to generate statistical values.
function [mu_true, sigma, rho] = getStatsMC(f, N)
    % IN:
    %   f - functions ordered by correlation to hi-fidelity model
    %   N - number of samples used to calculate the stats
    % OUT:
    %   mu_true - the mean of the truth function f1
    %   sigma - standard deviations of the function
    %   rho - correlation coefficients by model; rho_{1, i}
    z = sampleZ(N, 1);
    fvalues = f(z);
    mu_true = mean(fvalues(:, 1));
    sigma = std(fvalues)';
    rhos = corrcoef(fvalues);
    rho = rhos(1, :)';
end

