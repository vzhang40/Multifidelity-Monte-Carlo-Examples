% This function uses MC method to generate statistical values.
function stats = getStatsMC(models, N)
    % IN:
    %   models - functions ordered by correlation to hi-fidelity model
    %   N - number of samples used to calculate the stats
    % OUT:
    %   mfmc - a structure with fields: mu, sigma, and rho for means,
    %   standard deviations, and correlation coefficients
    f = @(x) [models(1).f(x), models(2).f(x)];
    z = sampleZ(N, 1);
    fvalues = f(z);
    stats.mus = mean(fvalues)';
    stats.sigma = std(fvalues)';
    rhos = corrcoef(fvalues);
    stats.rho = rhos(1, :)';
end

