% This function uses MC method to generate statistical values.
function stats = getStats(models)
    % IN:
    %   models - functions ordered by correlation to hi-fidelity model
    % OUT:
    %   stats - a structure with fields: mu, sigma, and rho for mean,
    %   standard deviation, and correlation coefficient
    stats.mus = mean(models)';
    stats.sigma = std(models)';
    rhos = corrcoef(models);
    stats.rho = rhos(1, :)';
end

