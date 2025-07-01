% This function uses the data to generate statistical values.
function stats = getStats(models)
    % IN:
    %   models - functions/models
    %   N - number of points used
    %   a - lower bound of domain
    %   b - upper bound of domain
    % OUT: 
    %   stats - a structure with fields:
    %       mus - mean values of models
    %       sigma - standard deviations for model values
    %       rho - correlation coefficient between surrogate and truth model
    stats.mus = mean(models)';
    stats.sigma = std(models)';
    rhos = corrcoef(models);
    stats.rho = rhos(1, :)';
end

