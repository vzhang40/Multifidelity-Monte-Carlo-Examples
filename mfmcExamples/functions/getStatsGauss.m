% This function uses gauss quadrature to generate statisical values
function stats = getStatsGauss(models, N, a, b)
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
    fs = {models.f};
    for i = 1:length(fs)
        stats.mus(i) = (1./(b-a)).*gauss(fs{i}, N, a, b); 
        stats.sigma(i) = sqrt((1./(b-a)).*gauss(@(x) (fs{i}(x) - stats.mus(i)).^2, N, a, b))';
        stats.rho(i) = (1./(stats.sigma(1).*stats.sigma(i))).*(1./(b-a)).*gauss(@(x) (fs{1}(x) - stats.mus(1)).*(fs{i}(x) - stats.mus(i)), N, a, b)';
    end
    stats.mus = stats.mus';
    stats.sigma = stats.sigma';
    stats.rho = stats.rho';
end
