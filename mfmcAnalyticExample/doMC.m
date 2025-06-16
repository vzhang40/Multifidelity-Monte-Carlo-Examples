% This function calculates the MC estimator R times
function s = doMC(f, m, R)
    % IN: 
    %   f - model
    %   m - a vector of sample size
    %   R - number of replicates
    % OUT:
    %   s - a vector of size R, MC estimators
    z = sampleZ(m, R);
    s = mean(f(z));
end