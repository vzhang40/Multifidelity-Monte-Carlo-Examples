% This function calculates the MC estimator R times
function s = doMC(f, m, R, flag)
    % IN: 
    %   f - model
    %   m - a vector of sample size
    %   R - number of replicates
    % OUT:
    %   s - a vector of size R, MC estimators
    if flag == "exp"
        z = sampleZ(m, R);
        s = mean(f(z));
    elseif flag == "cdr"
        z = randi(length(f), [m, R]);
        s = mean(reshape(f(z), [m, R]));
    end
end