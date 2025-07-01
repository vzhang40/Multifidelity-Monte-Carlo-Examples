% This function calculates the MC estimator R times
function s = doMC(f, m, R, flag)
    % IN: 
    %   f - truth model
    %   m - sample size
    %   R - number of replicates
    %   flag - "exp" for exponential example or "cdr" for CDR example
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