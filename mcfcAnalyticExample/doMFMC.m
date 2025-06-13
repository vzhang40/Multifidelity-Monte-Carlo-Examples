% This function calculates the MFMC estimator R times
function s = doMFMC(fs, alpha, m, R)
    % IN: 
    %   fs - must be a cell array of models
    %   alpha - a vector of weights
    %   m - a vector of sample sizes
    %   R - number of replicates
    % OUT:
    %   s - a vector of size R, MFMC estimators
    z = sampleZ(m(end), R);
    data = cell(1, 2);
    for j = 1:length(fs)
        data{j} = fs{j}(z(1:m(j), :));
    end
    s = mean(data{1});
    for i = 2:length(fs)
        s = s + alpha(i).*(mean(data{2}) - mean(data{2}(1:m(i-1), :)));
    end
end