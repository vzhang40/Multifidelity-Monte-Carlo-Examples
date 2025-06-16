% This function calculates the MFMC estimator R times
function s = doMFMC(models, mfmc, R, n)
    % IN: 
    %   fs - must be a cell array of models
    %   alpha - a vector of weights
    %   m - a vector of sample sizes
    %   R - number of replicates
    % OUT:
    %   s - a vector of size R, MFMC estimators
    z = sampleZ(mfmc.m(n, end), R);
    fs = {models.f};
    data = cell(1, length(fs));
    for j = 1:length(fs)
        data{j} = fs{j}(z(1:mfmc.m(n, j), :));
    end
    s = mean(data{1});
    for i = 2:length(fs)
        s = s + mfmc.alpha(i).*(mean(data{i}) - mean(data{i}(1:mfmc.m(n, i-1), :)));
    end
end