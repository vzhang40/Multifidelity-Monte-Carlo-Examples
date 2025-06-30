% This function calculates the MFMC estimator R times
function s = doMFMC(models, mfmc, R, n, flag)
    % IN: 
    %   fs - must be a cell array of models
    %   alpha - a vector of weights
    %   m - a vector of sample sizes
    %   R - number of replicates
    % OUT:
    %   s - a vector of size R, MFMC estimators
    if flag == "exp"
        z = sampleZ(mfmc.m(n, end), R);
        fs = {models.f};
        col = length(fs);
        data = cell(1, col);
        for j = 1:col
            data{j} = fs{j}(z(1:mfmc.m(n, j), :));
        end
    elseif flag == "cdr"
        [rows, col] = size(models.f);
        z = randi(rows, mfmc.m(n,end), R);
        data = cell(1, col);
        for j = 1:col
            data{j} = reshape(models.f(z(1:mfmc.m(n, j), :), j), [mfmc.m(n,j), R]);
        end
    end
    s = mean(data{1});
    for i = 2:col
        s = s + mfmc.alpha(i).*(mean(data{i}) - mean(data{i}(1:mfmc.m(n, i-1), :)));
    end
end