% This function does Multi-fidelity Monte Carlo Linear Regression for R replicates.
function [beta, Cxy] = doMFMCLR(models, mfmc, Cxx, d, R, n, flag, A)
    % IN: 
    %   models: functions in order of accuracy
    %   mfmc: structure with weights and sample sizes
    %   Cxx: exact value of Cxx matrix 
    %   d: degree of polynomial regression
    %   R: number of replicates
    %   n: index for cost
    %   flag: "exp" or "cdr" for the associated example
    %   A: if "exp", A = [], if "cdr" A is the input matrix
    % OUT:
    %   betaMFMC: a (# of features) x R matrix of R number of coefficients
    nFeat = size(Cxx, 1);
    if flag == "exp"
        fs = {models.f};
        col = length(fs);
        Xm = cell(1 ,col);
        Ym = cell(1, col);
        z = sampleZ(mfmc.m(n, end)*R, 1);
        Xm{col} = reshape((z).^(0:d)', [nFeat, mfmc.m(n, end), R]);
        z = reshape(z, [mfmc.m(n, end), R]);
        for j = 1:col
           Xm{j} = Xm{col}(:, 1:mfmc.m(n, j), :);
           Ym{j} = reshape(fs{j}(z(1:mfmc.m(n,j), :)), [mfmc.m(n,j), 1 ,R]);
        end
        
    elseif flag == "cdr"
        [rows, col] = size(models);
        Xm = cell(1 ,col);
        Ym = cell(1, col);
        z = randi(rows, [mfmc.m(n,end)*R, 1]);
        [temp, ~] = getX(A(z, :), d);
        Xm{col} = reshape(temp, [nFeat, mfmc.m(n, end), R]);
        z = reshape(z, [mfmc.m(n, end), R]);
        for j = 1:col
            Xm{j} = Xm{col}(:, 1:mfmc.m(n, j), :);
            Ym{j} = reshape(models(reshape(z(1:mfmc.m(n,j), :), [mfmc.m(n,j)*R, 1]), j), [mfmc.m(n,j), 1 ,R]);
        end
    end
    Cxy = (1./mfmc.m(n, 1)).*(pagemtimes(Xm{1},Ym{1}));
    for i = 2:col
        Cxy = Cxy + mfmc.alpha(i).*((1./mfmc.m(n, i)).*(pagemtimes(Xm{i},Ym{i})) - (1./mfmc.m(n, i-1)).*(pagemtimes(Xm{i-1},Ym{i}(1:mfmc.m(n, i-1), :, :))));
    end
    Cxy = reshape(Cxy, [nFeat, R]);
    beta = Cxx\Cxy;
end