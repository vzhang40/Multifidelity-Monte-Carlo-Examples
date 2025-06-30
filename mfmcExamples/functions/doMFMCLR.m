function [beta, Cxy] = doMFMCLR(models, mfmc, Cxx, d, R, n)
    % This function does MF Monte Carlo Linear Regression.
    % In: 
    % models: functions in order of accuracy
    % mfmc: structure with weights and sample sizes
    % Cxx: exact value of Cxx matrix 
    % d: degree of polynomial regression
    % R: number of replicates
    % n: index for cost
    %
    % Out:
    % betaMFMC: a (d+1) x R matrix of R number of coefficients 
    fs = {models.f};
    col = length(fs);
    Xm = cell(1 ,col);
    Ym = cell(1, col);
    z = sampleZ(mfmc.m(n, end)*R, 1);
    Xm{col} = reshape((z).^(0:d)', [d+1, mfmc.m(n, end), R]);
    z = reshape(z, [mfmc.m(n, end), R]);
    for j = 1:col
       Xm{j} = Xm{col}(:, 1:mfmc.m(n, j), :);
       Ym{j} = reshape(fs{j}(z(1:mfmc.m(n,j), :)), [mfmc.m(n,j), 1 ,R]);
    end
    Cxy = (1./mfmc.m(n, 1)).*(pagemtimes(Xm{1},Ym{1}));
    for i = 2:col
        Cxy = Cxy + mfmc.alpha(i).*((1./mfmc.m(n, i)).*(pagemtimes(Xm{i},Ym{i})) - (1./mfmc.m(n, i-1)).*(pagemtimes(Xm{i-1},Ym{i}(1:mfmc.m(n, i-1), :, :))));
    end
    Cxy = reshape(Cxy, [d+1, R]);
    beta = Cxx\Cxy;
end