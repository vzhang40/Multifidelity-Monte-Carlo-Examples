% This function does Monte Carlo Linear Regression for R replicates.
function [beta, Cxy] = doMCLR(f, Cxx, d, m, R)
    % IN: 
    %   f: function that we want to approximate
    %   Cxx: exact value of Cxx matrix 
    %   d: degree of polynomial regression
    %   m: Monte Carlo Sample Size
    %   R: number of replicates
    % 
    % OUT:
    %   betaMC: a (d+1) x R matrix of R number of coefficients 

    randX = sampleZ(m*R, 1);
    Xm = reshape((randX).^(0:d)', [d+1, m, R]);
    Ym = reshape(f(randX), [m, 1, R]);
    Cxy = reshape((1./m)*pagemtimes(Xm,Ym), [d+1, R]);
    beta = Cxx\Cxy;
end
