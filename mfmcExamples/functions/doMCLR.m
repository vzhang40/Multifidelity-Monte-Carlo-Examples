% This function does Monte Carlo Linear Regression for R replicates.
function [beta, Cxy] = doMCLR(f, Cxx, d, m, R, flag, A)
    % IN: 
    %   f: function that we want to approximate
    %   Cxx: exact value of Cxx matrix 
    %   d: degree of polynomial regression
    %   m: Monte Carlo Sample Size
    %   R: number of replicates
    %   flag: "exp" or "cdr" for the associated example
    %   A: if "exp", A = [], if "cdr" A is the input matrix
    % OUT:
    %   betaMC: a (# of features) x R matrix of R number of coefficients 
    nFeat = size(Cxx, 1);
    if flag == "exp"
        randX = sampleZ(m*R, 1);
        Xm = reshape((randX).^(0:d)', [nFeat, m, R]);
        Ym = reshape(f(randX), [m, 1, R]);
    elseif flag == "cdr"
        randX = randi(length(f), [m*R, 1]);
        [Xm, ~] = getX(A(randX, :), d);
        Xm = reshape(Xm, [nFeat, m, R]);
        Ym = reshape(f(randX, :), [m, 1, R]);
    end
    Cxy = reshape((1./m)*pagemtimes(Xm,Ym), [nFeat, R]);
    beta = Cxx\Cxy;
end
