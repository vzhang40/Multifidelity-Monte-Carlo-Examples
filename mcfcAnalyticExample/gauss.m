% This function evaluates the integral of f over [a,b] using Gaussian quadrature
function gaussInt = gauss(f, N, a, b)
    % IN:
    %   f - integrating function
    %   N - # of points
    %   a - lower bound of domain
    %   b - upper bound of domain
    % OUT:
    %   gaussInt - value of integral of f from a to b
    [x,w] = lgwt(N,a,b);
    summand = w.*f(x);
    gaussInt = sum(summand);
end
