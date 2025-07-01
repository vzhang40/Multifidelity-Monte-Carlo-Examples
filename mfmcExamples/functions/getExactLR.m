% This function gives the exact value for Cxx as well as Cxy and the
% best fit polynomial with order d over the interval [a,b]
function exactLR = getExactLR(f, d, a, b)
    % IN: 
    %   f: the truth model
    %   d: maximum polynomial degree for regression
    %   [a, b]: intergration bounds; domain of interest
    %
    % OUT:
    %   exactLR: structure with fields:
    %       Cxx: exact value of E[XX']
    %       Cxy: "exact" value of E[XY]
    %       beta: optimal beta chocie to approximate f
    %       poly: anonymous function of best fit polynomial
    
    exactLR.Cxx = zeros(d+1);
    N_gauss = 2*d + 2;
    [x,w] = lgwt(N_gauss,a,b);
    exactLR.Cxy = (1./(b-a)).*x.^(0:d)'*(f(x).*w); % values of Cxy are average values of x^(i-1)y
    vecX = (1./(b-a))*(x.^(0:d*2)'*w);
    for i = 1:d+1
        exactLR.Cxx(:, i) = vecX(i:i+d);
    end
    exactLR.beta = exactLR.Cxx\exactLR.Cxy;
    exactLR.poly = @(x) x.^(0:d)*exactLR.beta;
end