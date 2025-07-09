% This function gives the exact value for Cxx as well as Cxy and the
% best fit polynomial with order d over the interval [a,b]
function exactLR = getExactLR(f, d, flag, a, b)
    % IN: 
    %   f: the truth model
    %   d: maximum polynomial degree for regression
    %   flag: "exp" if exponential example, "cdr" if CDR example
    %   if flag == "exp"
    %   [a, b]: intergration bounds; domain of interest
    %   if flag == "cdr"
    %   a: matrix of input values
    %   b: []
    % OUT:
    %   exactLR: structure with fields:
    %       Cxx: exact value of E[XX']
    %       Cxy: "exact" value of E[XY]
    %       beta: optimal beta chocie to approximate f
    %       poly: anonymous function of best fit polynomial (only for
    %       "exp")
    if flag == "exp"
        exactLR.Cxx = zeros(d+1);
        N_gauss = 2*d + 2;
        [x,w] = lgwt(N_gauss,a,b);
        exactLR.Cxy = (1./(b-a)).*x.^(0:d)'*(f(x).*w); % values of Cxy are average values of x^(i-1)y
        vecX = (1./(b-a))*(x.^(0:d*2)'*w);
        for i = 1:d+1
            exactLR.Cxx(:, i) = vecX(i:i+d);
        end
    elseif flag == "cdr"
        n = size(a, 1);
        [X, ~] = getX(a,d);
        exactLR.Cxx = (1/n)*(X)*(X');
        exactLR.Cxy = (1/n)*X*f;
    end
    exactLR.beta = exactLR.Cxx\exactLR.Cxy;
    if flag == "exp"
        exactLR.poly = @(x) x.^(0:d)*exactLR.beta;
    end
end

