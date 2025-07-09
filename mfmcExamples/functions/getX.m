% This function takes in the input matrix and the maximum power of the
% polynomial regression to create the X matrix (f = X'beta) and a power
% matrix.
function [X, powers] = getX(A,d)
    % IN: 
    %   A: a matrix with N rows and q columns representing N sets of input
    %   data with q inputs/parameters
    %   d: maximum polynomial order of the polynomial regression
    % OUT: 
    %   X: a matrix (N by # of input features) so that f = X'beta
    %   powers: a matrix (# of input features by q) that contains exponents
    %   for the input values.
    [N, q] = size(A);
    [temp{1:q}] = deal((0:d)');
    powers = combinations(temp{:}); % creates a table with all 5-digit combinations of 0, 1, 2
    powers = powers{:, :}; % turns table into matrix 
    mask = sum(powers, 2) <= 2; % check total degree order
    powers = powers(mask, :); % power matrix
    pp(1,:,:) = powers'; % adds a third dimension to the power matrix
    pp = repmat(pp,[N 1 1]); % repeats the power matrix to match the input matrix
    p = size(pp,3); % this is the number of features (terms) in our input matrix 
    temp = repmat(A,[1 1 p]); % repeating the input matrix to the number of different terms
    X = squeeze(prod(temp.^pp,2))';
end