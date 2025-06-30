% This functions draws NxR samples from a uniform distribution over [0,5]
function z = sampleZ(N, R)
    % IN:
    %   N - row dimension
    %   R - column dimension
    % OUT:
    %   z - a NxR matrix with random samples
    z = 5.*rand([N,R]);
end

