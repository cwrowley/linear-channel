function [D,x] = cheb(n)
% CHEB  compute differentiation matrix and grid for Chebyshev polynomials
%   [D,X] = CHEB(N)
%
%   Collocation points in X consist of N+1 points between -1 and 1
%   From Trefethen, Spectral Methods in Matlab, p.54

if n == 0
    D = 0;
    x = 1;
    return
end

x = cos(pi * (0:n) / n)';
c = [2; ones(n-1, 1); 2] .* (-1).^(0:n)';
X = repmat(x, 1, n+1);
dX = X - X';
D = (c * (1./c)') ./ (dX + eye(n+1));
D = D - diag(sum(D'));
