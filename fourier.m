function [D, x] = fourier(n)
% FOURIER  compute differentiation matrix and grid for Fourier modes
%
%   [D,X] = FOURIER(N)
%
%   Collocation points in X consist of N+1 points between 0 and 2pi
%
% From Trefethen, Spectral Methods in Matlab, p.22 (program 4)

if n == 0
    D = 0;
    x = 1;
    return
end

h = 2*pi/n;
x = h * (1:n)';
column = [0 0.5*(-1).^(1:n-1).*cot((1:n-1)*h/2)]';
row = column([1 n:-1:2]);
D = toeplitz(column, column([1 n:-1:2]));
