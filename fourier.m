% [D,x] = FOURIER(N)  compute D = differentiation matrix, x = Fourier grid
%
% (N points between 0 and 2pi)
%
% From Trefethen, Spectral Methods in Matlab, p.22 (program 4)

function [D,x] = fourier(N)
if N==0, D=0; x=1; return, end
h = 2*pi/N
x = h*(1:N)'
column = [0 0.5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]'
D = toeplitz(column, column([1 N:-1:2]))