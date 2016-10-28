% [D,x] = CHEB(N)  compute D = differentiation matrix, x = Chebyshev grid
%
% (N+1 points between -1 and 1)
% From Trefethen, Spectral Methods in Matlab, p.54

function [D,x] = cheb(N)
if N==0, D=0; x=1; return, end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+eye(N+1));
D = D - diag(sum(D'));
