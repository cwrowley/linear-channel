function [T, sig, Tinv] = halfgrambal2(Xhalf, Yhalf, nmodes)
% [T, sig, Tinv] = halfgrambal2(Xhalf, Yhalf, nmodes)
% 
% Computes balancing transformation T and Hankel singular values sig,
% given a factorization of the gramians X and Y such that
%   Wc = Xhalf * Xhalf'
%   Wo = Yhalf * Yhalf'
%
% Also computes first (nmodes) columns of inverse, Tinv
%
% Reduces rank of Gramians before computing svd, if #snapshots > 2*n

% take inner products of snapshots
[n,nsnaps] = size(Xhalf);
if nsnaps > 2*n
    disp('Computing singular values of Xhalf')
    [U, sig] = svd(Xhalf);
    % r = rank(sig);
    r = size(find(diag(sig)), 1);
    Xhalf = U(:,1:r) * sig(1:r,1:r);
    disp(sprintf('    rank = %d',r))
end

[n, nsnaps] = size(Yhalf);
if nsnaps > 2*n
    disp('Computing singular values of Yhalf')
    [U, sig] = svd(Yhalf);
    % r = rank(sig);
    r = size(find(diag(sig)), 1);
    Yhalf = U(:,1:r) * sig(1:r,1:r);
    disp(sprintf('    rank = %d',r))
end


%for i = 1:size(Yhalf,2) ; for j = 1:size(Xhalf,2)
%   YMX(i,j) = mip(Yhalf(:,i),Xhalf(:,i),Ny,Nz)
%end; end

%[u,sig,v]= svd(YMX);
[u,sig,v] = svd(Yhalf' * Xhalf);
sig = diag(sig);
r = size(find(sig),1);  % number of nonzero singular values
if nargin > 2
    r = min([nmodes r]);
end
sig = sig(1:r);
u = u(:,1:r);
v = v(:,1:r);
sighalfinv = diag(sqrt(1./sig));
disp('Computing balancing transformation')
T = Xhalf * v * sighalfinv;
Tinv = sighalfinv * u' * Yhalf';
