function [T,sig,Tinv] = snapshotbal4mod(primal,dual,time,nmodes,nproj)
%%%
% MODIFIED TO 
%
% [T,sig,Tinv] = snapshotbal4(primal,dual,time,nmodes,nproj)
% 
% Computes balancing transformation T and Hankel singular values sig,
% given a collections of snapshots for primal and dual problems (stored in
% 'primal' and 'dual', resp).
%
% Also computes first (nmodes) columns of inverse, Tinv
%
% If 'nproj' is present, computes using Willcox methods, projecting
% each Gramian separately onto a subspace of dimension nproj.
%
% Snapshots are stored as columns of the matrices primal and dual, and
% the time vector is stored in tp (assumed the same for primal, dual data)

primal = dtweight(primal, time);
dual = dtweight(dual, time);

if nargin > 5
    [U,sig] = svd(primal);
    Xhalf = U(:,1:nproj) * sig(1:nproj,1:nproj);
    [U,sig] = svd(dual);
    Yhalf = U(:,1:nproj) * sig(1:nproj,1:nproj);
    [T, sig, Tinv] = halfgrambal2(Xhalf, Yhalf, nmodes);
else
    [T, sig, Tinv] = halfgrambal2(primal, dual, nmodes); %but this doesn't even include the orption for fewer inputs
end
