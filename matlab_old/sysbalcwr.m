% function [T, sig, Tinv] = sysbalcwr(A,B,C)
%
%   Finds a balancing transformation for the state-space model.
%   Eigenvalues of A must have negative real part. 
%
%  From MUTOOLS toolbox -- modified from sysbal by CWR, 14 Aug 2003

%   Copyright 1991-2002 MUSYN Inc. and The MathWorks, Inc. 
% $Revision: 1.7 $

function [T, sig, Tinv,A,B,C] = sysbalcwr(A,B,C)
   if nargin == 0
     disp(['usage: [T, sig, Tinv] = sysbalcwr(A,B,C)']);
     return
   end %if nargin<1

   [n,m]=size(B); [p,n]=size(C);
   [Ts,A]=schur(A);
   B = Ts'*B;
   C = C*Ts;
 % check that A is stable.
   if any(real(eig(A))>=0),
     disp('SYS must be stable')
     return
    end

 % find observability Gramian, S'*S (S upper triangular)
   S = sjh6(A,C);
 % find controllability Gramian R*R' (R upper triangular)
   perm = n:-1:1;
   R = sjh6(A(perm,perm)',B(perm,:)');
   R = R(perm,perm)';
 % calculate the Hankel-singular values
   [U,T,V] = svd(S*R);
   sig = diag(T);
   sighalfinv = 1./sqrt(sig);
   % balancing coordinates
   Tinv = diag(sighalfinv)*U'*S;
   B = Tinv*B; A = Tinv*A;
   T = R*V*diag(sighalfinv);
   C = C*T; A = A*T;

   T = Ts * T;
   Tinv = Tinv * Ts';
