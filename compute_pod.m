function [Tpod, PODsig, primal, tp] = compute_pod(sys, t)

if nargin > 1
    primal = impulse(sys, t);
    tp = t;
else
    [primal, tp] = impulse(sys);
end

%invert Laplacian for v part using Mmat
 
%primal = -Mmat*primal'; %why is this here? Converting back to Laplacian coords seemingly, perhaps to match with adjoint coords?
primal = primal'; % stay in velocity, vorticity coords

%primal = primal';
[Tpod, PODsig] = svd(primal,'econ');
PODsig = diag(PODsig);
