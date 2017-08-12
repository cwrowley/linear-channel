function sys = reduced_modelmod(A, B, C, T, Tinv, ntrunc,dt)
% MODIFIED TO ALLOW FOR PROJECTION ONTO POD MODES IN DIRECTION GIVEN BY ADJOINT 
% BPOD MODES (i.e., non-biorthogonal modes)
%
% modified to allow for output of discrete-time systems
T1 = T(:,1:ntrunc);
T1inv = Tinv(1:ntrunc,:);

M = T1inv*T1;

A1 = M*T1inv * A * T1;
B1 = M*T1inv * B;
C1 = C * T1;


    sys = ss(A1, B1, C1, 0);
if nargin == 7
    sys = ss(A1, B1, C1, 0,dt);
    % hack to ensure that discrete-time system has correct initial
    % condition for error analysis
    sys.b = B1;
end