function sys = reduced_model(A, B, C, T, Tinv, ntrunc,dt)
% modified to allow for output of discrete-time systems
T1 = T(:,1:ntrunc);
T1inv = Tinv(1:ntrunc,:);
A1 = T1inv * A * T1;
B1 = T1inv * B;
C1 = C * T1;


    sys = ss(A1, B1, C1, 0);
if nargin == 7
    sys = ss(A1, B1, C1, 0,dt);
    % hack to ensure that discrete-time system has correct initial
    % condition for error analysis
    sys.b = B1;
end