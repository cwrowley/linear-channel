
nz = 8; %16
ny = 16; %32
Re = 1000;
% 1000, for longer, larger transient growths

dt = 0.05;
Tf = 1000;
tp = 0 : dt : Tf;

[A,Aadj, B, Mmat, npts, y, z, MmatTot] = define_eqns(ny, nz, Re);
C = eye(2 * npts);
channel = ss(A, B, C, 0);
[ysim, t] = impulse(channel);
tic
[T, sig, Tinv] = sysbalcwr(A, B, C);
t = toc;

fprintf('WORK HERE')
