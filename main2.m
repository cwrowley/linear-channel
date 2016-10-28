[A, B, npts, y, z] = define_eqns(ny, nz, Re)
C = eye(2 * npts);
channel = ss(A, B, C, 0);
[y, t] = impulse(channel
tic
[T, sig, Tinv] = sysbalcwr(A, B, C);
t = toc;
fprintf('WORK HERE')
