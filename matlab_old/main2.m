ny = 20;
nz = 20;
Re = 2000;
[A, B, M, npts, y, z] = linear_channel(ny, nz, Re);
C = eye(2 * npts);
channel = ss(A, B, C, 0);
[y, t] = impulse(channel);
