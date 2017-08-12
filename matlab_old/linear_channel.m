function [A, B, M, npts, y, z] = linear_channel(ny, nz, Re)
% LINEAR_CHANNEL  Compute matrices for streamwise constant linear channel
%   [A, B, M, npts, y, z] = LINEAR_CHANNEL(ny, nz, Re)
%
%   ny: number of Chebyshev modes in y-dir
%   nz: number of Fourier modes in z-dir
%   Re: Reynolds number
%
%   Equations are xdot = Ax + Bu
%   npts: number of grid points ((ny-1) * nz; state dimension is 2 * npts)
%   y:    grid in wall-normal direction (ny+1 points)
%   z:    grid in spanwise direction (nz points)

npts = (ny-1) * nz;
testderivs = 0;  % flag: set > 0 to run derivative tests

% Chebyshev modes in y-dir
[Dy11,y] = cheb(ny);
Dyy1 = Dy11^2;
Dy1 = Dy11(2:ny,2:ny);
Dyy1 = Dyy1(2:ny,2:ny); % homogeneous boundary conditions
%
% fourth-order boundary conditions
S = diag([0; 1./(1-y(2:ny).^2); 0]);
D4y = (diag(1-y.^2)*Dy11^4 - 8*diag(y)*Dy11^3 - 12*Dy11^2)*S;
D4y1=D4y(2:ny,2:ny);
%clamped BC's on v
D2y = (diag(1-y.^2)*Dy11^2 - 4*diag(y)*Dy11 - 2*eye(ny+1))*S;
Dyy1 = D2y(2:ny,2:ny);



% stack into matrix
Dy = zeros(npts);
Dyy = zeros(npts);
for i=1:ny-1, for j=1:ny-1
        Dy((i-1)*nz+1:i*nz, (j-1)*nz+1:j*nz) = diag(Dy1(i,j)*ones(nz,1));
        Dyy((i-1)*nz+1:i*nz, (j-1)*nz+1:j*nz) = diag(Dyy1(i,j)*ones(nz,1));
end, end

% stack into matrix
D4y = zeros(npts);
%Dyy = zeros(Npts);
for i=1:ny-1, for j=1:ny-1
        D4y((i-1)*nz+1:i*nz, (j-1)*nz+1:j*nz) = diag(D4y1(i,j)*ones(nz,1));
%        Dyy((i-1)*Nz+1:i*Nz, (j-1)*Nz+1:j*Nz) = diag(Dyy1(i,j)*ones(Nz,1));
end, end


% Fourier modes in z-dir
[Dz1,z] = fourier(nz);
Dzz1 = Dz1^2;

% stack into matrix
Dz = zeros(npts);
Dzz = zeros(npts);
for i=1:ny-1
    Dz((i-1)*nz+1:i*nz, (i-1)*nz+1:i*nz) = Dz1;
    Dzz((i-1)*nz+1:i*nz, (i-1)*nz+1:i*nz) = Dzz1;
end

% Define grid
zvar = zeros(npts,1); yvar = zeros(npts,1);
for j=1:ny-1
    yvar((j-1)*nz+1:j*nz) = y(j+1);
    zvar((j-1)*nz+1:j*nz) = z;
end

% Mean flow
Uvar = 1-yvar.^2;
Uprime = -2*yvar;

if (testderivs)
    uvar = zeros(npts,1);
    uvar_y = zeros(npts,1);
    uvar_z = zeros(npts,1);
    uvar_yy = zeros(npts,1);
    uvar_zz = zeros(npts,1);
    for i=1:nz, for j=1:ny-1
            zz = z(i); yy = y(j+1);
            uvar(i + (j-1)*nz) = (1-yy^2) * sin(zz);
            uvar_y(i + (j-1)*nz) = -2*yy * sin(zz);
            uvar_z(i + (j-1)*nz) = (1-yy^2) * cos(zz);
            uvar_yy(i + (j-1)*nz) = -2 * sin(zz);
            uvar_zz(i + (j-1)*nz) = -(1-yy^2) * sin(zz);
        end, end
    
    erry = Dy*uvar - uvar_y;
    errz = Dz*uvar - uvar_z;
    erryy = Dyy*uvar - uvar_yy;
    errzz = Dzz*uvar - uvar_zz;
    
    figure(1)
    plotvar(z,y,uvar_y)
    figure(2)
    plotvar(z,y,Dy*uvar)
    
    figure(3)
    ind = 1:npts;
    plot(ind,erry,ind,errz,ind,erryy,ind,errzz)
    legend('y','z','yy','zz')
end


M = [-(Dyy+Dzz), zeros(npts); zeros(npts), eye(npts)];
% the A matrix for streamwise-constant perturbations
A = inv(M) * [-1/Re* (Dzz^2 + 2*Dyy*Dzz + D4y), zeros(npts);
       -diag(Uprime) * Dz, 1/Re * (Dyy + Dzz)];

M = (M + M')/2;

%Aadj = inv(Mmat)*A'*Mmat;
%Aadj=A'; %changed!

%B = write_ic(ny,nz,z,y,1);
%ic_width = 1;
%B = velocity_perturbation(y, z, ic_width);
T_final = 32.9;
B = optimal_perturbation(ny, nz, A, T_final);
