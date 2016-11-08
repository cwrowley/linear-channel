function [A,Aadj, B, Mmat, Npts, y, z, MmatTot] = define_eqns(Ny, Nz, Re)
% define_eqns - compute matrices for linear channel (streamwise const)
%
%  [A, Aadj, B, Mmat, y, z] = define_eqns(Ny, Nz, Re)
%
%  Ny: number of Chebyshev modes in y-dir
%  Nz: number of Fourier modes in z-dir
%  Re: Reynolds number
%
%  Equations are xdot = Ax + Bu
%  Npts: number of gridpoints in channel (dim(A) = 2*Npts)
%  y:    grid in wall-normal direction
%  z:    grid in spanwise direction
%
% modified by Milos Ilak
% returns forward and adjoint A matrix, M matrix

Npts = (Ny-1)*Nz;
testderivs = 0;  % flag: set > 0 to run derivative tests

% Chebyshev modes in y-dir
[Dy11,y] = cheb(Ny);
Dyy1 = Dy11^2;
Dy1 = Dy11(2:Ny,2:Ny);
Dyy1 = Dyy1(2:Ny,2:Ny); % homogeneous boundary conditions
%
% fourth-order boundary conditions
  S = diag([0; 1./(1-y(2:Ny).^2); 0]);
 D4y = (diag(1-y.^2)*Dy11^4 - 8*diag(y)*Dy11^3 - 12*Dy11^2)*S;
D4y1=D4y(2:Ny,2:Ny);
%clamped BC's on v
D2y = (diag(1-y.^2)*Dy11^2 - 4*diag(y)*Dy11 - 2*eye(Ny+1))*S;
  Dyy1 = D2y(2:Ny,2:Ny);



% stack into matrix
Dy = zeros(Npts);
Dyy = zeros(Npts);
for i=1:Ny-1, for j=1:Ny-1
        Dy((i-1)*Nz+1:i*Nz, (j-1)*Nz+1:j*Nz) = diag(Dy1(i,j)*ones(Nz,1));
        Dyy((i-1)*Nz+1:i*Nz, (j-1)*Nz+1:j*Nz) = diag(Dyy1(i,j)*ones(Nz,1));
end, end

% stack into matrix
D4y = zeros(Npts);
%Dyy = zeros(Npts);
for i=1:Ny-1, for j=1:Ny-1
        D4y((i-1)*Nz+1:i*Nz, (j-1)*Nz+1:j*Nz) = diag(D4y1(i,j)*ones(Nz,1));
%        Dyy((i-1)*Nz+1:i*Nz, (j-1)*Nz+1:j*Nz) = diag(Dyy1(i,j)*ones(Nz,1));
end, end


% Fourier modes in z-dir
[Dz1,z] = fourier(Nz);
Dzz1 = Dz1^2;

% stack into matrix
Dz = zeros(Npts);
Dzz = zeros(Npts);
for i=1:Ny-1
    Dz((i-1)*Nz+1:i*Nz, (i-1)*Nz+1:i*Nz) = Dz1;
    Dzz((i-1)*Nz+1:i*Nz, (i-1)*Nz+1:i*Nz) = Dzz1;
end

% Define grid
zvar = zeros(Npts,1); yvar = zeros(Npts,1);
for j=1:Ny-1
    yvar((j-1)*Nz+1:j*Nz) = y(j+1);
    zvar((j-1)*Nz+1:j*Nz) = z;
end

% Mean flow
Uvar = 1-yvar.^2;
Uprime = -2*yvar;

if (testderivs)
    uvar = zeros(Npts,1);
    uvar_y = zeros(Npts,1);
    uvar_z = zeros(Npts,1);
    uvar_yy = zeros(Npts,1);
    uvar_zz = zeros(Npts,1);
    for i=1:Nz, for j=1:Ny-1
            zz = z(i); yy = y(j+1);
            uvar(i + (j-1)*Nz) = (1-yy^2) * sin(zz);
            uvar_y(i + (j-1)*Nz) = -2*yy * sin(zz);
            uvar_z(i + (j-1)*Nz) = (1-yy^2) * cos(zz);
            uvar_yy(i + (j-1)*Nz) = -2 * sin(zz);
            uvar_zz(i + (j-1)*Nz) = -(1-yy^2) * sin(zz);
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
    ind = 1:Npts;
    plot(ind,erry,ind,errz,ind,erryy,ind,errzz)
    legend('y','z','yy','zz')
end

 % A = [1/Re * (Dyy + Dzz), zeros(Npts);
 %     -diag(Uprime) * Dz, 1/Re * (Dyy + Dzz)];
%  



Mmat = [-(Dyy+Dzz), zeros(Npts); zeros(Npts), eye(Npts)];
%the A matrix for streamwise-constant data
A = [-1/Re* (Dzz^2 + 2*Dyy*Dzz + D4y), zeros(Npts);
       -diag(Uprime) * Dz, 1/Re * (Dyy + Dzz)];
A=inv(Mmat)*A;

%figure(1)
%e = eig(A);
%plot(real(e),imag(e),'bx')
%title('Eigenvalues of A')

Mmat = (Mmat + Mmat')/2;
MmatHalf = Mmat^0.5;
Mint = defineIntegralWeights(y,z);
MmatTot = MmatHalf*diag([Mint;Mint])*MmatHalf;

Aadj = MmatTot^-1*A'*MmatTot;


%adjoint equations (streamwise constant, so only off-diag. term moves)
% Aadj = [1/Re* (Dzz^2 + 2*Dyy*Dzz + D4y), -diag(Uprime)*Dz;
%    zeros(Npts), 1/Re * (Dyy + Dzz)];
% Aadj=inv(Mmat)*Aadj;
 %Aadj = Aadj*inv(Mmat);
 %A=A';


%Aadj = inv(Mmat)*A'*Mmat;
Aadj=A'; %changed!


B=write_ic(Ny,Nz,z,y,1);




%B(1:2*Npts)=1;

%B(Npts + i + (j-1)*Nz) = 1;   % disturbance in eta
