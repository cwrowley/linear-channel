function plotvar(z,y,u,clev)
% plotvar(z,y,u,clev)

if (nargin < 4), clev = 14; end

ny = length(y)-1;
nz = length(z);

% hack to duplicate endpoints of periodic BCs for symmetric plots
zz = zeros(nz+1,1);
zz(1) = 0;
zz(2:end) = z;

[Z,Y] = meshgrid(zz,y);
unew = zeros(ny+1,nz+1);

for j = 2:ny
    offset = nz*(j-2);
    unew(j,2:nz+1) = u(offset+1:offset+nz);
    unew(j,1) = unew(j,nz+1);
end

contourf(Z,Y,unew,clev)
axis equal
xlabel('z');
ylabel('y');