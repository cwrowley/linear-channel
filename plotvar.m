function plotvar(z,y,u,clev)
% plotvar(z,y,u,clev)

if (nargin < 4), clev = 14; end

Ny = length(y)-1;
Nz = length(z);

% hack to duplicate endpoints of periodic BCs for symmetric plots
zz = zeros(Nz+1,1);
zz(1) = 0;
zz(2:end) = z;

[Z,Y] = meshgrid(zz,y);
unew = zeros(Ny+1,Nz+1);
unew(1,:) = 0;
unew(Ny+1,:) = 0;
for i=1:Nz, for j=2:Ny
    unew(j,i+1) = u(i + Nz * (j-2));
end, end
for j=2:Ny
    unew(j,1) = u(Nz + Nz * (j-2));
end

contourf(Z,Y,unew,clev)
%[c,h] = contour(Z,Y,unew,clev);
%set(h,'LineWidth',1)
set(gca,'FontSize',16)
axis equal
%xlabel('z');
%ylabel('y');