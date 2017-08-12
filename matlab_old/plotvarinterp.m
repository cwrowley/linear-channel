function plotvarinterp(z,y,u,clev,Nzinterp,yinterp)
% add in functionality to interpolate across more data points, to make
% contours smoother. Note that interpolation should be done to match the
% appropriate spectral function (Fourier in z, Chebychev in y)

% original syntax
% plotvar(z,y,u,clev)

%if (nargin < 4), 
    clev = 14; 
%end

Ny = length(y)-1;
Nz = length(z);
Nyinterp = length(yinterp)-1;
%Nzinterp = length(zinterp);

% hack to duplicate endpoints of periodic BCs for symmetric plots
% remove for initial gridding, for proper interpolation
%zz = zeros(Nz+1,1);
%zz(1) = 0;
%zz(2:end) = z;

%[Z,Y] = meshgrid(zz,y);

unewold = zeros(Ny+1,Nz+1);
unewold(1,:) = 0;
unewold(Ny+1,:) = 0;
for i=1:Nz, for j=2:Ny
    unewold(j,i+1) = u(i + Nz * (j-2));
end, end
for j=2:Ny
    unewold(j,1) = u(Nz + Nz * (j-2));
end

% cut out last row, so that endpoints are not duplicated for interpolation
% purposes
unew = unewold(:,1:(end-1));

%unewinterpz = zeros(Ny+1,Nz+1);

unewinterpz = interpft(unew,Nzinterp,2);

%figure
%plot(unew(4,:),'ro-')
%hold on
%plot(unewold(4,:),'gx--')

%size(unewinterpz)
%unewinterpz(:,1) = unewinterpz(:,end);
%interpolation

% hack to duplicate endpoints of periodic BCs for symmetric plots
zzinterp = linspace(0,2*pi,Nzinterp+1);
unewinterpz = [unewinterpz(:,end),unewinterpz]; %add additional column for z=0;

[Z,Y] = meshgrid(zzinterp,y);


[c,h] = contourf(Z,Y,unewinterpz,clev); shading flat
%[c,h] = contour(Z,Y,unew,clev);
set(h, 'EdgeColor','none')
%set(h,'LineWidth',0)
set(gca,'FontSize',16)

axis equal
%xlabel('z');
%ylabel('y');