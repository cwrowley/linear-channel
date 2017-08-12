function [init_vec] = velocity_perturbation(y, z, icwidth)

nz = length(z);
ny = length(y) - 1;
npts = nz * (ny-1);

if nargin < 3
    icwidth = 1;
end

% vertical velocity perturbation
vel = zeros(ny-1, nz);
y_mid = 0;
z_mid = z(nz/2);
for k = 1:(ny-1)
    for m = 1:nz
        vel(k,m) = 0.5 * (cos(pi*y(k+1))+1) * ...
            exp(-2 * ((y(k+1)-y_mid)^2 + (z(m)-z_mid)^2) / icwidth^2);
    end
end

% re-write as vector
vel_vec = zeros(npts,1);
for k=1:(ny-1)
    vel_vec((k-1)*nz+1:k*nz) = vel(k,:);
end

% vorticity perturbation is zero
init_vec = [vel_vec; zeros(npts,1)];
 
