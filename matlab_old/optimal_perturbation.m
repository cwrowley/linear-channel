function vec = optimal_perturbation(ny, nz, A, t)
% OPTIMAL_PERTURBATION  Initial perturbation that maximizes energy growth
%
%   vec = OPTIMAL_PERTURBATION(ny, nz, A, t)

expA = expm(A * t);
n = length(A);
[~,~,V] = svd(expA);

% return a linear combination of the first two singular vectors so that
% vorticity = zero at z = 0 or 2*pi
vec1 = V(:,1);
vec2 = V(:,2);

% index corresponding to vorticity at a point with z = 2*pi
jmid = ny/2;
npts = nz * (ny-1);
ind = npts + nz * jmid;

a = vec1(ind);
b = vec2(ind);
r = a^2 + b^2;
a = a / r;
b = b / r;

vec = b * vec1 - a * vec2;

end