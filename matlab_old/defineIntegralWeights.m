function M = defineIntegralWeights(y,z);

Ny = length(y)-1;
Nz = length(z);

[~,dy] = clencurt(Ny);
dy(1) = []; % don't need these, since the function is zero there anyway (careful though, vorticity is nonzero)
dy(end) = [];

dz = gradient(z); % will be constant, so just use this

for i=1:Nz, for j=1:(Ny-1)
        M(i+Nz*(j-1)) = abs(dz(i)*dy(j));
end, end
M = M';