function wsnapshots = dtweight(snapshots, time)
% wsnapshots = dtweight(snapshots, time)
%
% weight snapshots by quadrature coefficients, such that
%
% wsnapshots * wsnapshots' = \int_0^\infty (data)*(data)' dt

[n,ntsnap] = size(snapshots);
nt = length(time);
numruns = ntsnap/nt;

% weight snapshots by square roots of quadrature coefs (trapezoidal)
wsnapshots = snapshots;
for i=0:numruns-1
    wsnapshots(:, 1 + i*nt)  = snapshots(:, 1 + i*nt) * sqrt((time(2)-time(1))/2);
    wsnapshots(:, nt + i*nt) = snapshots(:, nt + i*nt) * sqrt((time(end)-time(end-1))/2);
    for j=2:nt-1
        wsnapshots(:, j + i*nt) = snapshots(:, j + i*nt) * sqrt((time(j+1)-time(j-1))/2);
    end
end