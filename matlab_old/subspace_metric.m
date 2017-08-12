function dist = subspace_metric(S1,S2)
% Compute distance between subspaces using Grassmanian metric
dist = sqrt(1 - abs(det(S1'*S2)));
end

