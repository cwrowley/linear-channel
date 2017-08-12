function dist = subspace_test(S1,S2)
% Compute distance between subspaces using Grassmanian metric
dist = abs(acos(det(S1'*S2)));
end

