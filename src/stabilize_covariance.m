function P = stabilize_covariance(P, floor_value)
%STABILIZE_COVARIANCE Symmetrize a covariance matrix and floor eigenvalues.

P = real((P + P') / 2);
[V, D] = eig(P);
eigenvalues = max(real(diag(D)), floor_value);
P = real(V * diag(eigenvalues) * V');
P = real((P + P') / 2);
end
