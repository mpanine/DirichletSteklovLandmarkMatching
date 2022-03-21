% Commute time distances

load('C:\Users\nastyad\Technion\Tools\3D\armadillos\al0.mat');
load('C:\Users\nastyad\Technion\Tools\3D\armadillos\al3.mat');

[W, A] = mshlp_matrix(shape);

% Eigendecomposition
Am = sparse([1:length(A)], [1:length(A)], A);
[evecs evals] = eigs(W, Am, 200, -1e-5);
evals = diag(evals);
% Remove zero eigenvalue and corresponding eigenvector
evals = evals(2:end);
evecs = evecs(:, 2:end);

% Commute time distance
LL = 1./abs(evals);
K  = evecs*diag(LL(:))*evecs';
D  = sqrt(repmat(diag(K),[1 size(K,2)]) + repmat(diag(K)',[size(K,2) 1]) - 2*K);
% D(D < eps) = eps; N = size(D, 1);
% D(1:(N + 1):N^2) = 0;

% Geodesic distance - for comparison
Dg = fastmarch(shape);

% MDS
addpath('C:\Users\nastyad\Technion\MSc\Prev courses - relevant\Numerical Geometry of Images\HWs\hw3');
% Y = classical_mds(D, 3);
Y = classical_mds(Dg, 3);
shape2 = struct('TRIV', shape.TRIV, 'X', Y(:, 1), 'Y', Y(:, 2), 'Z', Y(:, 3));
figure; draw_shape(shape2);



count = 1;
for t = [256 512 1024 4096],
    LL = 1./abs(evals);
    K  = evecs*diag(LL(:))*evecs';
    D  = sqrt(repmat(diag(K),[1 size(K,2)]) + repmat(diag(K)',[size(K,2) 1]) - 2*K);
    subplot(1,4,count); count = count+1;
    trisurf(shape.TRIV, shape.X, shape.Y, shape.Z, D(1,:)); axis image; axis off; shading interp;
end
