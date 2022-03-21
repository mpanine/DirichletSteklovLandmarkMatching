function [Basis, evals] = computeDirichletLB(Shape, boundary_indices, numEigs)
%computeDirichletLB computes the Laplace-Beltrami basis subject to
%Dirichlet boundary conditions on the boundary indices.
import basis.harmonic.harmonic_utils.*

[Wii, ~, Aii] = boundary_poisson_split(Shape, boundary_indices); % The boundary term Wib shouldn't be necessary with Dirichlet BC.

boundary_values = zeros(length(boundary_indices),1);
    
try
    [evecs, evals] = eigs(Wii, Aii, numEigs, 1e-6);
catch
    % In case of trouble make the laplacian definite
    [evecs, evals] = eigs(Wii - 1e-8*speye(size(Wii,1)), Aii, numEigs, 'sm');
end


evals = diag(evals);

[evals, order] = sort(abs(evals),'ascend');


evecs = evecs(:,order);

Basis = zeros(Shape.nv, numEigs);

for i = 1:numEigs
    Basis(:,i) = BoundaryReinsert(evecs(:,i), boundary_indices, boundary_values );
end






end