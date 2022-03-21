function [Basis, evals] = computeDirichletLB(Shape, boundary_indices, numEigs, numEigs_margin,eig_threshold)
%computeDirichletLB computes the Laplace-Beltrami basis subject to
%Dirichlet boundary conditions on the boundary indices.

[Wii, ~, Aii] = boundary_poisson_split(Shape, boundary_indices); % The boundary term Wib shouldn't be necessary with Dirichlet BC.

boundary_values = zeros(length(boundary_indices),1);
    
try
    [evecs, evals] = eigs(Wii, Aii, numEigs+numEigs_margin, 1e-6);
catch
    % In case of trouble make the laplacian definite
    [evecs, evals] = eigs(Wii - 1e-8*speye(size(Wii,1)), Aii, numEigs+numEigs_margin, 'sm');
end

prods = evecs'*Aii*evecs;
inv_pr = diag( diag(prods).^-0.5 ); %L2 Normalization, just in case (sometimes necessary)
evecs = evecs * inv_pr;

evals = diag(evals);

[evals, order] = sort(abs(evals),'ascend');


evecs = evecs(:,order);

to_remove = evals<eig_threshold;
num_removed = sum(to_remove);
assert(num_removed<=numEigs_margin);
if num_removed>0
    evals(to_remove) = [];
    evecs(:,to_remove) = [];
end

sign_of_first_ef = sign( sum( sign( evecs(:,1) ) ) ); % Fix a positive sign for the first eigenfunction.
evecs(:,1) = sign_of_first_ef*evecs(:,1);

Basis = zeros(Shape.nv, numEigs);

for i = 1:numEigs
    Basis(:,i) = BoundaryReinsert(evecs(:,i), boundary_indices, boundary_values );
end
evals = evals(1:numEigs);





end