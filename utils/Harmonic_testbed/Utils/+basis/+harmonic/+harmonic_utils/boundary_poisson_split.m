function [Wii, Wib, Aii] = boundary_poisson_split(Shape, boundary_indices)
% Splits the cotangent Laplacian on a mesh to interior-interior and
% boundary-interior parts in order to solve the Poisson equation.
% Also produces the interior lumped mass matrix.

%Shape: output of MeshInfo, contains the cotangent Laplacian
%boundary_indices: indices of the boundary elements

%Wii: interior-interior cotangent Laplacian
%Wib: interior-boundary cotangent Laplacian
%Aii: interior-interior lumped mass matrix

interior_vertices = 1:size(Shape.W,1);

int_ind = ~ismember(interior_vertices, boundary_indices) ;%True for all non-boundary vertices.
interior_vertices = interior_vertices(int_ind);

Wii = Shape.W(interior_vertices,interior_vertices);
Wib = Shape.W(interior_vertices, boundary_indices);

% A = sparse(diag( sum(Shape.A) ));

Aii = Shape.A(interior_vertices,interior_vertices);


end