function [ W, A, area_matrix, sqrt_area ] = compute_area_weights( S )
%COMPUTE_AREA_WEIGHTS computes the area weights of the shape
% Input: Shape S
% Output: [W, A, area_matrix, sqrt_area]
    % A: mixed voronoi area weights
    [W, A] = cotLaplacian(get_mesh_vtx_pos(S), S.surface.TRIV);
    % A: one-ring neighbor area weights
    % S.A = vertexAreas(get_mesh_vtx_pos(S), S.surface.TRIV);

    A = sparse(1:length(A), 1:length(A), A);
    
     area_matrix = diag(A);
     sqrt_area = sqrt(sum(area_matrix));

end

