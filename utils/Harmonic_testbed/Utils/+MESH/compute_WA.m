function [ S ] = compute_WA( S )
    %COMPUTE_WA Summary of this function goes here
    %   Detailed explanation goes here
    
    [S.W, A] = cotLaplacian(get_mesh_vtx_pos(S), S.surface.TRIV);
    S.A = sparse(1:S.nv,1:S.nv,A,S.nv,S.nv);
end

