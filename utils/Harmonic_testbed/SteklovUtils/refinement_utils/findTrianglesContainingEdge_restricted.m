function [tri_indices, edge_variants] = findTrianglesContainingEdge_restricted(TRIV, edge, start_tri, end_tri)
%FINDTRIANGLESCONTAININGEDGE finds the indices of the triangles containing a
%certain edge. Outputs an empty array if no such triangle exists.
%Only considers a restricted range of triangles from start_tri to end_tri
%TRIV: triangles
%edge: indices of the summits of the edge


tri_indices = []; %Future output
edge_variants = []; %Number of the 'edge realization' that is found in the triangle

edge_realizations = edgeRealizations(edge);

for ind = start_tri:end_tri

    check = sum( TRIV(ind,:) == edge_realizations ,2);
    [a,b] = max(check);
        
    if a == 2
        
        tri_indices = [tri_indices ind];
        edge_variants = [edge_variants b];
        
    end


end


end

