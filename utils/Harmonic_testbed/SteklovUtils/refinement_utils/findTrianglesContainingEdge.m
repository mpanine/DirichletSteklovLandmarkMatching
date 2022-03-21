function [tri_indices, edge_variants] = findTrianglesContainingEdge(surface, edge)
%FINDTRIANGLESCONTAININGEDGE finds the indices of the triangles containing a
%certain edge. Outputs an empty array if no such triangle exists.
%surface: structure with fields TRIV, X, Y, Z and VERT
%edge: indices of the summits of the edge


tri_indices = []; %Future output
edge_variants = []; %Number of the 'edge realization' that is found in the triangle

edge_realizations = edgeRealizations(edge);

for ind = 1:size(surface.TRIV, 1)

    check = sum( surface.TRIV(ind,:) == edge_realizations ,2);
    [a,b] = max(check);
        
    if a == 2
        
        tri_indices = [tri_indices ind];
        edge_variants = [edge_variants b];
        
    end


end


end

