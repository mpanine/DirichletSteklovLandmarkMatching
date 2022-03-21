function [tri_indices, edge_variants] = findTrianglesWithEdge(TRIV, edge)
%FINDTRIANGLESWITHEDGE finds the indices of the triangles containing a
%certain edge. Outputs an empty array if no such triangle exists.
%TRIV: list of triangles -- probably compatible with intrinsic triangulations
%edge: indices of the summits of the edge


tri_indices = []; %Future output
edge_variants = []; %Number of the 'edge realization' that is found in the triangle

edge_realizations = edgeRealizations(edge);

for ind = 1:size(TRIV, 1)

    check = sum( TRIV(ind,:) == edge_realizations ,2);
    [a,b] = max(check);
        
    if a == 2
        
        tri_indices = [tri_indices ind];
        edge_variants = [edge_variants b];
        
    end


end


end
