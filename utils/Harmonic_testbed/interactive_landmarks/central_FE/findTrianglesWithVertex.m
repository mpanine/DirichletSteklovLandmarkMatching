function tri_indices = findTrianglesWithVertex(TRIV, vertex)
%FINDTRIANGLESWITHEDGE finds the indices of the triangles containing a
%certain vertex. Outputs an empty array if no such triangle exists.
%TRIV: list of triangles -- probably compatible with intrinsic triangulations
%vertex: index of the vertex (can be a vector)

tri_indices = [];

for i = 1:length(vertex)

    [indices,~] = find(TRIV == vertex(i)); % This is not the most useful function ever devised.
    tri_indices = [tri_indices; indices];

end

end