function tri_list = findTrianglesContainingVertex(surface, vertex)
%FINDTRIANGLESCONTAININGVERTEX builds a list of triangles that contain a
%given vertex
%surface: structure with field TRIV
%vertex: index of the center vertex


[tri_list,~] = find(surface.TRIV == vertex);


end

