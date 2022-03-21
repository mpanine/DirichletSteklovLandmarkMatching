function neighbors = findVerticesAdjacentToVertex(S, vertex)
%findVerticesAdjacentToVertex produces a list of vertices adjacent to the
%specified vertex
%S: structure containing field surface (X,Y,Z, VERT, TRIV)
%vertex: index of vertex

tri_list = findTrianglesContainingVertex(S.surface, vertex);

all_vertices = S.surface.TRIV(tri_list,:);

neighbors = unique(all_vertices);


end