function edges = findEdgesWithVertex(TRIV, vertex)
%FINDEDGESWITHVERTEX finds the edges containing a certain vertex. 
%TRIV: list of triangles -- probably compatible with intrinsic triangulations
%vertex: index of the vertex (can be a vector)


tri_indices = findTrianglesWithVertex(TRIV, vertex);


edges = trianglesToEdges(TRIV(tri_indices, :));

edge_list = [];

for i = 1:length(vertex)
    
     [u,~]  = find(edges == vertex(i));
     edge_list = [edge_list; u];
    
end

edge_list = unique(edge_list);

edges = edges(edge_list,:);



end