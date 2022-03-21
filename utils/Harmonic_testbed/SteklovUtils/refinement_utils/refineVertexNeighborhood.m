function surface = refineVertexNeighborhood(surface, vertex, num_circles)
%REFINEVERTEXNEIGHBORHOOD refines a neighborhood of a vertex.
%surface: structure containing fields TRIV, VERT etc
%vertex: index of the vertex of interest
%num_circles: size of the considered neighborhood. Effectively the minimum
%is 1.


tri_list = findTriangleNeigborhood(surface, vertex, num_circles);
edge_list = findEdgesTriangleList(surface, tri_list);


%% Inserting new vertices

for t = 1:length(tri_list)

    surface = insertVertexInTriangle(surface, tri_list(t) );

end


%% Flipping the edges

for e = 1:size(edge_list,1)

    surface = edgeFlip(surface, edge_list(e,:) );

end


end