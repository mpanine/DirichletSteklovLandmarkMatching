function stiff = edgeStiffness(Shape, edge, landmark, radius)
%EDGESTIFFNNESS computes the stiffness associated to an edge
%TRIV: list of triangles

stiff = 0;


triangle_list = findTrianglesWithEdge(Shape.surface.TRIV, edge);

for i =  1:length(triangle_list)

    
    edge_lengths = triangleEdgeLengths(Shape, triangle_list(i));
    stiff =  stiff + triangleStiffness(Shape.surface.TRIV(triangle_list(i), :), edge_lengths, landmark, edge, radius);


end






end