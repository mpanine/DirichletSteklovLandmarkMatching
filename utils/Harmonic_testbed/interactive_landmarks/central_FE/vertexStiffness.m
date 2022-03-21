function stiff = vertexStiffness(Shape, vertex, landmark, radius)
%VERTEXSTIFFNESS computes the stiffness associated to a vertex (diagonal element)
%TRIV: list of triangles

stiff = 0;

triangle_list = findTrianglesWithVertex(Shape.surface.TRIV, vertex);

for i =  1:length(triangle_list)
    
    edge_lengths = triangleEdgeLengths(Shape, triangle_list(i));
    stiff =  stiff + triangleStiffness(Shape.surface.TRIV(triangle_list(i), :), edge_lengths, landmark, [vertex vertex], radius);


end






end