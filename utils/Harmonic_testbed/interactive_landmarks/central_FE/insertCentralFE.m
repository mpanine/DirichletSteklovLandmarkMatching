function Shape = insertCentralFE(Shape, landmark, radius)
%INSERTCENTRALFE modifes the stiffness matrix of Shape with a central finite
%element. The new stiffness matrix is called WC.
%Shape: structure with fields Shape.W (stiffness) Shape.surface.TRIV etc
%landmark: index of the landmark
%radius: radius of the center of the landmark

%IMPORTANT: Iterating this function over landmarks will not produce the
%desired result if the landmarks are two edges or less apart.


if isfield(Shape, 'WC') %Reuses the previously computed WC. Allows to iterate this function over mutliple landmarks (exterior loop).
    
    WC = Shape.WC;
    
else %Creates a new stiffness matrix to be modified below.
    
    WC = Shape.W;
    
end


%% Off-Diagonal Elements

adj_tri_list = findTrianglesWithVertex(Shape.surface.TRIV, landmark); % Triangles adjacent to the landmark
adj_tris = Shape.surface.TRIV(adj_tri_list,:);

edges_off_diag = trianglesToEdges(adj_tris); % Edges on which stiffness is to be computed. (off diagonal edges)

for i = 1:size(edges_off_diag,1)
    
    WC(edges_off_diag(i,1), edges_off_diag(i,2) ) = edgeStiffness(Shape, edges_off_diag(i,:), landmark, radius);
    WC(edges_off_diag(i,2), edges_off_diag(i,1) ) = WC(edges_off_diag(i,1), edges_off_diag(i,2) );
    
end


%% Diagonal Elements

vertices = trianglesToVertices(adj_tris);

for i = 1:length(vertices)

    WC(vertices(i), vertices(i)) = vertexStiffness(Shape, vertices(i), landmark, radius);
    

end


%% Make output

Shape.WC = WC;


end

