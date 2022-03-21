function surface = refineTriangle(surface, triangle_ind)
%REFINETRIANGLE refines a triangle via inserting a central point and edge
%flips
%surface: structure with fields TRIV, X, Y, Z and VERT
%triangle_ind: index of the triangle to be refined

triangle = surface.TRIV(triangle_ind,:);
edge1 = triangle(1:2);
edge2 = triangle(2:3);
edge3 = [triangle(3) triangle(1)];


surface = insertVertexInTriangle(surface, triangle_ind);


%% Flip the edges

surface = edgeFlip(surface, edge1);
surface = edgeFlip(surface, edge2);
surface = edgeFlip(surface, edge3);


end

