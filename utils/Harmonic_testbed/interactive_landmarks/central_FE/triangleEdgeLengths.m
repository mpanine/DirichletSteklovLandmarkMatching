function lengths = triangleEdgeLengths(Shape, tri_index)
%triangleEdgeLengths computes the lengths of the edges of a triangle of
%Shape given by tri_index

T = Shape.surface.TRIV(tri_index,:);

V1 = Shape.surface.VERT(T(1),:);
V2 = Shape.surface.VERT(T(2),:);
V3 = Shape.surface.VERT(T(3),:);

lengths = zeros(1,3);

lengths(1) = norm(V1 - V2);
lengths(2) = norm(V2 - V3);
lengths(3) = norm(V3 - V1);





end