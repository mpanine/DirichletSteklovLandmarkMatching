function lengths = edgeLengths(Shape, edges)
%edgeLengths computes the lengths of the specified length of the edges

lengths = sqrt( sum(( Shape.surface.VERT(edges(:,1),:) - Shape.surface.VERT(edges(:,2),:)  ).^2, 2) );



end