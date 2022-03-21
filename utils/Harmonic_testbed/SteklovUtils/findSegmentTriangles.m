function TRIV = findSegmentTriangles(Shape, segment)
% Finds the triangles that are entirely contained in "segment"

%Shape: structure containing fields .surface.TRIV and .surface.VERT
%segment: list of vertices of Shape that form the desired segment.

pre_our_triangles = sum( ismember(Shape.surface.TRIV, segment), 2); % Number of vertices in segment

TRIV = Shape.surface.TRIV( pre_our_triangles == 3,:);




end