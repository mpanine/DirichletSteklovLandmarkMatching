function [W, S, segment_TRIV, boundary_old, boundary_new, boundary_edges_old] = splitSegmentSteklov(Shape, segment)
%Produces the operators involved in the Steklov problem on a segment of a
%shape
%Shape: structure containing fields .surface.TRIV and .surface.VERT and .W
%segment: list of vertices of Shape that form the desired segment.

%W: the output stiffness matrix
%S: the output Steklov BC matrix (mass matrix on the boundary)
%segment_TRIV: TRIV of the segment. Uses the initial indexing
%boundary_old: boundary indices in the origial indices 
%boundary_new: boundary indices in the indices relative to the order of "segment"
%boundary_edges_old: defines the boundary connectivity

%% Computing W
segment_TRIV = findSegmentTriangles(Shape, segment);
preW = cotLaplacianSteklov(Shape.surface.VERT, segment_TRIV);
W = preW(segment, segment);

%% Computing S

[boundary_old, boundary_new, boundary_edges_old] = findSegmentBoundary(Shape, segment); % Finding the boundary


preS = boundarySteklov(Shape.surface.VERT, boundary_edges_old);
S = preS(segment, segment);



end