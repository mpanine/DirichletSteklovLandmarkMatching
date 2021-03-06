function [W, S, segment_TRIV, boundaries_new, boundary_edges_old] = splitMeshLandmarkSteklov_clean(Shape, landmarks, boundaries)
%Produces the operators involved in the Steklov problem with central landmarks on a segment of a
%shape that has already been refined at the landmarks with refineCircleAroundCenter

%Shape: REFINED structure containing fields .surface.TRIV and .surface.VERT and .W
%landmarks: list of landmarks
%boundary_old: output of refineCircleAroundCenter (indices of the landmark ) vertices of Shape that form the desired segment.

%W: the output stiffness matrix
%S: the output Steklov BC matrix (mass matrix on the boundary)
%segment_TRIV: TRIV of the segment. Uses the initial indexing
%boundary_old: boundary indices in the origial indices 
%boundary_new: boundary indices in the indices relative to the order of "segment"
%boundary_edges_old: defines the boundary connectivity

boundary_old = vertcat(boundaries{:});


%% Computing W

nv = size(Shape.SHAPE.surface.VERT,1);
segment = 1:nv;
segment(landmarks) = []; %List of all non-landmark vertices.

segment_TRIV = findSegmentTriangles(Shape.SHAPE, segment);
preW = cotLaplacianSteklov(Shape.SHAPE.surface.VERT, segment_TRIV);
W = preW(segment, segment);

%% Computing S

[boundary_new, boundary_edges_old] = findCentralLandmarkBoundary(Shape, segment, boundary_old);% Finding the boundary


preS = boundarySteklov(Shape.SHAPE.surface.VERT, boundary_edges_old);
S = preS(segment, segment);

%% Computing boundaries_new

boundaries_new = cell(length(landmarks),1); % Initializes the cell array with the old indexing scheme.

aa = ( 1:length(segment) )';

for i = 1:length(landmarks)
    
%     boundaries{i}
%     min(boundaries{i})
%     max(segment)
    boundaries_new{i} = aa( ismember( segment, boundaries{i}) );
%     boundaries_new{i}



end

end