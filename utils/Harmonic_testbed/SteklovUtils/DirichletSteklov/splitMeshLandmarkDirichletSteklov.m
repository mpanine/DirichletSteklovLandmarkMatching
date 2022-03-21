function [W, S, segment_TRIV, boundaries_new, boundary_edges_old] = splitMeshLandmarkDirichletSteklov(Shape, landmarks, boundaries)
%Produces the operators involved in the Dirichlet-Steklov problem with central landmarks on a segment of a
%shape that has already been refined at the landmarks with refineCircleAroundCenter

%The difference between this and the Steklov version is the way it treats
%the boundary components.

%Shape: REFINED structure containing fields .surface.TRIV and .surface.VERT and .W
%landmarks: list of landmarks
%boundaries: cell array in which each cell is the list of verices of esch boundary component.

%W: the output stiffness matrix
%S: the output Steklov BC matrix (mass matrix on the boundary)
%segment_TRIV: TRIV of the segment. Uses the initial indexing
%boundary_old: boundary indices in the origial indices 
%boundaries_new: cell array of boundary indices in the indices relative to the order of "segment"
%boundary_edges_old: defines the boundary connectivity


boundary_old = vertcat(boundaries{:});

%% Computing W

nv = size(Shape.SHAPE.surface.X,1);
segment = 1:nv;
segment(landmarks) = []; %List of all non-landmark vertices.

segment_TRIV = findSegmentTriangles(Shape.SHAPE, segment);
preW = cotLaplacianSteklov(Shape.SHAPE.surface.VERT, segment_TRIV);
W = preW(segment, segment);

%% Computing S

[~, boundary_edges_old] = findCentralLandmarkBoundary(Shape, segment, boundary_old);% Finding the boundary


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