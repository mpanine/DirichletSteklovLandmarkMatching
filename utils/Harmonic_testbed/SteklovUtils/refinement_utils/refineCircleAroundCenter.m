function [S, boundary, normalTriangles] = refineCircleAroundCenter(S, landmarks, radii, numsubdivisions)
%refineCircleAroundCenter refines the mesh around the landmarks in
%order to build a small circle of radius "radius" around it.

%shape: structure containing the field SHAPE
%landmarks: central vertex (a landmark)
%radii: radii of the circles
%numsubdivisions: number of subdivisions of each neighboring triangle

%boundary: cell array of vertices corresponding to the landmarks
%normalTriangles: triangles that have two sides

% figure
% subplot(1,2,1)
% trimesh(S.surface.TRIV, S.surface.X, S.surface.Y, S.surface.Z, ...
%         'FaceColor','interp', ...
%         'FaceAlpha', 1, 'EdgeColor', 'k'); axis equal;
%     hold on;
%     plot3(S.surface.X(landmarks),S.surface.Y(landmarks),S.surface.Z(landmarks),'.','MarkerSize',20,'Color','r');
%     title('Original')

normalTriangles = cell(length(landmarks),1);
boundary = cell(length(landmarks),1);

for l = 1:length(landmarks)

%     tri_list = findTrianglesContainingVertex(S.surface, landmarks(l) );
    tri_list = OrderedTrianglesNearLandmark(S.surface, landmarks(l)  );


    for t = 1:length(tri_list)

        [S, currentNormalTriangles] = refineSingleCentralTriangle(S, tri_list(t), landmarks(l) , radii(l), numsubdivisions(l) );
        normalTriangles{l} = [normalTriangles{l}; currentNormalTriangles];
    
    end
    
    boundary{l} = findVerticesAdjacentToVertex(S, landmarks(l));
    boundary{l}( boundary{l} == landmarks(l) ) = [];

end

shape_settings = shape.Settings('computeMHB', false,...
                        'computeGeoDist', false,...
                        'MHB_settings', basis.MHB_Settings(0),...
                        'findNeigh',false,...
                        'computeNormals', false);
cache_settings = cache.Settings();


% shape = recomputeShape(shape);
S = shape.compute(S, shape_settings, cache_settings);


% subplot(1,2,2)
% trimesh(S.surface.TRIV, S.surface.X, S.surface.Y, S.surface.Z, ...
%         'FaceColor','interp', ...
%         'FaceAlpha', 1, 'EdgeColor', 'k'); axis equal;
%     hold on;
%     plot3(S.surface.X(landmarks),S.surface.Y(landmarks),S.surface.Z(landmarks),'.','MarkerSize',20,'Color','r');
%     title('Cricle Refined')









end