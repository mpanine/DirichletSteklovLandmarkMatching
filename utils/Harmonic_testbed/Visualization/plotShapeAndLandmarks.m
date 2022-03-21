function plotShapeAndLandmarks(SHAPE, landmarks_a, landmarks_b)
%This should be self-explanatory

figure
trimesh(SHAPE.surface.TRIV, SHAPE.surface.X, SHAPE.surface.Y, SHAPE.surface.Z, ...
'FaceColor','interp', ...
'FaceAlpha', 1, 'EdgeColor', [0.9 0.9 0.9]); axis equal;

hold on

plot3(SHAPE.surface.X(landmarks_a),SHAPE.surface.Y(landmarks_a),SHAPE.surface.Z(landmarks_a),'.','Color','r','MarkerSize',10)

plot3(SHAPE.surface.X(landmarks_b),SHAPE.surface.Y(landmarks_b),SHAPE.surface.Z(landmarks_b),'.','Color','b','MarkerSize',20)










end