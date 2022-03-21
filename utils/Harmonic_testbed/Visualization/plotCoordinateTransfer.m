function plotCoordinateTransfer(S1, S2, p2p_S2_to_S1, landmarks)
%This should be self-explanatory

% figure
trimesh(S2.surface.TRIV, S1.surface.X(p2p_S2_to_S1), S1.surface.Y(p2p_S2_to_S1), S1.surface.Z(p2p_S2_to_S1), ...
'FaceColor','interp', ...
'FaceAlpha', 1, 'EdgeColor', 'k'); axis equal;

if nargin > 3
    
    hold on
    
    plot3(S1.surface.X(landmarks),S1.surface.Y(landmarks),S1.surface.Z(landmarks),'.','Color','r','MarkerSize',20)
    
    
    
end







end