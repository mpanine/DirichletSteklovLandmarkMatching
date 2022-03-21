function shortest_edges = findShortestEdgesAtLandmarks(Shape, landmarks)
%findShortestEdgeAtLandmarks find the lengths of the shortest edges adjacent to the landmarks


shortest_edges = zeros(size(landmarks));




for i = 1:length(landmarks)
    
    [a,~] = find(Shape.Elist == landmarks(i) );    
    edges = Shape.Elist(a,:);
    
    edge_lengths =  sqrt(  sum( ( Shape.surface.VERT(edges(:,1),:) - Shape.surface.VERT(edges(:,2),:) ).^2 , 2 )    );
    
    shortest_edges(i) = min(edge_lengths);

end









end