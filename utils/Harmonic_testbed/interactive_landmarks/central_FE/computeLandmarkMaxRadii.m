function radii = computeLandmarkMaxRadii(Shape, landmarks)
%COMPUTELANDMARKMAXRADII computes the maximal internal radius (length of
%shortest adjacent edge) of a list of landmarks.

radii = zeros(length(landmarks),1);


for i = 1:length(landmarks)

    radii(i) = min( edgeLengths(Shape, findEdgesWithVertex(Shape.surface.TRIV, landmarks(i) ) ) );


end




end