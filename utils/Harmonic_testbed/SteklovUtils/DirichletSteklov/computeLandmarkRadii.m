function landmarks_radii = computeLandmarkRadii(Src, Src_landmarks, Tar, Tar_landmarks, Steklov_settings)
%computeLandmarkRadii: computes the radii used to refine the circles around
%the landmarks.
%Src: Source structure
%Src_landmarks: List of used source landmarks
% Tar: Target structure
% Tar_landmarks: List of used target landmarks
% Steklov_settings: Structure containing the general settings for the Steklov approach

Src_shortest_edges = findShortestEdgesAtLandmarks(Src, Src_landmarks);
Tar_shortest_edges = findShortestEdgesAtLandmarks(Tar, Tar_landmarks);
    
landmarks_radii = zeros(size(Src_shortest_edges));
for land_ind = 1:length(landmarks_radii)
    
    landmarks_radii(land_ind) = min(Src_shortest_edges(land_ind), Tar_shortest_edges(land_ind) );
    
end

landmarks_radii = Steklov_settings.radii_factor * landmarks_radii;











end