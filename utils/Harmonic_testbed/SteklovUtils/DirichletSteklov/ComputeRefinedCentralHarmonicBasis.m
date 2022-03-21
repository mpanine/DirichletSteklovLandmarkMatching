function h_basis = ComputeRefinedCentralHarmonicBasis(S_refined, S_refined_boundaries, landmarks)
%Computes the harmonic basis for a given shape and landmarks with small circles of radius radii at the landmarks.
%Shape: shape structure (MeshInfo output)
%landmarks: indices of the landmarks
%radii: radii of the circles at the landmarks
%numsubdivisions: number of subdivisions of each triangle at a landmark

nv = size(S_refined.SHAPE.W,1);
h_basis = zeros(size(S_refined.SHAPE.W,1),length(landmarks));

boundary_long = [];
boundary_section_length = zeros(size(landmarks));

for i = 1:length(S_refined_boundaries)
    boundary_long = [boundary_long; S_refined_boundaries{i}];
    boundary_section_length(i) = length(S_refined_boundaries{i});
end

boundary_section_length = boundary_section_length(:);
boundary_section_ends = [0; cumsum(boundary_section_length)];

[Wii, Wib, ~] = boundary_poisson_split(S_refined.SHAPE, boundary_long);


for i = 1:length(landmarks)
    
    landmark_values = zeros(length(boundary_long),1);    
    landmark_values(boundary_section_ends(i) + 1: boundary_section_ends(i+1) ) = 1;
    
    h_basis(:,i) = BoundaryReinsert(mldivide(Wii, -Wib*landmark_values), boundary_long, landmark_values );
       
    
end


end