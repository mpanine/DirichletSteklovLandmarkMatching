function h_basis = ComputeCentralHarmonicBasis(S, landmarks, radii, numsubdivisions)
%Computes the harmonic basis for a given shape and landmarks with small circles of radius radii at the landmarks.
%Shape: shape structure (MeshInfo output)
%landmarks: indices of the landmarks
%radii: radii of the circles at the landmarks
%numsubdivisions: number of subdivisions of each triangle at a landmark

nv = size(S.W,1);
h_basis = zeros(size(S.W,1),length(landmarks));

[S, boundary] = refineCircleAroundCenter(S, landmarks, radii, numsubdivisions);

boundary_long = [];
boundary_section_length = zeros(size(landmarks));

for i = 1:length(boundary)
    boundary_long = [boundary_long; boundary{i}];
    boundary_section_length(i) = length(boundary{i});
end

boundary_section_length = boundary_section_length(:);
boundary_section_ends = [0; cumsum(boundary_section_length)];

[Wii, Wib, ~] = boundary_poisson_split(S, boundary_long);


for i = 1:length(landmarks)
    
    landmark_values = zeros(length(boundary_long),1);    
    landmark_values(boundary_section_ends(i) + 1: boundary_section_ends(i+1) ) = 1;
    
    pre_h_basis = BoundaryReinsert(mldivide(Wii, -Wib*landmark_values), boundary_long, landmark_values );
    h_basis(:,i) = pre_h_basis(1:nv); %Restrict to original vertices
    
    
end


end