function h_basis = ComputeCentralHarmonicBasis(Shape, landmarks)
%Computes the harmonic basis for a given shape using central FE and landmarks by linear solve.
%Shape: shape structure (MeshInfo output)
%landmarks: indices of the landmarks

h_basis = zeros(size(Shape.W,1),length(landmarks));

[Wii, Wib, ~] = boundary_poisson_split_central(Shape, landmarks);


for i = 1:length(landmarks)
    
    landmark_values = zeros(length(landmarks),1);
    landmark_values(i) = 1;
    
    pre_h_basis = mldivide(Wii, -Wib*landmark_values);
    h_basis(:,i) = BoundaryReinsert(pre_h_basis, landmarks, landmark_values );
    
    
end


end