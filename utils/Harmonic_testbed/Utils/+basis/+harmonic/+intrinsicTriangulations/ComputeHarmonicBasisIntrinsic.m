function h_basis = ComputeHarmonicBasisIntrinsic(W, landmarks)
%Computes the harmonic basis for a given shape and landmarks by linear solve.
%Shape: shape structure (MeshInfo output)
%landmarks: indices of the landmarks
import basis.harmonic.intrinsicTriangulations.*
import basis.harmonic.harmonic_utils.*

h_basis = zeros(size(W,1),length(landmarks));

[Wii, Wib, ~] = boundary_poisson_split_intrinsic(W, W, landmarks);


for i = 1:length(landmarks)
    
    landmark_values = zeros(length(landmarks),1);
    landmark_values(i) = 1;
    
    pre_h_basis = mldivide(Wii, -Wib*landmark_values);
    h_basis(:,i) = BoundaryReinsert(pre_h_basis, landmarks, landmark_values );
    
    
end


end