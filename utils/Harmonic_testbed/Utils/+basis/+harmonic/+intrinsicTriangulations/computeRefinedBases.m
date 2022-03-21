function [HarmonicBasis, LaplaceBasis] = computeRefinedBases(shape_path, numvertices, landmarks, landmarks_radii, numEigs) 
%computeRefinedBases computes the harmonic and LB bases with boundary
%conditions at the landmarks. 
%shape_path: path of the .obj version of the object
%numvertices: number of vertices in the original shape.
%landmarks: list of landmarks
%landmarks_radii
import basis.harmonic.intrinsicTriangulations.*

[W, A] = getIntrinsicOperators( shape_path, landmarks, landmarks_radii );

HarmonicBasis = ComputeHarmonicBasisIntrinsic(W, landmarks);

[LaplaceBasis, ~] = computeIntrinsicDirichletLB(W, A, landmarks, numEigs);


HarmonicBasis = HarmonicBasis(1:numvertices,:); % Remove the additional intrinsic vertices.
LaplaceBasis = LaplaceBasis(1:numvertices,:);


end