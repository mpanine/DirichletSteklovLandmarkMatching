function [fullp2pTarSrc, fullp2pSrcTar] = PrincipledLandmarkHarmonicInitialGuess(Src_refined,  Tar_refined)
% Computes the initial guess of the p2p and Fmaps using the guess matching
% the landmark harmonic functions (kronecker type stuff).

num_bounds = length(Src_refined.STEKLOV.MassBoundary);

H_Src = zeros(size(Src_refined.STEKLOV.FullBasis_segment, 1), num_bounds);
H_Tar = zeros(size(Tar_refined.STEKLOV.FullBasis_segment, 1), num_bounds);

for i = 1:num_bounds
    
    H_Src(:,i) = sum( Src_refined.STEKLOV.MassBoundary{i} * Src_refined.STEKLOV.FullBasis_segment(Src_refined.STEKLOV.boundaries{i},:), 1) * Src_refined.STEKLOV.FullBasis_segment';
    H_Tar(:,i) = sum( Tar_refined.STEKLOV.MassBoundary{i} * Tar_refined.STEKLOV.FullBasis_segment(Tar_refined.STEKLOV.boundaries{i},:), 1) * Tar_refined.STEKLOV.FullBasis_segment';
    
    
    
    
end

    fullp2pTarSrc = annquery(H_Src', H_Tar', 1);
    fullp2pSrcTar = annquery(H_Tar', H_Src', 1);


    



end