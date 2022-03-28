function [Src_refined,Tar_refined,fullp2pTarSrc, fullp2pSrcTar,fullp2pTarSrc_ZO, fullp2pSrcTar_ZO] = compute_steklov_map(Src, Src_landmarks, Tar, Tar_landmarks, Steklov_settings)
%COMPUTE_STEKLOV_MAP Summary of this function goes here
%   Detailed explanation goes here

    %% Refine Shapes (landmark circles) and compute all relevant bases
    Src_refined = ComputePrincipledSteklovAll(Src, Src_landmarks, Steklov_settings);
    Tar_refined = ComputePrincipledSteklovAll(Tar, Tar_landmarks, Steklov_settings);

    %% Compute the boundary correspondence -- Mutltiple methods available
    %fprintf('refinement done\n');

    switch Steklov_settings.InitialGuess

        case 'trivial' %Trivial circle matching.
            [fullp2pTarSrc, fullp2pSrcTar] = PrincipledSteklovBoundaryMatch_trivial(Src_refined, Tar_refined, Steklov_settings); 

        case 'landmark_harmonics' % Guess using landmark harmonics.
            [fullp2pTarSrc, fullp2pSrcTar] = PrincipledLandmarkHarmonicInitialGuess(Src_refined,  Tar_refined);

        case 'normal_derivatives' %Normal derivatives of landmark harmonics.
            [fullp2pTarSrc, fullp2pSrcTar] = PrincipledSteklovBoundaryMatch_normal_derivatives(Src_refined,  Tar_refined, Steklov_settings);

        case 'conformal_energy'
            [fullp2pTarSrc, fullp2pSrcTar] = PrincipledSteklovBoundaryMatch_conformal_energy(Src_refined,  Tar_refined, Steklov_settings);

    end




    %% ZoomOut Procedure -- Produces the desired output

    [fullp2pTarSrc_ZO, fullp2pSrcTar_ZO] = Principled_findFullP2P_ZO(Src_refined, Src_landmarks, Tar_refined, Tar_landmarks, fullp2pTarSrc, fullp2pSrcTar, Steklov_settings);
end

