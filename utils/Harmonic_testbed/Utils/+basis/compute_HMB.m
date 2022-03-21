function [BASIS] = compute_HMB(Src, Tar, HMB_settings)
    %Compute the eigen vectors of the cotangent Laplacian matrix

    assert(isa(HMB_settings, 'basis.HMB_Settings'), ...
        sprintf('HMB_settings has to be of type basis.HMB_Settings. Got %s instead.', ...
        class(HMB_settings)));

    numBasisFun = HMB_settings.numBasisFun;
    NumControlLandmarks = HMB_settings.parameters.numControlLandmarks;
    intrinsic = HMB_settings.parameters.intrinsic;
    numLBFun = HMB_settings.parameters.numLBFun;
    
    landmarks_Src = Src.landmarkIndices;
    landmarks_Tar = Tar.landmarkIndices;
    
    Src.landmarkIndices(numBasisFun+NumControlLandmarks+1:end) = [];
    Tar.landmarkIndices(numBasisFun+NumControlLandmarks+1:end) = [];
    
    if intrinsic
        Src_shortest_edges = basis.harmonic.findShortestEdgesAtLandmarks(Src, Src.landmarkIndices(1:end-NumControlLandmarks));
        Tar_shortest_edges = basis.harmonic.findShortestEdgesAtLandmarks(Tar, Tar.landmarkIndices(1:end-NumControlLandmarks));
        landmarks_radii = zeros(size(Src_shortest_edges));
        for land_ind = 1:length(landmarks_radii)

            landmarks_radii(land_ind) = min(Src_shortest_edges(land_ind), Tar_shortest_edges(land_ind) );

        end
        landmarks_radii = 0.1 * landmarks_radii;
        
%         try % for shapes with holes
            [Src_Harmonic_PreRefined, Src_LB_PreRefined] = basis.harmonic.intrinsicTriangulations.computeRefinedBases(Src.objpath, Src.nv, Src.landmarkIndices(1:end-NumControlLandmarks), landmarks_radii, numLBFun);
            [Tar_Harmonic_PreRefined, Tar_LB_PreRefined] = basis.harmonic.intrinsicTriangulations.computeRefinedBases(Tar.objpath, Tar.nv, Tar.landmarkIndices(1:end-NumControlLandmarks), landmarks_radii, numLBFun);
            SRC.basis = Src_Harmonic_PreRefined; 
            TAR.basis = Tar_Harmonic_PreRefined;
            SRC.basis_LB = Src_LB_PreRefined;
            TAR.basis_LB = Tar_LB_PreRefined;
%         catch % fall back to intrinsic=false
%             warning(['Shape pair ' Src.name '-' Tar.name ' probably has holes. The intrinsic has been replaced with the none-intrinsic.'])
%             SRC.basis = basis.harmonic.harmonic_utils.ComputeHarmonicBasis(Src, Src.landmarkIndices(1:end-NumControlLandmarks) );
%             TAR.basis = basis.harmonic.harmonic_utils.ComputeHarmonicBasis(Tar, Tar.landmarkIndices(1:end-NumControlLandmarks) );
%         end
    else
        SRC.basis = basis.harmonic.harmonic_utils.ComputeHarmonicBasis(Src, Src.landmarkIndices(1:end-NumControlLandmarks) );
        TAR.basis = basis.harmonic.harmonic_utils.ComputeHarmonicBasis(Tar, Tar.landmarkIndices(1:end-NumControlLandmarks) );
        [LB_Src, ~] = basis.harmonic.isometricZoomOut.computeDirichletLB(Src, Src.landmarkIndices(1:end-NumControlLandmarks), numLBFun);
        [LB_Tar, ~] = basis.harmonic.isometricZoomOut.computeDirichletLB(Tar, Tar.landmarkIndices(1:end-NumControlLandmarks), numLBFun); 
        SRC.basis_LB = LB_Src;
        TAR.basis_LB = LB_Tar;
    end
    
    Src.landmarkIndices = landmarks_Src;
    Tar.landmarkIndices = landmarks_Tar;
    
    SRC.basis_inverse = SRC.basis';
    TAR.basis_inverse = TAR.basis';
    BASIS.SRC = SRC;
    BASIS.TAR = TAR;
end




% try %Temporary fix for shapes with holes
% 
%     [Src_Harmonic_PreRefined, Src_LB_PreRefined] = computeRefinedBases(Src_objpath, Src.SHAPE.nv, Src.SHAPE.landmarkIndices(1:end-NumControlLandmarks), landmarks_radii, numEigs);
%     [Tar_Harmonic_PreRefined, Tar_LB_PreRefined] = computeRefinedBases(Tar_objpath, Tar.SHAPE.nv, Tar.SHAPE.landmarkIndices(1:end-NumControlLandmarks), landmarks_radii, numEigs);
% 
% 
%     pF_H1_PreRef = annquery(Src_Harmonic_PreRefined', Tar_Harmonic_PreRefined', 1);
%     pF_H2_PreRef = annquery(Tar_Harmonic_PreRefined', Src_Harmonic_PreRefined', 1);
% 
%     [pF_IsoZO2_PreRef, pF_IsoZO1_PreRef] = isometricZoomOut(Src.SHAPE, Tar.SHAPE, Src_Harmonic_PreRefined, Tar_Harmonic_PreRefined, Src_LB_PreRefined, Tar_LB_PreRefined,...
%             pF_H2_PreRef, pF_H1_PreRef, standard_ZO_start, numEigs, 1, Src.SHAPE.landmarkIndices(1:end-NumControlLandmarks), Tar.SHAPE.landmarkIndices(1:end-NumControlLandmarks));
% 
% 
% catch
% 
%     warning(['Shape pair ' num2str(i) ': ' shapePairList{i} ' probably has holes. The PreRef has been replaced with standard.'])
% 
%     pF_H1_PreRef = pF_H1;
%     pF_H2_PreRef = pF_H2;
%     pF_IsoZO1_PreRef = pF_IsoZO1;
%     pF_IsoZO2_PreRef = pF_IsoZO2;
% 
% end