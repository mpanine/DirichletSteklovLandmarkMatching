function [ pF_TarSrc, pF_SrcTar ] = refine_harmonic( Src, Tar, BASIS, pF_SrcTar_PreRef, pF_TarSrc_PreRef, varargin )
    %REFINE_HARMONIC integration in the testbed of the Harmonic ZoomOut
    %refinement.
    import basis.harmonic.isometricZoomOut.*
    
    basis_settings = BASIS.settings;
    numBasisFun = basis_settings.numBasisFun;
    NumControlLandmarks = basis_settings.parameters.numControlLandmarks;
    
    p = inputParser;
    addOptional(p,'start_size',1);
    addOptional(p,'end_size',numBasisFun);
    addOptional(p,'step_size',1);
    parse(p,varargin{:});
    start_size = p.Results.start_size;
    end_size = p.Results.end_size;
    step_size = p.Results.step_size;
    
    Src_Harmonic_PreRefined = BASIS.SRC.basis;
    Tar_Harmonic_PreRefined = BASIS.TAR.basis;
    Src_LB_PreRefined = BASIS.SRC.basis_LB;
    Tar_LB_PreRefined = BASIS.TAR.basis_LB;
    
    landmarks_Src = Src.landmarkIndices;
    landmarks_Tar = Tar.landmarkIndices;
    
    Src.landmarkIndices(numBasisFun+NumControlLandmarks+1:end) = [];
    Tar.landmarkIndices(numBasisFun+NumControlLandmarks+1:end) = [];
    
    [pF_SrcTar, pF_TarSrc] = isometricZoomOut(Src, Tar, Src_Harmonic_PreRefined, ...
        Tar_Harmonic_PreRefined, Src_LB_PreRefined, Tar_LB_PreRefined,...
        pF_SrcTar_PreRef, pF_TarSrc_PreRef, start_size, end_size, step_size, ...
        Src.landmarkIndices(1:end-NumControlLandmarks), ...
        Tar.landmarkIndices(1:end-NumControlLandmarks));
    
    Src.landmarkIndices = landmarks_Src;
    Tar.landmarkIndices = landmarks_Tar;
    
end

