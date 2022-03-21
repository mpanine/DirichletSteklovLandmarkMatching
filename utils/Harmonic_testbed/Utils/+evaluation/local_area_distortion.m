function [ area_distortion ] = local_area_distortion( Src, Tar )
    %LOCAL_AREA_DISTORTION Computes the local area distortion of the source
    %shape compared to the target shape at landmarks.
    
    assert(isfield(Src, 'landmarkIndices'),...
            sprintf('landmarkIndices needs to be defined on the Source shape %s to allow the computation of the LDB.',...
                    Src.name));
                
    assert(isfield(Tar, 'landmarkIndices'),...
            sprintf('landmarkIndices needs to be defined on the Target shape %s to allow the computation of the LDB.',...
                    Tar.name));
                
    areas_Src = diag(Src.D);
    areas_Tar = diag(Tar.D);
    
    area_distortion = zeros(1, Src.nv);
    area_distortion(Src.landmarkIndices) = max(areas_Src(Src.landmarkIndices)./areas_Tar(Tar.landmarkIndices),...
                          areas_Tar(Tar.landmarkIndices)./areas_Src(Src.landmarkIndices))';
    
    
end

