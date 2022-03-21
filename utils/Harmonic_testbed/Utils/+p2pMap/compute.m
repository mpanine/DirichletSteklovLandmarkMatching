function [ pF ] = compute( BASIS_Src, BASIS_Tar, F, p2pMap_settings )
%COMPUTE Computes the p2p map given a shape and a p2pMap.Settings object
    fprintf('Computing p2p map with %s...',p2pMap_settings.type); tic;
    
    switch p2pMap_settings.type
        case 'KNN'
            pF = knnsearch((F*BASIS_Src.basis')', BASIS_Tar.basis);
        case 'ANN'
            if isa(p2pMap_settings.basis_settings_Src,'basis.HMB_Settings')
                pF = annquery(BASIS_Src.basis', BASIS_Tar.basis', 1);
            else
                pF = annquery(F*BASIS_Src.basis', BASIS_Tar.basis', 1);
            end
        otherwise
            pF = knnsearch((F*BASIS_Src.basis')', BASIS_Tar.basis);
    end
    
    t =toc; fprintf('done:%.4fs\n',t);
end

