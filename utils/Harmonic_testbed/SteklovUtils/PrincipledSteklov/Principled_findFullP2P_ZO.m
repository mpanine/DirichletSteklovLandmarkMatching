function [fullp2pTarSrc, fullp2pSrcTar] = Principled_findFullP2P_ZO(Src_refined, Src_landmarks, Tar_refined, Tar_landmarks, fullp2pTarSrc, fullp2pSrcTar, Steklov_settings)
%Principled_findFullP2P_ZO_clean_v2: Improves the previous guess of p2p map
%by a ZoomOut prodecure on the Laplace-Betrami portion of the basis.


%% Settings: make this into arguments, eventually

original_only_at_end = true; % If true, the final p2p maps are between the original meshes, including the known landmark correspondence.


%% Extract settings from Steklov_settings

num_bounds = length(Src_landmarks);

ProperWeight = sqrt(Steklov_settings.weight_Proper);
OrthoWeight = sqrt(Steklov_settings.weight_Orthonormality);
BijectWeight = sqrt(Steklov_settings.weight_Bijectivity);

ZO_start = num_bounds * Steklov_settings.DS_num_eigs + Steklov_settings.LB_ZO_start; 
ZO_end = num_bounds * Steklov_settings.DS_num_eigs + Steklov_settings.LB_ZO_end; 
ZO_step = Steklov_settings.LB_ZO_step;

%% Preliminaries

num_bounds = size(Src_refined.STEKLOV.DS_evecs,3);


%% Main iteration

used_zo = ZO_start:ZO_step:ZO_end; % List of ZoomOut sizes iterated over
if used_zo(end) < ZO_end
    used_zo = [used_zo ZO_end]; % Ensures that the end size is always used.
end

% figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])

for zo = used_zo
    
    %Computing the Fmap restricted to map corresponding subspaces to each
    %other
    
    FmapSrcTar = buildRestrictedFmapFromP2P_segment(Src_refined, Tar_refined, fullp2pTarSrc, Steklov_settings, zo);
    FmapTarSrc = buildRestrictedFmapFromP2P_segment(Tar_refined, Src_refined, fullp2pSrcTar, Steklov_settings, zo);
    
    
  
    

    %% Computing the new p2p map
    
    
                                                           

         [Src_FTarSrc, Tar_FTarSrc] = build_Principled_EmbeddingsForS2_to_S1_p2p(Src_refined.STEKLOV.FullBasis_segment(:,1:zo) , Src_refined.STEKLOV.FullBasisProduct(1:zo,1:zo), Src_refined.STEKLOV.W_segment,...
                                                                                 Tar_refined.STEKLOV.FullBasis_segment(:,1:zo) , Tar_refined.STEKLOV.FullBasisProduct(1:zo,1:zo), Tar_refined.STEKLOV.W_segment,...
                                                                                 FmapSrcTar, FmapTarSrc, OrthoWeight, ProperWeight, BijectWeight, Steklov_settings);       
                                                                              
         [Tar_FSrcTar, Src_FSrcTar] = build_Principled_EmbeddingsForS2_to_S1_p2p(Tar_refined.STEKLOV.FullBasis_segment(:,1:zo) , Tar_refined.STEKLOV.FullBasisProduct(1:zo,1:zo), Tar_refined.STEKLOV.W_segment,...
                                                                                 Src_refined.STEKLOV.FullBasis_segment(:,1:zo) , Src_refined.STEKLOV.FullBasisProduct(1:zo,1:zo), Src_refined.STEKLOV.W_segment,...
                                                                                 FmapTarSrc, FmapSrcTar, OrthoWeight, ProperWeight, BijectWeight, Steklov_settings); 
         
        if (zo ~= ZO_end)||(~original_only_at_end) %Normal behavior for all iterations, except possibly the last one.                                                                     
                                                                             
             fullp2pTarSrc = annquery(Src_FTarSrc' , Tar_FTarSrc', 1);
             fullp2pSrcTar = annquery(Tar_FSrcTar' , Src_FSrcTar', 1);
        
        elseif original_only_at_end % Restrict the final map to the original shape.
            
             Src_relevant = Src_refined.SHAPE.nv_unrefined - num_bounds; % Index of the maximal relevant inde
             Tar_relevant = Tar_refined.SHAPE.nv_unrefined - num_bounds;
            
             fullp2pTarSrc = annquery(Src_FTarSrc(1:Src_relevant, :)' , Tar_FTarSrc(1:Tar_relevant, :)', 1);
             fullp2pSrcTar = annquery(Tar_FSrcTar(1:Tar_relevant, :)' , Src_FSrcTar(1:Src_relevant ,:)', 1);
             
             fullp2pTarSrc = reinsertBoundaryCorrespondence_Principled(fullp2pTarSrc, Tar_landmarks, Src_landmarks, Src_refined.STEKLOV.Interior);
             fullp2pSrcTar = reinsertBoundaryCorrespondence_Principled(fullp2pSrcTar, Src_landmarks, Tar_landmarks, Tar_refined.STEKLOV.Interior);
            
            
        else % Replace the landmarks into the map at the last iteration (into the refined shape)
             
             fullp2pTarSrc = reinsertBoundaryCorrespondence_Principled(fullp2pTarSrc, Tar_landmarks, Src_landmarks, Src_refined.STEKLOV.segment);
             fullp2pSrcTar = reinsertBoundaryCorrespondence_Principled(fullp2pSrcTar, Src_landmarks, Tar_landmarks, Tar_refined.STEKLOV.segment);
             
        end
        
       

         
    
end
    

    
 




    
    
    
end


















