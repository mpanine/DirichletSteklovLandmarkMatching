function [fullp2pTarSrc, fullp2pSrcTar] = findFullP2P_ZO_clean(Src_refined, Src_landmarks, Tar_refined, Tar_landmarks, FmapSrcTar, FmapTarSrc, Steklov_settings)
%findFullP2P_LB_V4: Converts the functional maps obtained at the boundary into point to point
%maps in the bulk. Uses both Dirichlet-Steklov eigenvunctions and Dirichlet LB eigenfunctions.
% This version works entirely in the bulk of the shape in order to
% guarantee an interaction between the DS and LB parts of the map.



%% Settings: make this into arguments, eventually

original_only_at_end = true; % If true, the final p2p maps are between the original meshes, including the known landmark correspondence.

W_normalize_DS_ef = true; %False DOES SOEMTHIGN ELSE NOW> %Set to true in order to use W-normalized DS-eigenfunctions. Currently "true" is recommended.


%% Extract settings from Steklov_settings

SteklovWeight = sqrt(Steklov_settings.weight_Steklov);
ProperWeight = sqrt(Steklov_settings.weight_Proper);
OrthoWeight = sqrt(Steklov_settings.weight_Orthonormality);
BijectWeight = sqrt(Steklov_settings.weight_Bijectivity);

LBProperWeight = sqrt(Steklov_settings.weight_LB_Proper);
LBOrthoWeight = sqrt(Steklov_settings.weight_LB_Orthonormality);
LBBijectWeight = sqrt(Steklov_settings.weight_LB_Bijectivity);


ZO_start = Steklov_settings.LB_ZO_start; 
ZO_end = Steklov_settings.LB_ZO_end; 
ZO_step = Steklov_settings.LB_ZO_step;

%% Preliminaries

num_bounds = size(Src_refined.STEKLOV.evecs,3);
% num_eigs = size(Src_refined.STEKLOV.evecs,2);


% Src_non_landmarks = 1:Src_refined.SHAPE.nv;
% Src_non_landmarks(Src_landmarks) = [];
% 
% Tar_non_landmarks = 1:Tar_refined.SHAPE.nv;
% Tar_non_landmarks(Tar_landmarks) = [];


%% H1 - normalizing the DS eigenfunctions and correcting the Fmaps accordingly

Src_sqrt_evals = zeros(size(Src_refined.STEKLOV.evals));
Tar_sqrt_evals = zeros(size(Tar_refined.STEKLOV.evals));

Src_inv_sqrt_evals = zeros(size(Src_refined.STEKLOV.evals));
Tar_inv_sqrt_evals = zeros(size(Tar_refined.STEKLOV.evals));

Src_evecs_normalized = Src_refined.STEKLOV.evecs;
Tar_evecs_normalized = Tar_refined.STEKLOV.evecs;

for i = 1:num_bounds %Normalizing the bulk DS eigenvectors and correcting the relevant Fmaps.
    
    Src_sqrt_evals(:,:,i) = diag( diag( Src_refined.STEKLOV.evals(:,:,i) ).^0.5 );
    Src_inv_sqrt_evals(:,:,i) = diag( diag( Src_refined.STEKLOV.evals(:,:,i) ).^-0.5 );
    
    Tar_sqrt_evals(:,:,i) = diag( diag( Tar_refined.STEKLOV.evals(:,:,i) ).^0.5 );
    Tar_inv_sqrt_evals(:,:,i) = diag( diag( Tar_refined.STEKLOV.evals(:,:,i) ).^-0.5 );
    
    
    if W_normalize_DS_ef
    
        Src_evecs_normalized(:,:,i) = Src_refined.STEKLOV.evecs(:,:,i) * Src_inv_sqrt_evals(:,:,i); %DO NOT FORGET TO USE THESE IN THE CODE BELOW!!!
        Tar_evecs_normalized(:,:,i) = Tar_refined.STEKLOV.evecs(:,:,i) * Tar_inv_sqrt_evals(:,:,i);

        FmapSrcTar{i} = Tar_sqrt_evals(:,:,i) * FmapSrcTar{i} * Src_inv_sqrt_evals(:,:,i);
        FmapTarSrc{i} = Src_sqrt_evals(:,:,i) * FmapTarSrc{i} * Tar_inv_sqrt_evals(:,:,i);
    
    else %A do nothing quick-fix. NOT FOR NOW
        
%         Src_evecs_normalized(:,:,i) = Src_refined.STEKLOV.evecs(:,:,i); %DO NOT FORGET TO USE THESE IN THE CODE BELOW!!!
%         Tar_evecs_normalized(:,:,i) = Tar_refined.STEKLOV.evecs(:,:,i);
% 
%         FmapSrcTar{i} = FmapSrcTar{i};
%         FmapTarSrc{i} = FmapTarSrc{i};
        
        Src_evecs_normalized(:,:,i) = Src_refined.STEKLOV.evecs(:,:,i) * Src_inv_sqrt_evals(:,:,i); %DO NOT FORGET TO USE THESE IN THE CODE BELOW!!!
        Tar_evecs_normalized(:,:,i) = Tar_refined.STEKLOV.evecs(:,:,i) * Tar_inv_sqrt_evals(:,:,i);

        FmapSrcTar{i} = Tar_sqrt_evals(:,:,i) * FmapSrcTar{i} * Src_inv_sqrt_evals(:,:,i);
        FmapTarSrc{i} = Src_sqrt_evals(:,:,i) * FmapTarSrc{i} * Tar_inv_sqrt_evals(:,:,i);
        
    end
    
end


    
%% Computing the first guess of the Tar -> Src bulk point-to-point maps map.

[SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Src_evecs_normalized, Src_refined.STEKLOV.evals,...
                                                                   Tar_evecs_normalized, Tar_refined.STEKLOV.evals,...
                                                                   FmapSrcTar, FmapTarSrc, SteklovWeight, ProperWeight, BijectWeight); %Function defined below.
fullp2pTarSrc = annquery(SrcFTarSrc(:,:)', TarFTarSrc(:,:)', 1);

%% Computing the first guess of the Src -> Tar point-to-point maps map.
    
[TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Tar_evecs_normalized, Tar_refined.STEKLOV.evals,...
                                                                   Src_evecs_normalized, Src_refined.STEKLOV.evals,...
                                                                   FmapTarSrc, FmapSrcTar, SteklovWeight, ProperWeight, BijectWeight);
fullp2pSrcTar = annquery(TarFSrcTar(:,:,:)', SrcFSrcTar(:,:,:)', 1);  


%% Main iteration

used_zo = ZO_start:ZO_step:ZO_end; % List of ZoomOut sizes iterated over
if used_zo(end) < ZO_end
    used_zo = [used_zo ZO_end]; % Ensures that the end size is always used.
end

for zo = used_zo
    
    %Computing the LB Fmap parts (bulk)
    
    FmapLBSrcTar = Tar_refined.STEKLOV.LB_evecs(:,1:zo)' * Tar_refined.SHAPE.W * Src_refined.STEKLOV.LB_evecs(fullp2pTarSrc,1:zo); %LB part of the map. H1 INNER PRIODUCT!!!
    FmapLBTarSrc = Src_refined.STEKLOV.LB_evecs(:,1:zo)' * Src_refined.SHAPE.W * Tar_refined.STEKLOV.LB_evecs(fullp2pSrcTar,1:zo);

   
    for i = 1:num_bounds
        
        if W_normalize_DS_ef
        
            FmapSrcTar{i} = Tar_refined.STEKLOV.evals(:,:,i) * Tar_evecs_normalized(Tar_refined.STEKLOV.boundaries{i},:,i)'...
                            * Tar_refined.STEKLOV.MassBoundary{i} *...
                            Src_evecs_normalized(fullp2pTarSrc(Tar_refined.STEKLOV.boundaries{i}),:,i);

            FmapTarSrc{i} = Src_refined.STEKLOV.evals(:,:,i) * Src_evecs_normalized(Src_refined.STEKLOV.boundaries{i},:,i)' ...
                            * Src_refined.STEKLOV.MassBoundary{i} * ...
                            Tar_evecs_normalized(fullp2pSrcTar(Src_refined.STEKLOV.boundaries{i}),:,i);
                    
        else
            
            FmapSrcTar{i} = Tar_evecs_normalized(Tar_refined.STEKLOV.boundaries{i},:,i)'...
                            * Tar_refined.STEKLOV.MassBoundary{i} *...
                            Src_evecs_normalized(fullp2pTarSrc(Tar_refined.STEKLOV.boundaries{i}),:,i);

            FmapTarSrc{i} = Src_evecs_normalized(Src_refined.STEKLOV.boundaries{i},:,i)' ...
                            * Src_refined.STEKLOV.MassBoundary{i} * ...
                            Tar_evecs_normalized(fullp2pSrcTar(Src_refined.STEKLOV.boundaries{i}),:,i);
            
        end
        
    end
    

    
    %% Computing the new p2p map
    
    if (zo ~= ZO_end)||(~original_only_at_end) %Normal behavior for all iterations, except possibly the last one.
        
        [SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Src_evecs_normalized, Src_refined.STEKLOV.evals,...
                                                                           Tar_evecs_normalized, Tar_refined.STEKLOV.evals,...
                                                                           FmapSrcTar, FmapTarSrc, SteklovWeight, ProperWeight, BijectWeight);
        [TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Tar_evecs_normalized, Tar_refined.STEKLOV.evals,...
                                                                           Src_evecs_normalized, Src_refined.STEKLOV.evals,...
                                                                           FmapTarSrc, FmapSrcTar, SteklovWeight, ProperWeight, BijectWeight);


        [Src_LB_TarSrc, Tar_LB_TarSrc] = build_LB_EmbeddingsForS2_to_S1_p2p(Src_refined.STEKLOV.LB_evecs(:,1:zo), Tar_refined.STEKLOV.LB_evecs(:,1:zo),...
                                                                            FmapLBSrcTar, FmapLBTarSrc, LBProperWeight, LBOrthoWeight, LBBijectWeight);
        [Tar_LB_SrcTar, Src_LB_SrcTar] = build_LB_EmbeddingsForS2_to_S1_p2p(Tar_refined.STEKLOV.LB_evecs(:,1:zo), Src_refined.STEKLOV.LB_evecs(:,1:zo),...
                                                                            FmapLBTarSrc, FmapLBSrcTar, LBProperWeight, LBOrthoWeight, LBBijectWeight);

        fullp2pTarSrc = annquery([SrcFTarSrc, Src_LB_TarSrc]', [TarFTarSrc, Tar_LB_TarSrc]', 1);
        fullp2pSrcTar = annquery([TarFSrcTar, Tar_LB_SrcTar]', [SrcFSrcTar, Src_LB_SrcTar]', 1);  
        
    else %Compute the final p2p map on the interior. Activates when zo = ZO_end AND original_only_at_end is true (changes the last iteration)
        
        [SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Src_evecs_normalized(Src_refined.STEKLOV.Interior,:,:), Src_refined.STEKLOV.evals,...
                                                                           Tar_evecs_normalized(Tar_refined.STEKLOV.Interior,:,:), Tar_refined.STEKLOV.evals,...
                                                                           FmapSrcTar, FmapTarSrc, SteklovWeight, ProperWeight, BijectWeight);
        [TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Tar_evecs_normalized(Tar_refined.STEKLOV.Interior,:,:), Tar_refined.STEKLOV.evals,...
                                                                           Src_evecs_normalized(Src_refined.STEKLOV.Interior,:,:), Src_refined.STEKLOV.evals,...
                                                                           FmapTarSrc, FmapSrcTar, SteklovWeight, ProperWeight, BijectWeight);


        [Src_LB_TarSrc, Tar_LB_TarSrc] = build_LB_EmbeddingsForS2_to_S1_p2p(Src_refined.STEKLOV.LB_evecs(Src_refined.STEKLOV.Interior,1:zo),...
                                                                            Tar_refined.STEKLOV.LB_evecs(Tar_refined.STEKLOV.Interior,1:zo),...
                                                                            FmapLBSrcTar, FmapLBTarSrc, LBProperWeight, LBOrthoWeight, LBBijectWeight);
        [Tar_LB_SrcTar, Src_LB_SrcTar] = build_LB_EmbeddingsForS2_to_S1_p2p(Tar_refined.STEKLOV.LB_evecs(Tar_refined.STEKLOV.Interior,1:zo),...
                                                                            Src_refined.STEKLOV.LB_evecs(Src_refined.STEKLOV.Interior,1:zo),...
                                                                            FmapLBTarSrc, FmapLBSrcTar, LBProperWeight, LBOrthoWeight, LBBijectWeight);

        fullp2pTarSrc = annquery([SrcFTarSrc, Src_LB_TarSrc]', [TarFTarSrc, Tar_LB_TarSrc]', 1);
        fullp2pSrcTar = annquery([TarFSrcTar, Tar_LB_SrcTar]', [SrcFSrcTar, Src_LB_SrcTar]', 1); 
        
        
        %Reinsert the known boundary correspondence.
        fullp2pSrcTar = reinsertBoundaryCorrespondence(fullp2pSrcTar, Src_landmarks, Tar_landmarks, Tar_refined.STEKLOV.Interior);                                                         
        fullp2pTarSrc = reinsertBoundaryCorrespondence(fullp2pTarSrc, Tar_landmarks, Src_landmarks, Src_refined.STEKLOV.Interior);  
        
 
    
    end
    

    
end
    




    
    
    
end






function [S1FS2S1, S2FS2S1] = build_DS_EmbeddingsForS2_to_S1_p2p_boundary(S1_evecs, S1_evals, S2_evecs, S2_evals, FmapS1S2, FmapS2S1, SteklovWeight, ProperWeight, OrthoWeight, BijectWeight)
%Builds the DS embeddings for a S2_to_S1 p2p map for the Dirichlet-Steklov
%part only. (Boundary matching)
%THIS WORKS FOR ONE BOUNDARY COMPONENT AT A TIME


        S1FS2S1 = SteklovWeight * S1_evecs*S1_evals +... %Steklov term
               ProperWeight * S1_evecs+... %"Proper" functional map term (pullback of p2p map)
               OrthoWeight * S1_evecs*FmapS1S2'+... %Orthonormality term.
               BijectWeight * S1_evecs*FmapS2S1; %Bijectivity term

        S2FS2S1 = SteklovWeight * S2_evecs*S2_evals*FmapS1S2 +... %Steklov term
               ProperWeight * S2_evecs*FmapS1S2+... %"Proper" functional map term (pullback of p2p map)
               OrthoWeight * S2_evecs+...%Orthonormality term.
               BijectWeight * S2_evecs; %Bijectivity term


end





function [S1FS2S1, S2FS2S1] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(S1_evecs, S1_evals, S2_evecs, S2_evals, FmapS1S2, FmapS2S1, SteklovWeight, ProperWeight, BijectWeight)
%Builds the DS embeddings for a S2_to_S1 p2p map for the Dirichlet-Steklov
%part only. (Bulk matching, implements bulk orthonormal version of the basis)

% TODO: (potentially) implement the corrected version of the orthonormality term

    num_bounds = size(S1_evecs,3);
    num_eigs = size(S1_evecs,2);
    

    S1FS2S1 = zeros(size(S1_evecs,1),num_bounds*num_eigs); %Initialize the arrays for holding the embeddings over which the nearest neighbor search is performed.
    S2FS2S1 = zeros(size(S2_evecs,1),num_bounds*num_eigs);


    for i = 1:num_bounds

        start_ind = (i-1)*num_eigs + 1;
        end_ind = i*num_eigs;

        S1FS2S1(:,start_ind:end_ind) = SteklovWeight * S1_evecs(:,:,i)*S1_evals(:,:,i) +... %Steklov term
               ProperWeight * S1_evecs(:,:,i)+... %"Proper" functional map term (pullback of p2p map)
               BijectWeight * S1_evecs(:,:,i)*FmapS2S1{i}; %Bijectivity term

        S2FS2S1(:,start_ind:end_ind) = SteklovWeight * S2_evecs(:,:,i)*S2_evals(:,:,i)*FmapS1S2{i} +... %Steklov term
               ProperWeight * S2_evecs(:,:,i)*FmapS1S2{i}+... %"Proper" functional map term (pullback of p2p map)
               BijectWeight * S2_evecs(:,:,i); %Bijectivity term


    end




end







function [S1FS2S1, S2FS2S1] = build_LB_EmbeddingsForS2_to_S1_p2p(S1_LB, S2_LB, FmapLBS1S2, FmapLBS2S1, LBProperWeight, LBOrthoWeight, LBBijectWeight)

             
        S1FS2S1 = LBProperWeight * S1_LB+... %"Proper" functional map term (pullback of p2p map)
                  LBOrthoWeight * S1_LB*FmapLBS1S2'+... %Orthonormality term. -- Should promote conformality (CHECK).
                  LBBijectWeight * S1_LB * FmapLBS2S1; %Bijectivity term

        S2FS2S1 = LBProperWeight * S2_LB * FmapLBS1S2+... %"Proper" functional map term (pullback of p2p map)
                  LBOrthoWeight * S2_LB+...%Orthonormality term.
                  LBBijectWeight * S2_LB; %Bijectivity term

end












