function [fullp2pTarSrc, fullp2pSrcTar] = findFullP2P_LB_V4(Src_evecs, Src_evals, Src_landmarks, Src_LB_eF, Src_refined_boundaries, Src, Src_refined, Src_Interior, Src_Mass_Boundary, Tar_evecs, Tar_evals, Tar_landmarks, Tar_LB_eF, Tar_refined_boundaries, Tar, Tar_refined, Tar_Interior, Tar_Mass_Boundary, FmapSrcTar, FmapTarSrc)
%findFullP2P_LB_V4: Converts the functional maps obtained at the boundary into point to point
%maps in the bulk. Uses both Dirichlet-Steklov eigenvunctions and Dirichlet LB eigenfunctions.
% This version works entirely in the bulk of the shape in order to
% guarantee an interaction between the DS and LB parts of the map.

% The output maps work on the REFINED SHAPES


%The LB-eigenfunctions are used in a Zoom-Out type procedure.
%IMPORTANT: THE LB EIGENFUNCTIONS MUST BE H1 (NOT L2) NORMALIZED. (i.e. L2 normalized and divided by the sqrt of the eigenvalue)


%% Settings: make this into arguments, eventually

interior_only_at_end = true; % If true, the final p2p maps are only between original interior vertices.

ZO_start = 5;
ZO_end = 100;%size(Src_LB_eF,2);
ZO_step = 5;

LB_weight = 5; % A weight of 1 is the default version. 5 works better.

%% Preliminaries

num_bounds = size(Src_evecs,3);
num_eigs = size(Src_evecs,2);


% Src_non_landmarks = 1:Src_refined.SHAPE.nv;
% Src_non_landmarks(Src_landmarks) = [];
% 
% Tar_non_landmarks = 1:Tar_refined.SHAPE.nv;
% Tar_non_landmarks(Tar_landmarks) = [];

%% Restricting the LB eigenfunctions to the interior

% Src_LB_eF_original = Src_LB_eF(1:Src.SHAPE.nv, :); % Restriction to original vertices. Used to build the LB Fmap.
% Tar_LB_eF_original = Tar_LB_eF(1:Tar.SHAPE.nv, :);
% 
% Src_LB_eF = Src_LB_eF(Src_Interior,:); % Restriction to the interior. Used for p2p matching.
% Tar_LB_eF = Tar_LB_eF(Tar_Interior,:);



%% H1 - normalizing the DS eigenfunctions and correcting the Fmaps accordingly

Src_sqrt_evals = zeros(size(Src_evals));
Tar_sqrt_evals = zeros(size(Tar_evals));

Src_inv_sqrt_evals = zeros(size(Src_evals));
Tar_inv_sqrt_evals = zeros(size(Tar_evals));

for i = 1:num_bounds %Normalizing the bulk DS eigenvectors and correcting the relevant Fmaps.
    
    Src_sqrt_evals(:,:,i) = diag( diag( Src_evals(:,:,i) ).^0.5 );
    Src_inv_sqrt_evals(:,:,i) = diag( diag( Src_evals(:,:,i) ).^-0.5 );
    
    Tar_sqrt_evals(:,:,i) = diag( diag( Tar_evals(:,:,i) ).^0.5 );
    Tar_inv_sqrt_evals(:,:,i) = diag( diag( Tar_evals(:,:,i) ).^-0.5 );
    
    Src_evecs(:,:,i) = Src_evecs(:,:,i) * Src_inv_sqrt_evals(:,:,i);
    Tar_evecs(:,:,i) = Tar_evecs(:,:,i) * Tar_inv_sqrt_evals(:,:,i);
    
    FmapSrcTar{i} = Tar_sqrt_evals(:,:,i) * FmapSrcTar{i} * Src_inv_sqrt_evals(:,:,i);
    FmapTarSrc{i} = Src_sqrt_evals(:,:,i) * FmapTarSrc{i} * Tar_inv_sqrt_evals(:,:,i);
    
end


    
%% Computing the first guess of the Tar -> Src bulk point-to-point maps map.

[SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Src_evecs, Src_evals, Tar_evecs, Tar_evals, FmapSrcTar, FmapTarSrc); %Function defined below.
fullp2pTarSrc = annquery(SrcFTarSrc(:,:)', TarFTarSrc(:,:)', 1);

%% Computing the first guess of the Src -> Tar point-to-point maps map.
    
[TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Tar_evecs, Tar_evals, Src_evecs, Src_evals, FmapTarSrc, FmapSrcTar);
fullp2pSrcTar = annquery(TarFSrcTar(:,:,:)', SrcFSrcTar(:,:,:)', 1);  

    
%% ZoomOut-like refinement

% figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])

% for i = 1:num_bounds
%    
%    subplot(1,2,1)
%    plot(boundary_p2pSrcTar); title('Src -> Tar boundary p2p')
%    hold on
%     
%    subplot(1,2,2)
%    plot(boundary_p2pTarSrc); title('Tar -> Src boundary p2p')
%    hold on
%    
%    pause(5)
%    
% end

for zo = ZO_start:ZO_step:ZO_end
    
    %Computing the LB Fmap parts (bulk)
    
    FmapLBSrcTar = Tar_LB_eF(:,1:zo)' * Tar_refined.SHAPE.W * Src_LB_eF(fullp2pTarSrc,1:zo); %LB part of the map. H1 INNER PRIODUCT!!!
    FmapLBTarSrc = Src_LB_eF(:,1:zo)' * Src_refined.SHAPE.W * Tar_LB_eF(fullp2pSrcTar,1:zo);

  
    %Computing the Dirichlet-Steklov Fmap part (bulk, rather than boundary)
    
%     for i = 1:num_bounds
%         
%         FmapSrcTar{i} = Tar_evecs(Tar_non_landmarks,:,i)' * Tar_refined.SHAPE.W(Tar_non_landmarks,Tar_non_landmarks) * Src_evecs(fullp2pTarSrc(Tar_non_landmarks),:,i);
%         FmapTarSrc{i} = Src_evecs(Src_non_landmarks,:,i)' * Src_refined.SHAPE.W(Src_non_landmarks,Src_non_landmarks) * Tar_evecs(fullp2pSrcTar(Src_non_landmarks),:,i);
%         
%     end
    
    %Computing the Dirichlet-Steklov Fmap part (Boundary-based bulk H1 product)
    
    for i = 1:num_bounds
        
        FmapSrcTar{i} = Tar_evals(:,:,i) * Tar_evecs(Tar_refined_boundaries{i},:,i)' * Tar_Mass_Boundary{i} * Src_evecs(fullp2pTarSrc(Tar_refined_boundaries{i}),:,i);
        FmapTarSrc{i} = Src_evals(:,:,i) * Src_evecs(Src_refined_boundaries{i},:,i)' * Src_Mass_Boundary{i} * Tar_evecs(fullp2pSrcTar(Src_refined_boundaries{i}),:,i);
        
    end
    
    
%% Useful plot (inv check etc)
    
%     clf
% 
%     subplot(2,3,1)
%     imagesc(combineFmapCellToMatrix(FmapSrcTar) ); title('Steklov Src->Tar'); colorbar
%     subplot(2,3,2)
%     imagesc(combineFmapCellToMatrix(FmapTarSrc) ); title('Steklov Tar->Src'); colorbar
%     subplot(2,3,3)
%     imagesc(combineProductFmapCellToMatrix(FmapTarSrc, FmapSrcTar, false )); title('Steklov Inv check'); colorbar
%     
%     subplot(2,3,4)
%     imagesc(FmapLBSrcTar); title('LB Src->Tar'); colorbar
%     subplot(2,3,5)
%     imagesc(FmapLBTarSrc); title('LB Tar->Src'); colorbar
%     subplot(2,3,6)
%     imagesc(FmapLBSrcTar * FmapLBTarSrc); title('LB Inv check'); colorbar
%     
%     if zo == ZO_start
%         pause(1)
%     else
%         pause(0.1)
%     end
    
    %% Computing the new p2p map
    
    if (zo ~= ZO_end)||(~interior_only_at_end) %Normal behavior for all iterations, except possibly the last one.
        
        [SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Src_evecs, Src_evals, Tar_evecs, Tar_evals, FmapSrcTar, FmapTarSrc);
        [TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Tar_evecs, Tar_evals, Src_evecs, Src_evals, FmapTarSrc, FmapSrcTar);


        [Src_LB_TarSrc, Tar_LB_TarSrc] = build_LB_EmbeddingsForS2_to_S1_p2p(Src_LB_eF(:,1:zo), Tar_LB_eF(:,1:zo), FmapLBSrcTar, FmapLBTarSrc);
        [Tar_LB_SrcTar, Src_LB_SrcTar] = build_LB_EmbeddingsForS2_to_S1_p2p(Tar_LB_eF(:,1:zo), Src_LB_eF(:,1:zo), FmapLBTarSrc, FmapLBSrcTar);

        fullp2pTarSrc = annquery([SrcFTarSrc, LB_weight*Src_LB_TarSrc]', [TarFTarSrc, LB_weight*Tar_LB_TarSrc]', 1);
        fullp2pSrcTar = annquery([TarFSrcTar, LB_weight*Tar_LB_SrcTar]', [SrcFSrcTar, LB_weight*Src_LB_SrcTar]', 1);  
        
    else %Compute the final p2p map on the interior. Activates when zo = ZO_end AND interior_only_at_end is true (changes the last iteration)
        
        [SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Src_evecs(Src_Interior,:,:), Src_evals, Tar_evecs(Tar_Interior,:,:), Tar_evals, FmapSrcTar, FmapTarSrc);
        [TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Tar_evecs(Tar_Interior,:,:), Tar_evals, Src_evecs(Src_Interior,:,:), Src_evals, FmapTarSrc, FmapSrcTar);


        [Src_LB_TarSrc, Tar_LB_TarSrc] = build_LB_EmbeddingsForS2_to_S1_p2p(Src_LB_eF(Src_Interior,1:zo), Tar_LB_eF(Tar_Interior,1:zo), FmapLBSrcTar, FmapLBTarSrc);
        [Tar_LB_SrcTar, Src_LB_SrcTar] = build_LB_EmbeddingsForS2_to_S1_p2p(Tar_LB_eF(Tar_Interior,1:zo), Src_LB_eF(Src_Interior,1:zo), FmapLBTarSrc, FmapLBSrcTar);

        fullp2pTarSrc = annquery([SrcFTarSrc, LB_weight*Src_LB_TarSrc]', [TarFTarSrc, LB_weight*Tar_LB_TarSrc]', 1);
        fullp2pSrcTar = annquery([TarFSrcTar, LB_weight*Tar_LB_SrcTar]', [SrcFSrcTar, LB_weight*Src_LB_SrcTar]', 1); 
        
    
    end
    

    
end
    
    
    
    
    
    
end




     
    
    











function [S1FS2S1, S2FS2S1] = build_DS_EmbeddingsForS2_to_S1_p2p_boundary(S1_evecs, S1_evals, S2_evecs, S2_evals, FmapS1S2, FmapS2S1)
%Builds the DS embeddings for a S2_to_S1 p2p map for the Dirichlet-Steklov
%part only. (Boundary matching)
%THIS WORKS FOR ONE BOUNDARY COMPONENT AT A TIME

    num_bounds = size(S1_evecs,3);
    num_eigs = size(S1_evecs,2);

%     S1FS2S1 = zeros(size(S1_evecs,1),num_eigs); %Initialize the arrays for holding the embeddings over which the nearest neighbor search is performed.
%     S2FS2S1 = zeros(size(S2_evecs,1),num_eigs);


        S1FS2S1 = S1_evecs*S1_evals +... %Steklov term
               S1_evecs+... %"Proper" functional map term (pullback of p2p map)
               S1_evecs*FmapS1S2'+... %Orthonormality term.
               S1_evecs*FmapS2S1; %Bijectivity term

        S2FS2S1 = S2_evecs*S2_evals*FmapS1S2 +... %Steklov term
               S2_evecs*FmapS1S2+... %"Proper" functional map term (pullback of p2p map)
               S2_evecs+...%Orthonormality term.
               S2_evecs; %Bijectivity term


   




end





function [S1FS2S1, S2FS2S1] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(S1_evecs, S1_evals, S2_evecs, S2_evals, FmapS1S2, FmapS2S1)
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

        S1FS2S1(:,start_ind:end_ind) = S1_evecs(:,:,i)*S1_evals(:,:,i) +... %Steklov term
               S1_evecs(:,:,i)+... %"Proper" functional map term (pullback of p2p map)
               S1_evecs(:,:,i)*FmapS2S1{i}; %Bijectivity term

        S2FS2S1(:,start_ind:end_ind) = S2_evecs(:,:,i)*S2_evals(:,:,i)*FmapS1S2{i} +... %Steklov term
               S2_evecs(:,:,i)*FmapS1S2{i}+... %"Proper" functional map term (pullback of p2p map)
               S2_evecs(:,:,i); %Bijectivity term


    end




end











function [S1FS2S1, S2FS2S1] = build_LB_EmbeddingsForS2_to_S1_p2p(S1_LB, S2_LB, FmapLBS1S2, FmapLBS2S1)

             
        S1FS2S1 = S1_LB+... %"Proper" functional map term (pullback of p2p map)
                  S1_LB*FmapLBS1S2'+... %Orthonormality term. -- Should promote conformality (CHECK).
                  S1_LB * FmapLBS2S1; %Bijectivity term

        S2FS2S1 = S2_LB * FmapLBS1S2+... %"Proper" functional map term (pullback of p2p map)
                  S2_LB+...%Orthonormality term.
                  S2_LB; %Bijectivity term

end












