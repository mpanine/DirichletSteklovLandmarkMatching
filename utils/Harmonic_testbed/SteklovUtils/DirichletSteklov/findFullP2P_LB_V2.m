function [fullp2pTarSrc, fullp2pSrcTar] = findFullP2P_LB_V2(Src_evecs, Src_evals, Src_landmarks, Src_LB_eF, Src_refined_boundaries, Src, Src_refined, Src_Interior, Src_Mass_Boundary, Tar_evecs, Tar_evals, Tar_landmarks, Tar_LB_eF, Tar_refined_boundaries, Tar, Tar_refined, Tar_Interior, Tar_Mass_Boundary, FmapSrcTar, FmapTarSrc)
%findFullP2P_LB_V2: Converts the functional maps obtained at the boundary into point to point
%maps in the bulk. Uses both Dirichlet-Steklov eigenvunctions and Dirichlet LB eigenfunctions.

% The output maps work on the INTERIOR of the shapes (all original non-landmark vertices)


%The LB-eigenfunctions are used in a Zoom-Out type procedure.
%IMPORTANT: THE LB EIGENFUNCTIONS MUST BE H1 (NOT L2) NORMALIZED. (i.e. L2 normalized and divided by the sqrt of the eigenvalue)

% THIS VERSION IMPLEMENTS SEPARATE BOUNDARY VS INTERIOR MATCHING

%% Settings: make this into arguments, eventually

ZO_start = 20;
ZO_end = size(Src_LB_eF,2);
ZO_step = 20;


%% Preliminaries

num_bounds = size(Src_evecs,3);
num_eigs = size(Src_evecs,2);



%% Restricting the LB eigenfunctions to the interior

Src_LB_eF_original = Src_LB_eF(1:Src.SHAPE.nv, :); % Restriction to original vertices. Used to build the LB Fmap.
Tar_LB_eF_original = Tar_LB_eF(1:Tar.SHAPE.nv, :);

Src_LB_eF = Src_LB_eF(Src_Interior,:); % Restriction to the interior. Used for p2p matching.
Tar_LB_eF = Tar_LB_eF(Tar_Interior,:);



% figure
% subplot(1,2,1)
% imagesc(Src_LB_eF_original' * Src.SHAPE.W * Src_LB_eF_original);
% subplot(1,2,2)
% imagesc(Tar_LB_eF_original' * Tar.SHAPE.W * Tar_LB_eF_original);




%% Building the boundary restriction of the DS eigenfunctions

Src_evecs_boundary = cell(num_bounds,1);
Tar_evecs_boundary = cell(num_bounds,1);

for i = 1:num_bounds
    
    Src_evecs_boundary{i} = Src_evecs(Src_refined_boundaries{i},:,:);
    Tar_evecs_boundary{i} = Tar_evecs(Tar_refined_boundaries{i},:,:);
    
end





    
%% Computing the first guess of the point-to-point maps map.

[SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p(Src_evecs, Src_evals, Tar_evecs, Tar_evals, FmapSrcTar, FmapTarSrc); %Function defined below.
fullp2pTarSrc = annquery(SrcFTarSrc(Src_Interior,:)', TarFTarSrc(Tar_Interior,:)', 1);

boundary_p2pTarSrc = cell(num_bounds,1);

for i = 1:num_bounds
    boundary_p2pTarSrc{i} = annquery(SrcFTarSrc(Src_refined_boundaries{i},:)', TarFTarSrc(Tar_refined_boundaries{i},:)', 1);
end


%% Computing the first guess of the point-to-point maps map.
    
[TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p(Tar_evecs, Tar_evals, Src_evecs, Src_evals, FmapTarSrc, FmapSrcTar);
fullp2pSrcTar = annquery(TarFSrcTar(Tar_Interior,:,:)', SrcFSrcTar(Src_Interior,:,:)', 1);  
    
boundary_p2pSrcTar = cell(num_bounds,1);

for i = 1:num_bounds
    boundary_p2pSrcTar{i} = annquery(TarFSrcTar(Tar_refined_boundaries{i},:)', SrcFSrcTar(Src_refined_boundaries{i},:)', 1);
end

    
%% ZoomOut-like refinement

% figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])

for zo = ZO_start:ZO_step:ZO_end
    
    %Computing the LB Fmap parts (bulk)
    
   
%     FmapLBSrcTar = Tar_LB_eF(:,1:zo)' * Tar.SHAPE.W(Tar_Interior,Tar_Interior) * Src_LB_eF(fullp2pTarSrc,1:zo); %LB part of the map. H1 INNER PRIODUCT!!!
%     FmapLBTarSrc = Src_LB_eF(:,1:zo)' * Src.SHAPE.W(Src_Interior,Src_Interior) * Tar_LB_eF(fullp2pSrcTar,1:zo);
%     
    
    FmapLBSrcTar = Tar_LB_eF_original(:,1:zo)' * Tar.SHAPE.W * BoundaryReinsertBatch(Src_LB_eF(fullp2pTarSrc,1:zo), Tar_landmarks, zeros(size(Tar_landmarks)) ) ; %LB part of the map. H1 INNER PRIODUCT!!!
    FmapLBTarSrc = Src_LB_eF_original(:,1:zo)' * Src.SHAPE.W * BoundaryReinsertBatch(Tar_LB_eF(fullp2pSrcTar,1:zo), Src_landmarks, zeros(size(Src_landmarks)) );


     
    %Computing the Dirichlet-Steklov Fmap part (boundary)
    
    FmapSrcTar = cell(num_bounds,1);
    FmapTarSrc = cell(num_bounds,1);
    
    
    
    
    
    
    for i = 1:num_bounds
        
        TarSS = diag(sum( Tar_refined.SHAPE.S(Tar_refined_boundaries{i}, Tar_refined_boundaries{i}) ));
        SrcSS = diag(sum( Src_refined.SHAPE.S(Src_refined_boundaries{i}, Src_refined_boundaries{i}) ));
        
        TarSS = TarSS/trace(TarSS);
        SrcSS = SrcSS/trace(SrcSS);
        
        FmapSrcTar{i} = Tar_evecs_boundary{i}(:,:,i)' * TarSS * Src_evecs_boundary{i}(boundary_p2pTarSrc{i},:,i);
        FmapTarSrc{i} = Src_evecs_boundary{i}(:,:,i)' * SrcSS * Tar_evecs_boundary{i}(boundary_p2pSrcTar{i},:,i);
        
    end
    
%     for i = 1:num_bounds
%         
%         FmapSrcTar{i} = Tar_evecs(Tar_refined_boundaries{i},:,i)' * Tar_Mass_Boundary{i} * Src_evecs(boundary_p2pTarSrc{i},:,i);
%         FmapTarSrc{i} = Src_evecs(Src_refined_boundaries{i},:,i)' * Src_Mass_Boundary{i} * Tar_evecs(boundary_p2pSrcTar{i},:,i);
%         
%     end
        
    
    
%     clf
    
%     subplot(1,2,1)
%     plot(Src_evecs(fullp2pTarSrc(Tar_refined_boundaries{1}),:,1)); title('Src evecs boundary on Tar'); colorbar
%     subplot(1,2,2)
%     plot(Tar_evecs(fullp2pSrcTar(Src_refined_boundaries{1}),:,1)); title('Tar evecs boundary on Src'); colorbar
%     
%     pause(1)
    
%     subplot(1,2,1)
%     imagesc(Src_refined.SHAPE.S(Src_refined_boundaries{1}, Src_refined_boundaries{1})); title('Src boundary mass'); colorbar
%     subplot(1,2,2)
%     imagesc(Tar_refined.SHAPE.S(Tar_refined_boundaries{1}, Tar_refined_boundaries{1})); title('Tar boundary mass'); colorbar
    
%     subplot(2,3,1)
%     imagesc(FmapSrcTar{1}); title('Steklov Src->Tar'); colorbar
%     subplot(2,3,2)
%     imagesc(FmapTarSrc{1}); title('Steklov Tar->Src'); colorbar
%     subplot(2,3,3)
%     imagesc(FmapTarSrc{1}*FmapSrcTar{1}); title('Steklov Inv check'); colorbar
%     
%     subplot(2,3,4)
%     imagesc(FmapLBSrcTar); title('LB Src->Tar'); colorbar
%     subplot(2,3,5)
%     imagesc(FmapLBTarSrc); title('LB Tar->Src'); colorbar
%     subplot(2,3,6)
%     imagesc(FmapLBSrcTar * FmapLBTarSrc); title('LB Inv check'); colorbar
    
    if i == 1
        pause(5)
    else
        pause(0.1)
    end
    
    % Computing the new p2p map
    
    [SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p(Src_evecs, Src_evals, Tar_evecs, Tar_evals, FmapSrcTar, FmapTarSrc);
    [TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p(Tar_evecs, Tar_evals, Src_evecs, Src_evals, FmapTarSrc, FmapSrcTar);
    
        
    [Src_LB_TarSrc, Tar_LB_TarSrc] = build_LB_EmbeddingsForS2_to_S1_p2p(Src_LB_eF(:,1:zo), Tar_LB_eF(:,1:zo), FmapLBSrcTar, FmapLBTarSrc);
    [Tar_LB_SrcTar, Src_LB_SrcTar] = build_LB_EmbeddingsForS2_to_S1_p2p(Tar_LB_eF(:,1:zo), Src_LB_eF(:,1:zo), FmapLBTarSrc, FmapLBSrcTar);
    
    
    fullp2pTarSrc = annquery([SrcFTarSrc(Src_Interior,:), Src_LB_TarSrc]', [TarFTarSrc(Tar_Interior,:), Tar_LB_TarSrc]', 1);
    fullp2pSrcTar = annquery([TarFSrcTar(Tar_Interior,:), Tar_LB_SrcTar]', [SrcFSrcTar(Src_Interior,:), Src_LB_SrcTar]', 1);  
    
%     fullp2pTarSrc = annquery([SrcFTarSrc(Src_Interior,:)]', [TarFTarSrc(Tar_Interior,:)]', 1);
%     fullp2pSrcTar = annquery([TarFSrcTar(Tar_Interior,:)]', [SrcFSrcTar(Src_Interior,:)]', 1);  

    zo
    
end
    
    
    
    
    
    
end




     
    
    











function [S1FS2S1, S2FS2S1] = build_DS_EmbeddingsForS2_to_S1_p2p(S1_evecs, S1_evals, S2_evecs, S2_evals, FmapS1S2, FmapS2S1)
%Builds the DS embeddings for a S2_to_S1 p2p map for the Dirichlet-Steklov
%part only
%(Confirmed correct)

    num_bounds = size(S1_evecs,3);
    num_eigs = size(S1_evecs,2);

    S1FS2S1 = zeros(size(S1_evecs,1),num_bounds*num_eigs); %Initialize the arrays for holding the embeddings over which the nearest neighbor search is performed.
    S2FS2S1 = zeros(size(S2_evecs,1),num_bounds*num_eigs);


    for i = 1:num_bounds

        start_ind = (i-1)*num_eigs + 1;
        end_ind = i*num_eigs;

        S1FS2S1(:,start_ind:end_ind) = S1_evecs(:,:,i)*S1_evals(:,:,i) +... %Steklov term
               S1_evecs(:,:,i)+... %"Proper" functional map term (pullback of p2p map)
               S1_evecs(:,:,i)*FmapS1S2{i}'+... %Orthonormality term.
               S1_evecs(:,:,i)*FmapS2S1{i}; %Bijectivity term

        S2FS2S1(:,start_ind:end_ind) = S2_evecs(:,:,i)*S2_evals(:,:,i)*FmapS1S2{i} +... %Steklov term
               S2_evecs(:,:,i)*FmapS1S2{i}+... %"Proper" functional map term (pullback of p2p map)
               S2_evecs(:,:,i)+...%Orthonormality term.
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












