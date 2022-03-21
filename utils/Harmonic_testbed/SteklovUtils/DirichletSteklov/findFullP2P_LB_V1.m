function [fullp2pTarSrc, fullp2pSrcTar] = findFullP2P_LB_V1(Src_evecs, Src_evals, Src_LB_eF, Src_refined_boundaries, Src_refined, Tar_evecs, Tar_evals, Tar_LB_eF, Tar_refined_boundaries, Tar_refined, FmapSrcTar, FmapTarSrc)
%findFullP2P_LB_V1: Converts the functional maps obtained at the boundary into point to point
%maps in the bulk. Uses both Dirichlet-Steklov eigenvunctions and Dirichlet LB eigenfunctions.
%The LB-eigenfunctions are used in a Zoom-Out type procedure.
%IMPORTANT: THE LB EIGENFUNCTIONS MUST BE H1 (NOT L2) NORMALIZED. (i.e. L2 normalized and divided by the sqrt of the eigenvalue)

% THIS VERSION DOES NOT IMPLEMENT SEPARATE BOUNDARY VS INTERIOR MATCHING

%% Settings: make this into arguments, eventually

ZO_start = 10;
ZO_end = size(Src_LB_eF,2);
ZO_step = 10;


%% Combining the boundary Dirichlet-Steklov data into one matrix.


num_bounds = size(Src_evecs,3);
num_eigs = size(Src_evecs,2);



    
%% Computing the first guess of the fullp2pTarSrc map.

[SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p(Src_evecs, Src_evals, Tar_evecs, Tar_evals, FmapSrcTar, FmapTarSrc); %Function defined below.
fullp2pTarSrc = annquery(SrcFTarSrc', TarFTarSrc', 1);

% boundary_p2pTarSrc = annquery(SrcFTarSrc', TarFTarSrc', 1);

%% Computing the first guess of the fullp2pSrcTar map.
    
[TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p(Tar_evecs, Tar_evals, Src_evecs, Src_evals, FmapTarSrc, FmapSrcTar);
fullp2pSrcTar = annquery(TarFSrcTar', SrcFSrcTar', 1);  
    
% boundary_fullp2pSrcTar = annquery(TarFSrcTar', SrcFSrcTar', 1); 
    
%% ZoomOut-like refinement

figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])

for zo = ZO_start:ZO_step:ZO_end
    
    %Computing the LB Fmap parts (bulk)
    
   
    FmapLBSrcTar = Tar_LB_eF(:,1:zo)' * Tar_refined.SHAPE.A * Src_LB_eF(fullp2pTarSrc,1:zo); %LB part of the map
    FmapLBTarSrc = Src_LB_eF(:,1:zo)' * Src_refined.SHAPE.A * Tar_LB_eF(fullp2pSrcTar,1:zo);
    
    
    %Computing the Dirichlet-Steklov Fmap part (boundary)
    
    FmapSrcTar = cell(num_bounds,1);
    FmapTarSrc = cell(num_bounds,1);
    
    for i = 1:num_bounds
        
        FmapSrcTar{i} = Tar_evecs(Tar_refined_boundaries{i},:,i)' * Tar_refined.SHAPE.S(Tar_refined_boundaries{i}, Tar_refined_boundaries{i}) * Src_evecs(fullp2pTarSrc(Tar_refined_boundaries{i}),:,i);
        FmapTarSrc{i} = Src_evecs(Src_refined_boundaries{i},:,i)' * Src_refined.SHAPE.S(Src_refined_boundaries{i}, Src_refined_boundaries{i}) * Tar_evecs(fullp2pSrcTar(Src_refined_boundaries{i}),:,i);
        
    end
        
    
    clf
    
    subplot(1,2,1)
    plot(Src_evecs(fullp2pTarSrc(Tar_refined_boundaries{1}),:,1)); title('Src evecs boundary on Tar'); colorbar
    subplot(1,2,2)
    plot(Tar_evecs(fullp2pSrcTar(Src_refined_boundaries{1}),:,1)); title('Tar evecs boundary on Src'); colorbar
    
    pause(1)
    
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
%     
%     pause(0.1)
    
    
    % Computing the new p2p map
    
    [SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p(Src_evecs, Src_evals, Tar_evecs, Tar_evals, FmapSrcTar, FmapTarSrc);
    [TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p(Tar_evecs, Tar_evals, Src_evecs, Src_evals, FmapTarSrc, FmapSrcTar);
    
        
    [Src_LB_TarSrc, Tar_LB_TarSrc] = build_LB_EmbeddingsForS2_to_S1_p2p(Src_LB_eF(:,1:zo), Tar_LB_eF(:,1:zo), FmapLBSrcTar, FmapLBTarSrc);
    [Tar_LB_SrcTar, Src_LB_SrcTar] = build_LB_EmbeddingsForS2_to_S1_p2p(Tar_LB_eF(:,1:zo), Src_LB_eF(:,1:zo), FmapLBTarSrc, FmapLBSrcTar);
    
    
%     fullp2pTarSrc = annquery([SrcFTarSrc, Src_LB_TarSrc]', [TarFTarSrc, Tar_LB_TarSrc]', 1);
%     fullp2pSrcTar = annquery([TarFSrcTar, Tar_LB_SrcTar]', [SrcFSrcTar, Src_LB_SrcTar]', 1);  
%     
    fullp2pTarSrc = annquery([SrcFTarSrc]', [TarFTarSrc]', 1);
    fullp2pSrcTar = annquery([TarFSrcTar]', [SrcFSrcTar]', 1);  

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
                  S1_LB * FmapLBS2S1; %Bijectivity term

        S2FS2S1 = S2_LB * FmapLBS1S2+... %"Proper" functional map term (pullback of p2p map)
                  S2_LB; %Bijectivity term

end












