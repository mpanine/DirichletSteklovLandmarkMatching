function [fullp2pTarSrc, fullp2pSrcTar] = CentralSteklovFindFullP2P_ZO_clean(Src_refined, Src_landmarks, Tar_refined, Tar_landmarks, FmapSrcTar, FmapTarSrc, Steklov_settings)
%findFullP2P_LB_V4: Converts the functional maps obtained at the boundary into point to point
%maps in the bulk. Uses both Dirichlet-Steklov eigenvunctions and Dirichlet LB eigenfunctions.
% This version works entirely in the bulk of the shape in order to
% guarantee an interaction between the DS and LB parts of the map.



%% Settings: make this into arguments, eventually

% USED_S_EVALS = 1:Steklov_settings.num_eigs; %TEMP SETTING.

toplot = true; % Set to true for some plots.

original_only_at_end = true; % If true, the final p2p maps are between the original meshes, including the known landmark correspondence.


%% Extract settings from Steklov_settings

SteklovWeight = sqrt(Steklov_settings.weight_Steklov);
ProperWeight = sqrt(Steklov_settings.weight_Proper);
OrthoWeight = sqrt(Steklov_settings.weight_Orthonormality);
BijectWeight = sqrt(Steklov_settings.weight_Bijectivity);

LBProperWeight = sqrt(Steklov_settings.weight_LB_Proper);
LBOrthoWeight = sqrt(Steklov_settings.weight_LB_Orthonormality);
LBBijectWeight = sqrt(Steklov_settings.weight_LB_Bijectivity);


ZO_start = Steklov_settings.LB_ZO_start;  %For now also used for the Steklov part
ZO_end = Steklov_settings.LB_ZO_end; 
ZO_step = Steklov_settings.LB_ZO_step;

%% Preliminaries

num_bounds = length(Src_refined.STEKLOV.boundaries);
% num_eigs = size(Src_refined.STEKLOV.evecs,2);


% Src_non_landmarks = 1:Src_refined.SHAPE.nv;
% Src_non_landmarks(Src_landmarks) = [];
% 
% Tar_non_landmarks = 1:Tar_refined.SHAPE.nv;
% Tar_non_landmarks(Tar_landmarks) = [];

% Src_non_Interior = 1:Src_refined.SHAPE.nv;
% Src_non_Interior(Src_refined.STEKLOV.Interior) = [];
% 
% Tar_non_Interior = 1:Tar_refined.SHAPE.nv;
% Tar_non_Interior(Tar_refined.STEKLOV.Interior) = [];


%% H1 - normalizing the DS eigenfunctions and correcting the Fmaps accordingly

% Src_sqrt_evals = zeros(size(Src_refined.STEKLOV.evals));
% Tar_sqrt_evals = zeros(size(Tar_refined.STEKLOV.evals));
% 
% Src_inv_sqrt_evals = zeros(size(Src_refined.STEKLOV.evals));
% Tar_inv_sqrt_evals = zeros(size(Tar_refined.STEKLOV.evals));
% 
% Src_evecs_normalized = Src_refined.STEKLOV.evecs;
% Tar_evecs_normalized = Tar_refined.STEKLOV.evecs;

% for i = 1:num_bounds %Normalizing the bulk DS eigenvectors and correcting the relevant Fmaps.
    
    Src_sqrt_evals = diag( diag( Src_refined.STEKLOV.evals ).^0.5 );
    Src_inv_sqrt_evals = diag( diag( Src_refined.STEKLOV.evals ).^-0.5 );
%     Src_inv_sqrt_evals(1,:) = []; %Remove the constant 
%     Src_inv_sqrt_evals(:,1) = [];
    
    Tar_sqrt_evals = diag( diag( Tar_refined.STEKLOV.evals ).^0.5 );
    Tar_inv_sqrt_evals = diag( diag( Tar_refined.STEKLOV.evals ).^-0.5 );
%     Src_inv_sqrt_evals(1,:) = [];
%     Src_inv_sqrt_evals(:,1) = [];
%     Tar_inv_sqrt_evals
    
    Src_evecs_normalized = Src_refined.STEKLOV.evecs * Src_inv_sqrt_evals; %DO NOT FORGET TO USE THESE IN THE CODE BELOW!!!
    Tar_evecs_normalized = Tar_refined.STEKLOV.evecs * Tar_inv_sqrt_evals;
    
    FmapSrcTar = Tar_sqrt_evals * FmapSrcTar * Src_inv_sqrt_evals;
    FmapTarSrc = Src_sqrt_evals * FmapTarSrc * Tar_inv_sqrt_evals;
    
% end

% figure
% subplot(1,3,1)
% imagesc(abs( Src_evecs_normalized' * Src_refined.SHAPE.W * Src_evecs_normalized ))
% title('Src evecs W self-product')
% colorbar
% 
% subplot(1,3,2)
% imagesc(abs( Src_refined.STEKLOV.LB_evecs' * Src_refined.SHAPE.W * Src_refined.STEKLOV.LB_evecs ))
% title('Src LB evecs W self-product')
% colorbar
% 
% subplot(1,3,3)
% imagesc(abs( Src_evecs_normalized' * Src_refined.SHAPE.W * Src_refined.STEKLOV.LB_evecs ))
% title('Src evecs vs Src LB evecs W product')
% colorbar
    
%% Computing the first guess of the Tar -> Src bulk point-to-point maps map.

[SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Src_evecs_normalized(:,1:ZO_start), Src_refined.STEKLOV.evals(1:ZO_start,1:ZO_start),...
                                                                   Tar_evecs_normalized(:,1:ZO_start), Tar_refined.STEKLOV.evals(1:ZO_start,1:ZO_start),...
                                                                   FmapSrcTar(1:ZO_start,1:ZO_start), FmapTarSrc(1:ZO_start,1:ZO_start),...
                                                                   SteklovWeight, ProperWeight, BijectWeight); %Function defined below.
                                                               
% fullp2pTarSrc = annquery(SrcFTarSrc(:,:,:)', TarFTarSrc(:,:,:)', 1);

% fullp2pTarSrc = SeparateInteriorExteriorNearestNeighbor(Src_refined, SrcFTarSrc, Tar_refined, TarFTarSrc);
fullp2pTarSrc = SeparateInteriorBoundaryLandmarkNearestNeighbor(Src_refined, SrcFTarSrc, Src_landmarks, Tar_refined, TarFTarSrc, Tar_landmarks);
                                                   

%% Computing the first guess of the Src -> Tar point-to-point maps map.
    
[TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Tar_evecs_normalized(:,1:ZO_start), Tar_refined.STEKLOV.evals(1:ZO_start,1:ZO_start),...
                                                                   Src_evecs_normalized(:,1:ZO_start), Src_refined.STEKLOV.evals(1:ZO_start,1:ZO_start),...
                                                                   FmapTarSrc(1:ZO_start,1:ZO_start), FmapSrcTar(1:ZO_start,1:ZO_start),...
                                                                   SteklovWeight, ProperWeight, BijectWeight);
% fullp2pSrcTar = annquery(TarFSrcTar(:,:,:)', SrcFSrcTar(:,:,:)', 1);  
% fullp2pSrcTar = SeparateInteriorExteriorNearestNeighbor(Tar_refined, TarFSrcTar, Src_refined, SrcFSrcTar);

fullp2pSrcTar = SeparateInteriorBoundaryLandmarkNearestNeighbor(Tar_refined, TarFSrcTar, Tar_landmarks, Src_refined, SrcFSrcTar, Src_landmarks);


%% Initial plot

if toplot
    
    figure('Units','normalized','Position',[0.05 0.05 0.85 0.85])
    
    subplot(2,6,1)
        imagesc(FmapSrcTar); colorbar;
        title('W-Steklov Fmap Src->Tar')

        subplot(2,6,2)
        imagesc(FmapTarSrc); colorbar;
        title('W-Steklov Fmap Tar->Src')

        subplot(2,6,3)
        imagesc(FmapTarSrc*FmapSrcTar); colorbar;
        title('W-Steklov Fmap Inv test')
        
        subplot(2,6,4)
        imagesc(zeros(size(FmapSrcTar))); colorbar;
        title('W-LB Fmap Src->Tar')

        subplot(2,6,5)
        imagesc(zeros(size(FmapSrcTar))); colorbar;
        title('W-LB Fmap Tar->Src')

        subplot(2,6,6)
        imagesc(zeros(size(FmapSrcTar))); colorbar;
        title('W-LB Fmap Inv test')
              

        subplot(2,6,7)
        plot(fullp2pTarSrc(fullp2pSrcTar),'.')        
        title('p2p Inv test Src -> Src')

        subplot(2,6,8)
        plot(fullp2pSrcTar(fullp2pTarSrc),'.')
        title('p2p Inv test Tar -> tar')
    
        
        subplot(2,6,10)
        plot( sqrt(sum(...
                        ( Src_refined.SHAPE.surface.VERT(fullp2pTarSrc(fullp2pSrcTar),:)...
                            - Src_refined.SHAPE.surface.VERT).^2 ...
                        ,2))      ,'.');
        ylim([0 1])
        title('p2p Inv Eucl.Err Src -> Src')

        subplot(2,6,11)
        plot( sqrt(sum(...
                        ( Tar_refined.SHAPE.surface.VERT(fullp2pSrcTar(fullp2pTarSrc),:)...
                            - Tar_refined.SHAPE.surface.VERT).^2 ...
                        ,2))      ,'.');
        ylim([0 1])
        title('p2p Inv Eucl.Err Tar -> tar')
        
    pause(2)
    
end


%% Main iteration

used_zo = ZO_start:ZO_step:ZO_end; % List of ZoomOut sizes iterated over
if used_zo(end) < ZO_end
    used_zo = [used_zo ZO_end]; % Ensures that the end size is always used.
end

for zo = used_zo
    
    %Computing the LB Fmap parts (bulk)
    
    FmapLBSrcTar = Tar_refined.STEKLOV.LB_evecs(:,1:zo)' * Tar_refined.SHAPE.W * Src_refined.STEKLOV.LB_evecs(fullp2pTarSrc,1:zo); %LB part of the map. H1 INNER PRIODUCT!!!
    FmapLBTarSrc = Src_refined.STEKLOV.LB_evecs(:,1:zo)' * Src_refined.SHAPE.W * Tar_refined.STEKLOV.LB_evecs(fullp2pSrcTar,1:zo);
    
    
    FmapSrcTar = zeros(zo,zo); %This Fmap is a sum of integals over boundary components.
    FmapTarSrc = zeros(zo,zo);
   
    for b = 1:num_bounds
        
        FmapSrcTar = FmapSrcTar +  Tar_refined.STEKLOV.evals(1:zo,1:zo) * Tar_evecs_normalized(Tar_refined.STEKLOV.boundaries{b},1:zo)'...
                        * Tar_refined.STEKLOV.MassBoundary{b} *...
                        Src_evecs_normalized(fullp2pTarSrc(Tar_refined.STEKLOV.boundaries{b}),1:zo);
                    
        FmapTarSrc = FmapTarSrc + Src_refined.STEKLOV.evals(1:zo,1:zo) * Src_evecs_normalized(Src_refined.STEKLOV.boundaries{b},1:zo)' ...
                        * Src_refined.STEKLOV.MassBoundary{b} * ...
                        Tar_evecs_normalized(fullp2pSrcTar(Src_refined.STEKLOV.boundaries{b}),1:zo);
        
    end
    

    
    %% Computing the new p2p map
    
    if (zo ~= ZO_end)||(~original_only_at_end) %Normal behavior for all iterations, except possibly the last one.
        
        [SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Src_evecs_normalized(:,1:zo), Src_refined.STEKLOV.evals(1:zo,1:zo),...
                                                                           Tar_evecs_normalized(:,1:zo), Tar_refined.STEKLOV.evals(1:zo,1:zo),...
                                                                           FmapSrcTar, FmapTarSrc, SteklovWeight, ProperWeight, BijectWeight);
        [TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Tar_evecs_normalized(:,1:zo), Tar_refined.STEKLOV.evals(1:zo,1:zo),...
                                                                           Src_evecs_normalized(:,1:zo), Src_refined.STEKLOV.evals(1:zo,1:zo),...
                                                                           FmapTarSrc, FmapSrcTar, SteklovWeight, ProperWeight, BijectWeight);


        [Src_LB_TarSrc, Tar_LB_TarSrc] = build_LB_EmbeddingsForS2_to_S1_p2p(Src_refined.STEKLOV.LB_evecs(:,1:zo), Tar_refined.STEKLOV.LB_evecs(:,1:zo),...
                                                                            FmapLBSrcTar, FmapLBTarSrc, LBProperWeight, LBOrthoWeight, LBBijectWeight);
        [Tar_LB_SrcTar, Src_LB_SrcTar] = build_LB_EmbeddingsForS2_to_S1_p2p(Tar_refined.STEKLOV.LB_evecs(:,1:zo), Src_refined.STEKLOV.LB_evecs(:,1:zo),...
                                                                            FmapLBTarSrc, FmapLBSrcTar, LBProperWeight, LBOrthoWeight, LBBijectWeight);

%         fullp2pTarSrc = annquery([SrcFTarSrc, Src_LB_TarSrc]', [TarFTarSrc, Tar_LB_TarSrc]', 1);
%         fullp2pSrcTar = annquery([TarFSrcTar, Tar_LB_SrcTar]', [SrcFSrcTar, Src_LB_SrcTar]', 1);  
        
        fullp2pTarSrc = SeparateInteriorBoundaryLandmarkNearestNeighbor(Src_refined, [SrcFTarSrc, Src_LB_TarSrc], Src_landmarks,...
                                                                        Tar_refined, [TarFTarSrc, Tar_LB_TarSrc], Tar_landmarks);
        fullp2pSrcTar = SeparateInteriorBoundaryLandmarkNearestNeighbor(Tar_refined, [TarFSrcTar, Tar_LB_SrcTar], Tar_landmarks,...
                                                                        Src_refined, [SrcFSrcTar, Src_LB_SrcTar], Src_landmarks);
        
    else %Compute the final p2p map on the interior. Activates when zo = ZO_end AND original_only_at_end is true (changes the last iteration)
        
        [SrcFTarSrc, TarFTarSrc] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Src_evecs_normalized(Src_refined.STEKLOV.Interior,1:zo), Src_refined.STEKLOV.evals(1:zo,1:zo),...
                                                                           Tar_evecs_normalized(Tar_refined.STEKLOV.Interior,1:zo), Tar_refined.STEKLOV.evals(1:zo,1:zo),...
                                                                           FmapSrcTar, FmapTarSrc, SteklovWeight, ProperWeight, BijectWeight);
        [TarFSrcTar, SrcFSrcTar] = build_DS_EmbeddingsForS2_to_S1_p2p_bulk(Tar_evecs_normalized(Tar_refined.STEKLOV.Interior,1:zo), Tar_refined.STEKLOV.evals(1:zo,1:zo),...
                                                                           Src_evecs_normalized(Src_refined.STEKLOV.Interior,1:zo), Src_refined.STEKLOV.evals(1:zo,1:zo),...
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
    
    
    if toplot
    
        subplot(2,6,1)
        imagesc(FmapSrcTar); colorbar;
        title('W-Steklov Fmap Src->Tar')

        subplot(2,6,2)
        imagesc(FmapTarSrc); colorbar;
        title('W-Steklov Fmap Tar->Src')

        subplot(2,6,3)
        imagesc(FmapTarSrc*FmapSrcTar); colorbar;
        title('W-Steklov Fmap Inv test')
        
        subplot(2,6,4)
        imagesc(FmapLBSrcTar); colorbar;
        title('W-LB Fmap Src->Tar')

        subplot(2,6,5)
        imagesc(FmapLBTarSrc); colorbar;
        title('W-LB Fmap Tar->Src')

        subplot(2,6,6)
        imagesc(FmapLBTarSrc*FmapLBSrcTar); colorbar;
        title('W-LB Fmap Inv test')
              

        
        subplot(2,6,7)
        plot(fullp2pTarSrc(fullp2pSrcTar),'.')        
        title('p2p Inv test Src -> Src')

        subplot(2,6,8)
        plot(fullp2pSrcTar(fullp2pTarSrc),'.')
        title('p2p Inv test Tar -> tar')
    
        
        if (zo ~= ZO_end)||(~original_only_at_end)
        
            subplot(2,6,10)
            plot( sqrt(sum(...
                            ( Src_refined.SHAPE.surface.VERT(fullp2pTarSrc(fullp2pSrcTar),:)...
                                - Src_refined.SHAPE.surface.VERT).^2 ...
                            ,2))      ,'.');
            ylim([0 1])
            title('p2p Inv Eucl.Err Src -> Src')

            subplot(2,6,11)
            plot( sqrt(sum(...
                            ( Tar_refined.SHAPE.surface.VERT(fullp2pSrcTar(fullp2pTarSrc),:)...
                                - Tar_refined.SHAPE.surface.VERT).^2 ...
                            ,2))      ,'.');
            ylim([0 1])
            title('p2p Inv Eucl.Err Tar -> tar')

        else
            
            subplot(2,6,10)
            plot( sqrt(sum(...
                            ( Src_refined.SHAPE.surface.VERT(fullp2pTarSrc(fullp2pSrcTar),:)...
                                - Src_refined.SHAPE.surface.VERT(1:length(fullp2pSrcTar),:) ).^2 ...
                            ,2))      ,'.');
            ylim([0 1])
            title('p2p Inv Eucl.Err Src -> Src')

            subplot(2,6,11)
            plot( sqrt(sum(...
                            ( Tar_refined.SHAPE.surface.VERT(fullp2pSrcTar(fullp2pTarSrc),:)...
                                - Tar_refined.SHAPE.surface.VERT(1:length(fullp2pTarSrc),:)).^2 ...
                            ,2))      ,'.');
            ylim([0 1])
            title('p2p Inv Eucl.Err Tar -> tar')
            
            
        end
        
        
        
        pause(0.00001)
    
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

        S1FS2S1(:,start_ind:end_ind) = SteklovWeight * S1_evecs*S1_evals +... %Steklov term
               ProperWeight * S1_evecs+... %"Proper" functional map term (pullback of p2p map)
               BijectWeight * S1_evecs*FmapS2S1; %Bijectivity term

        S2FS2S1(:,start_ind:end_ind) = SteklovWeight * S2_evecs*S2_evals*FmapS1S2 +... %Steklov term
               ProperWeight * S2_evecs*FmapS1S2+... %"Proper" functional map term (pullback of p2p map)
               BijectWeight * S2_evecs; %Bijectivity term


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












