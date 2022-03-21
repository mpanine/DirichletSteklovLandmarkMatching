% Tests the Landmark Dirichlet-Steklov Approach


%% Technical stuff

clear all;
close all;

addpath(genpath(pwd))

shape_settings = shape.Settings('computeMHB', true,...
                        'computeGeoDist', false,...
                        'MHB_settings', basis.MHB_Settings(20),...
                        'computeNormals', false);
cache_settings = cache.Settings('a','b', false);




%% Load source
Src = struct;
% 
% Src.SHAPE = shape.compute('./Data/bunny_conformal_short/sphere_10242.off', ...
%                         shape_settings, cache_settings);
% load('./Data/bunny_conformal_short/Landmarks/sphere_10242')

% Src.SHAPE = shape.compute('./Data/bunny_conformal_short/flow00000000.off', ...
%                         shape_settings, cache_settings);
% load('./Data/bunny_conformal_short/Landmarks/flow00000000')

Src.SHAPE = shape.compute('./Data/Shrec07_fourleg/381.off', ...
                        shape_settings, cache_settings);
load('./Data/Shrec07_fourleg/Landmarks/381') 

% Src.SHAPE = shape.compute('./Data/Shrec07_fourleg/382.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Shrec07_fourleg/Landmarks/382') 

% Src.SHAPE = shape.compute('./Data/Faust_remeshed/vtx_5k/tr_reg_000.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Faust_remeshed/Landmarks/tr_reg_000') 


Src_landmarks = landmarks(1:21);
% Src_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13]);      
% Src_landmarks = landmarks([1, 6, 8, 10, 12, 14]); 
% Src_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13, 21, 15]); 
% Src_landmarks = landmarks(1:20)
% Src_landmarks = landmarks([1, 6, 8]); 
                    
% Src_basis = ComputeHarmonicBasis(Src.SHAPE, Src_landmarks);

% [~, Src_predom] = sort(Src_basis,2,'descend');

% Src_dom = Src_predom(:,1);

%% Load target
Tar = struct;
% 
% Tar.SHAPE = shape.compute('./Data/bunny_conformal_short/flow00000217.off', ...
%                         shape_settings, cache_settings);
% load('./Data/bunny_conformal_short/Landmarks/flow00000217') 

% Tar.SHAPE = shape.compute('./Data/bunny_conformal_short/flow00000000.off', ...
%                         shape_settings, cache_settings);
% load('./Data/bunny_conformal_short/Landmarks/flow00000000') 

% Tar.SHAPE = shape.compute('./Data/bunny_conformal_short/sphere_10242.off', ...
%                         shape_settings, cache_settings);
% load('./Data/bunny_conformal_short/Landmarks/sphere_10242')


Tar.SHAPE = shape.compute('./Data/Shrec07_fourleg/388.off', ...
                        shape_settings, cache_settings);
load('./Data/Shrec07_fourleg/Landmarks/388')   

% Tar.SHAPE = shape.compute('./Data/Shrec07_fourleg/400.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Shrec07_fourleg/Landmarks/400')   

% 
% Tar.SHAPE = shape.compute('./Data/Faust_remeshed/vtx_5k/tr_reg_001.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Faust_remeshed/Landmarks/tr_reg_001') 


Tar_landmarks = landmarks(1:21);
% Tar_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13]);  
% Tar_landmarks = landmarks([1, 6, 8, 10, 12, 14]); 
% Tar_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13, 21, 15]); 
% Tar_landmarks = landmarks(1:20);
% Tar_landmarks = landmarks([1, 6, 8]); 
                    
% Tar_basis = ComputeHarmonicBasis(Tar.SHAPE, Tar_landmarks);

% [~, Tar_predom] = sort(Tar_basis,2,'descend');

% Tar_dom = Tar_predom(:,1);



%% Internal Settings

num_eigs = 20;

numsubdivisions = num_eigs*ones(size(Src_landmarks)); % 50*ones(size(Src_landmarks)); %Number of subdivisions per triangle at landmark.
radii_factor = 0.5;




%% Landmark radii 

Src_shortest_edges = findShortestEdgesAtLandmarks(Src, Src_landmarks);
Tar_shortest_edges = findShortestEdgesAtLandmarks(Tar, Tar_landmarks);
    
        landmarks_radii = zeros(size(Src_shortest_edges));
        for land_ind = 1:length(landmarks_radii)
            
            landmarks_radii(land_ind) = min(Src_shortest_edges(land_ind), Tar_shortest_edges(land_ind) );
            
        end
    
        landmarks_radii = radii_factor * landmarks_radii;




%% Compute Bases


[Src_refined, Src_segment_TRIV, Src_evecs, Src_evals, Src_boundary_edges_refined, Src_boundaries_new, Src_Mass_Boundary, Src_Interior, Src_refined_boundaries, Src_normalTriangles, Src_W_segment] = ComputeCentralDirichletSteklovBasis(Src, Src_landmarks, landmarks_radii, numsubdivisions, num_eigs);

[Tar_refined, Tar_segment_TRIV, Tar_evecs, Tar_evals, Tar_boundary_edges_refined, Tar_boundaries_new, Tar_Mass_Boundary, Tar_Interior, Tar_refined_boundaries, Tar_normalTriangles, Tar_W_segment] = ComputeCentralDirichletSteklovBasis(Tar, Tar_landmarks, landmarks_radii, numsubdivisions, num_eigs);

Src_refined_boundary_list = unique(Src_boundary_edges_refined);
Tar_refined_boundary_list = unique(Tar_boundary_edges_refined);


%% Plot W inner products between Dirichlet-Steklov functions

aaa = 2;
bbb = 5;

figure
imagesc( abs( (Src_evecs(:,:,aaa) * diag( diag(Src_evals(:,:,aaa)).^-.5  ))' * Src_W_segment * ( Src_evecs(:,:,bbb) * diag( diag(Src_evals(:,:,bbb)).^-.5  ) ) ) )
colorbar

%% Test Interior and boundary
% figure
% subplot(1,2,1)
% trimesh(Src_refined.SHAPE.surface.TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'k'); axis equal;
% hold on
% plot3(Src_refined.SHAPE.surface.X(Src_Interior),Src_refined.SHAPE.surface.Y(Src_Interior),Src_refined.SHAPE.surface.Z(Src_Interior),'.','MarkerSize',20,'color','r')
% plot3(Src_refined.SHAPE.surface.X(Src_boundary_list),Src_refined.SHAPE.surface.Y(Src_boundary_list),Src_refined.SHAPE.surface.Z(Src_boundary_list),'.','MarkerSize',20,'color','b')
% title('Interior and boundary test')
% 
% subplot(1,2,2)
% trimesh(Tar_refined.SHAPE.surface.TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'k'); axis equal;
% hold on
% plot3(Tar_refined.SHAPE.surface.X(Tar_Interior),Tar_refined.SHAPE.surface.Y(Tar_Interior), Tar_refined.SHAPE.surface.Z(Tar_Interior),'.','MarkerSize',20,'color','r')
% plot3(Tar_refined.SHAPE.surface.X(Tar_boundary_list),Tar_refined.SHAPE.surface.Y(Tar_boundary_list), Tar_refined.SHAPE.surface.Z(Tar_boundary_list),'.','MarkerSize',20,'color','b')
% title('Interior and boundary test')
% 



%% Normalize spectra

for i = 1:length(Src_landmarks)
    
    Src_evals(:,:,i) = Src_evals(:,:,i)/max(max(Src_evals(:,:,i)));
    Tar_evals(:,:,i) = Tar_evals(:,:,i)/max(max(Tar_evals(:,:,i)));
    
end

% for i = 1:length(Src_landmarks)
%     
%     Src_evals(:,:,i) = Src_evals(:,:,i)/(Src_evals(10,10,i));
%     Tar_evals(:,:,i) = Tar_evals(:,:,i)/(Tar_evals(10,10,i));
%     
% end

%% Extend to full refined mesh by padding with NaNs   --- THIS IS REPLACED BELOW

Src_nv = size(Src_refined.SHAPE.surface.VERT,1);
Src_segment = 1:Src_nv;
Src_segment(Src_landmarks) = [];

Src_eF = extendSegmentFunctions(Src_evecs, Src_segment, Src_nv);


Tar_nv = size(Tar_refined.SHAPE.surface.VERT,1);
Tar_segment = 1:Tar_nv;
Tar_segment(Tar_landmarks) = [];

Tar_eF = extendSegmentFunctions(Tar_evecs, Tar_segment, Tar_nv);






%% Plot eigenfunctions



%     nef = 1;
%     n_bound = 1;
% 
%     figure
%     subplot(1,2,1)
%     trimesh(Src_segment_TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, abs(Src_eF(:,nef,n_bound))  , ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
%     title(sprintf('Lmk %d, DS eigenfunction %d',n_bound,nef))
% 
%     subplot(1,2,2)
%     trimesh(Tar_segment_TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, abs(Tar_eF(:,nef,n_bound) ), ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
%     title(sprintf('Lmk %d, DS eigenfunction %d',n_bound,nef))




%% Plotting functions on boundary components

% b_component = 1;
% ef_order = 4;
% landmark_num = 1;
% 
% figure
% subplot(1,2,1)
% plot(Src_evecs(Src_boundaries_new{b_component}, ef_order ,landmark_num))
% title('Source')
% 
% subplot(1,2,2)
% plot(Tar_evecs(Tar_boundaries_new{b_component}, ef_order ,landmark_num))
% title('Target')


%% Bases restricted to the boundary components.

for i = 1:length(Src_landmarks) % Both bases restricted to the boundary.
    
    Src_evecs_bound{i} = Src_evecs(Src_boundaries_new{i}, : ,i);
    Tar_evecs_bound{i} = Tar_evecs(Tar_boundaries_new{i}, : ,i);
    

    
%     figure
%     imagesc(log(abs(Fmap_guess{i}) + 1e-12) )
%     imagesc( Fmap_guess{i} )
    
end


%% Replace the NaNs by integrals at the boundaries

Src_eF(Src_landmarks,:,:) = 0; % Assigns all landmark values to 0.
Tar_eF(Tar_landmarks,:,:) = 0;

for i = 1:length(Src_evecs_bound) %Cycle over the boundaries with nonzero integrals
    
    Src_eF(Src_landmarks(i),:,i) = sum(Src_Mass_Boundary{i}*Src_evecs_bound{i},1); % Values at the landmarks are integrals over the boundary.
    Tar_eF(Tar_landmarks(i),:,i) = sum(Tar_Mass_Boundary{i}*Tar_evecs_bound{i},1);

end


%% Plot the extended eigenfunctions

% 
% 
%     nef = 1;
%     n_bound = 5;
% 
%     figure
%     subplot(1,2,1)
%     trimesh(Src_refined.SHAPE.surface.TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, Src_eF(:,nef,n_bound)  , ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'k'); axis equal; axis off; colorbar;
%     title(sprintf('Lmk %d, DS eigenfunction %d',n_bound,nef))
% 
%     subplot(1,2,2)
%     trimesh(Tar_refined.SHAPE.surface.TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, Tar_eF(:,nef,n_bound) , ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'k'); axis equal; axis off; colorbar;
%     title(sprintf('Lmk %d, DS eigenfunction %d',n_bound,nef))






%% Computing the correspondence

% for i = 1:length(Src_evecs_bound) %Guess the initial Fmap
% 
%     Fmap_guess{i} = Tar_guess{i}' * Src_guess{i};
% 
% end
% 
% p2p = DirichletSteklovFmapToP2P(Src_evecs_bound, Src_evals, Tar_evecs_bound, Tar_evals, Fmap_guess);



% [Fmap, p2p] = DirichletSteklovZoomOutV1(Src_evecs_bound, Src_evals, Src_Mass_Boundary, Tar_evecs_bound, Tar_evals, Tar_Mass_Boundary);





%% Computing the boundary correspondence - Bijective version

% [FmapSrcTar, p2pTarSrc, FmapTarSrc, p2pSrcTar] = DirichletSteklovZoomOutV2(Src_evecs_bound, Src_evals, Src_Mass_Boundary, Tar_evecs_bound, Tar_evals, Tar_Mass_Boundary);


%% Computing the boundary correspondence - Exhaustive initial search
% clc
% [FmapSrcTar, p2pTarSrc, FmapTarSrc, p2pSrcTar] = DirichletSteklovZoomOutV3(Src_evecs_bound, Src_evals, Src_Mass_Boundary, Tar_evecs_bound, Tar_evals, Tar_Mass_Boundary);

%% Computing the boundary correspondence - Exhaustive initial search wrt 

Src_harmonic_basis_refined = ComputeRefinedCentralHarmonicBasis(Src_refined, Src_refined_boundaries, Src_landmarks);
Tar_harmonic_basis_refined = ComputeRefinedCentralHarmonicBasis(Tar_refined, Tar_refined_boundaries, Tar_landmarks);

% [FmapSrcTar, p2pTarSrc, FmapTarSrc, p2pSrcTar] = DirichletSteklovZoomOutV4(Src, Src_refined, Src_evecs_bound, Src_evals, Src_Mass_Boundary, Src_landmarks, Src_harmonic_basis_refined, Src_refined_boundaries, Tar, Tar_refined, Tar_evecs_bound, Tar_evals, Tar_Mass_Boundary, Tar_landmarks, Tar_harmonic_basis_refined, Tar_refined_boundaries);

[FmapSrcTar, p2pTarSrc, FmapTarSrc, p2pSrcTar] = DirichletSteklovZoomOutV4(Src_refined, Src_evecs_bound, Src_evals, Src_Mass_Boundary, Src_refined_boundaries, Src_normalTriangles, Src_harmonic_basis_refined, Tar_refined, Tar_evecs_bound, Tar_evals, Tar_Mass_Boundary, Tar_refined_boundaries, Tar_normalTriangles, Tar_harmonic_basis_refined);



%% Find the full p2p map -- REFINED SHAPES
% 
% [fullp2pTarSrc, fullp2pSrcTar] = findFullP2P(Src_eF, Src_evals, Tar_eF, Tar_evals, FmapSrcTar, FmapTarSrc);

%% Find the full p2p map -- ORIGINAL SHAPES

[fullp2pTarSrc, fullp2pSrcTar] = findFullP2P(Src_eF(1:Src.SHAPE.nv,:,:), Src_evals, Tar_eF(1:Tar.SHAPE.nv,:,:), Tar_evals, FmapSrcTar, FmapTarSrc);

%% Find the full p2p map -- refined shapes, LB included

[Src_LB_eF, Src_LB_evals] = computeDirichletLB(Src_refined.SHAPE, [Src_landmarks(:); Src_refined_boundary_list(:)], 200);
[Tar_LB_eF, Tar_LB_evals] = computeDirichletLB(Tar_refined.SHAPE, [Tar_landmarks(:); Tar_refined_boundary_list(:)], 200);

Src_LB_eF_interior = Src_LB_eF(1:Src.SHAPE.nv,:); % Used to compare to HD ZoomOut (unpublished paper)
Tar_LB_eF_interior = Tar_LB_eF(1:Tar.SHAPE.nv,:);

Src_LB_eF_unmod = Src_LB_eF;
Tar_LB_eF_unmod = Tar_LB_eF;

Src_LB_eF = Src_LB_eF * diag( Src_LB_evals.^-0.5 );
Tar_LB_eF = Tar_LB_eF * diag( Tar_LB_evals.^-0.5 );

%% Find the full p2p map -- refined shapes, LB included part 2
% clc
% [fullp2pTarSrc_test, fullp2pSrcTar_test] = findFullP2P_LB_V1(Src_eF, Src_evals, Src_LB_eF, Src_refined_boundaries, Src_refined,...
%                                                              Tar_eF, Tar_evals, Tar_LB_eF, Tar_refined_boundaries, Tar_refined,...
%                                                              FmapSrcTar, FmapTarSrc);
% 
%                                                          
%% Find the full p2p map -- refined shapes, LB included part 2 -- INTERIOR - EXTERIOR SPLIT -- the output only works in the interior
% 
% [fullp2pTarSrc_test, fullp2pSrcTar_test] = findFullP2P_LB_V2(Src_eF, Src_evals, Src_landmarks, Src_LB_eF, Src_refined_boundaries, Src, Src_refined, Src_Interior, Src_Mass_Boundary,...
%                                                              Tar_eF, Tar_evals, Tar_landmarks, Tar_LB_eF, Tar_refined_boundaries, Tar, Tar_refined, Tar_Interior, Tar_Mass_Boundary,...
%                                                              FmapSrcTar, FmapTarSrc);
% 
% %                                                          
% %% Find the full p2p map -- refined shapes, LB included part 2 -- INTERIOR - EXTERIOR SPLIT -- the output only works in the interior -- V3 matching
% clc
% [fullp2pTarSrc_test, fullp2pSrcTar_test] = findFullP2P_LB_V3(Src_eF, Src_evals, Src_landmarks, Src_LB_eF, Src_refined_boundaries, Src, Src_refined, Src_Interior, Src_Mass_Boundary,...
%                                                              Tar_eF, Tar_evals, Tar_landmarks, Tar_LB_eF, Tar_refined_boundaries, Tar, Tar_refined, Tar_Interior, Tar_Mass_Boundary,...
%                                                              FmapSrcTar, FmapTarSrc);
                                                         
%% Find the full p2p map -- refined shapes, LB included part 4 -- INTERIOR - EXTERIOR SPLIT -- the output only works in the interior -- V4 matching
% clc
[fullp2pTarSrc_test, fullp2pSrcTar_test] = findFullP2P_LB_V4(Src_eF, Src_evals, Src_landmarks, Src_LB_eF, Src_refined_boundaries, Src, Src_refined, Src_Interior, Src_Mass_Boundary,...
                                                             Tar_eF, Tar_evals, Tar_landmarks, Tar_LB_eF, Tar_refined_boundaries, Tar, Tar_refined, Tar_Interior, Tar_Mass_Boundary,...
                                                             FmapSrcTar, FmapTarSrc);
                                                                                                             


%% Map to the original shape

fullp2pSrcTar_test = reinsertBoundaryCorrespondence(fullp2pSrcTar_test, Src_landmarks, Tar_landmarks, Tar_Interior);                                                         
fullp2pTarSrc_test = reinsertBoundaryCorrespondence(fullp2pTarSrc_test, Tar_landmarks, Src_landmarks, Src_Interior);                                                            


%% Find the full p2p map -- ORIGINAL SHAPES

% [fullp2pTarSrc_test, fullp2pSrcTar_test] = findFullP2P_LB_V1(Src_eF(1:Src.SHAPE.nv,:,:), Src_evals, Tar_eF(1:Tar.SHAPE.nv,:,:), Tar_evals, FmapSrcTar, FmapTarSrc);

%% Find full p2p map -- Original shapes, EXPEIMENTAL NORMALIZATION-- DOES NOT WORK

% [fullp2pTarSrc_normalized, fullp2pSrcTar_normalized] = findFullP2Pnormalized(Src_eF(1:Src.SHAPE.nv,:,:), Src_evals, Tar_eF(1:Tar.SHAPE.nv,:,:), Tar_evals, FmapSrcTar, FmapTarSrc);

%% Find the full p2p map -- ORIGINAL SHAPES, Iterative, DOES NOT WORK!

% IDEA: ENFORCE BIJECTIVITY IN LB BASIS (not implemented)

% [fullp2pTarSrc_ZO, fullp2pSrcTar_ZO, FmapSrcTar_combined, FmapTarSrc_combined] = findFullP2P_Iterative(Src_eF(1:Src.SHAPE.nv,:,:), Src_evals, Src.SHAPE.A, Tar_eF(1:Tar.SHAPE.nv,:,:), Tar_evals, Tar.SHAPE.A, FmapSrcTar, FmapTarSrc);
% 

%% Plot coord transfer -- REFINED SHAPES

% figure
% subplot(1,2,1)
% plotCoordinateTransfer(Tar_refined.SHAPE, Src_refined.SHAPE, fullp2pSrcTar)
% title('Coordinate transfer Tar->Src')
% 
% subplot(1,2,2)
% plotCoordinateTransfer(Src_refined.SHAPE, Tar_refined.SHAPE, fullp2pTarSrc)
% title('Coordinate transfer Src->Tar')

%% Plot coord transfer -- ORIGINAL SHAPES

figure
subplot(1,2,1)
plotCoordinateTransfer(Tar.SHAPE, Src.SHAPE, fullp2pSrcTar)
hold on
plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
title('Steklov Coordinate transfer Tar->Src')

subplot(1,2,2)
plotCoordinateTransfer(Src.SHAPE, Tar.SHAPE, fullp2pTarSrc)
hold on
plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
title('Steklov Coordinate transfer Src->Tar')


%% Plot coord transfer -- ZoomOut, refined shape

% figure
% subplot(1,2,1)
% plotCoordinateTransfer(Tar_refined.SHAPE, Src_refined.SHAPE, fullp2pSrcTar_test)
% hold on
% plot3(Tar_refined.SHAPE.surface.X(Tar_landmarks),Tar_refined.SHAPE.surface.Y(Tar_landmarks),Tar_refined.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
% title('ZO, refined, Steklov Coordinate transfer Tar->Src')
% 
% subplot(1,2,2)
% plotCoordinateTransfer(Src_refined.SHAPE, Tar_refined.SHAPE, fullp2pTarSrc_test)
% hold on
% plot3(Src_refined.SHAPE.surface.X(Src_landmarks),Src_refined.SHAPE.surface.Y(Src_landmarks),Src_refined.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
% title('ZO, refined, Steklov Coordinate transfer Src->Tar')

%% Plot coord transfer -- ZoomOut, original shape

figure
subplot(1,2,1)
plotCoordinateTransfer(Tar.SHAPE, Src.SHAPE, fullp2pSrcTar_test)
hold on
plot3(Tar_refined.SHAPE.surface.X(Tar_landmarks),Tar_refined.SHAPE.surface.Y(Tar_landmarks),Tar_refined.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
title('ZO, original, Steklov Coordinate transfer Tar->Src')

subplot(1,2,2)
plotCoordinateTransfer(Src.SHAPE, Tar.SHAPE, fullp2pTarSrc_test)
hold on
plot3(Src_refined.SHAPE.surface.X(Src_landmarks),Src_refined.SHAPE.surface.Y(Src_landmarks),Src_refined.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
title('ZO, original, Steklov Coordinate transfer Src->Tar')

%% Plot map lines -- ORIGINAL SHAPES

figure
subplot(1,2,1)
visualize_map_lines(Src.SHAPE.surface, Tar.SHAPE.surface, fullp2pSrcTar_test, Src_landmarks)
% hold on
% plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
title('Steklov ZO Function transfer Tar->Src')

subplot(1,2,2)
visualize_map_lines(Tar.SHAPE.surface, Src.SHAPE.surface, fullp2pTarSrc_test, Tar_landmarks)
% hold on
% plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
title('Steklov ZO Function transfer Src->Tar')

%% Plot map lines -- REFINED SHAPES

% figure
% subplot(1,2,1)
% visualize_map_lines(Src_refined.SHAPE.surface, Tar_refined.SHAPE.surface, fullp2pSrcTar_test, Src_landmarks)
% % hold on
% % plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
% title('Steklov ZO Function transfer Tar->Src')
% 
% subplot(1,2,2)
% visualize_map_lines(Tar_refined.SHAPE.surface, Src_refined.SHAPE.surface, fullp2pTarSrc_test, Tar_landmarks)
% % hold on
% % plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
% title('Steklov ZO Function transfer Src->Tar')



%% Compare to central harmonic matching

% Src_harmonic_basis = ComputeCentralHarmonicBasis(Src.SHAPE, Src_landmarks, landmarks_radii, numsubdivisions);
% Tar_harmonic_basis = ComputeCentralHarmonicBasis(Tar.SHAPE, Tar_landmarks, landmarks_radii, numsubdivisions);
% 
% pF_Src_to_Tar = annquery( Tar_harmonic_basis', Src_harmonic_basis', 1);
% pF_Tar_to_Src = annquery( Src_harmonic_basis', Tar_harmonic_basis', 1);

%% Compare to central harmonic matching

Src_harmonic_basis = ComputeCentralHarmonicBasis(Src.SHAPE, Src_landmarks, landmarks_radii, numsubdivisions);
Tar_harmonic_basis = ComputeCentralHarmonicBasis(Tar.SHAPE, Tar_landmarks, landmarks_radii, numsubdivisions);

% pF_Src_to_Tar = annquery( Tar_harmonic_basis_refined', Src_harmonic_basis_refined', 1);
% pF_Tar_to_Src = annquery( Src_harmonic_basis_refined', Tar_harmonic_basis_refined', 1);


pF_Src_to_Tar = annquery( Tar_harmonic_basis', Src_harmonic_basis', 1);
pF_Tar_to_Src = annquery( Src_harmonic_basis', Tar_harmonic_basis', 1);

%% Compare to central harmonic plus HD-ZoomOut (unpublished paper)
% 
[pF_Src_to_Tar_HDZO, pF_Tar_to_Src_HDZO] = isometricZoomOut(Src.SHAPE, Tar.SHAPE, Src_harmonic_basis, Tar_harmonic_basis, Src_LB_eF_interior, Tar_LB_eF_interior,...
                                                pF_Src_to_Tar, pF_Tar_to_Src, 20, 100, 5, Src_landmarks, Tar_landmarks);

                                            
%% Compare to central harmonic plus HD-ZoomOut (unpublished paper) -- REFINED SHAPES

% [pF_Src_to_Tar_HDZO, pF_Tar_to_Src_HDZO] = isometricZoomOut(Src_refined.SHAPE, Tar_refined.SHAPE, Src_harmonic_basis_refined, Tar_harmonic_basis_refined, Src_LB_eF_unmod, Tar_LB_eF_unmod,...
%                                                 pF_Src_to_Tar, pF_Tar_to_Src, 20, 100, 5, Src_landmarks, Tar_landmarks);



%% Harmonic coordinate transfer
% figure
% subplot(1,2,1)
% plotCoordinateTransfer(Tar_refined.SHAPE, Src_refined.SHAPE, pF_Src_to_Tar)
% hold on
% plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
% title('Harmonic Coordinate transfer Tar->Src')
% 
% subplot(1,2,2)
% plotCoordinateTransfer(Src_refined.SHAPE, Tar_refined.SHAPE, pF_Tar_to_Src)
% hold on
% plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
% title('Harmonic Coordinate transfer Src->Tar')
% 
% 
% 
%% Harmonic + HD ZoomOut coordinate transfer
% figure
% subplot(1,2,1)
% plotCoordinateTransfer(Tar_refined.SHAPE, Src_refined.SHAPE, pF_Src_to_Tar_HDZO)
% hold on
% plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
% title('Harmonic + HD ZO Coordinate transfer Tar->Src')
% 
% subplot(1,2,2)
% plotCoordinateTransfer(Src_refined.SHAPE, Tar_refined.SHAPE, pF_Tar_to_Src_HDZO)
% hold on
% plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
% title('Harmonic + HD ZO Coordinate transfer Src->Tar')


%% Harmonic coordinate transfer
figure
subplot(1,2,1)
plotCoordinateTransfer(Tar.SHAPE, Src.SHAPE, pF_Src_to_Tar)
hold on
plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
title('Harmonic Coordinate transfer Tar->Src')

subplot(1,2,2)
plotCoordinateTransfer(Src.SHAPE, Tar.SHAPE, pF_Tar_to_Src)
hold on
plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
title('Harmonic Coordinate transfer Src->Tar')



%% Harmonic + HD ZoomOut coordinate transfer
figure
subplot(1,2,1)
plotCoordinateTransfer(Tar.SHAPE, Src.SHAPE, pF_Src_to_Tar_HDZO)
hold on
plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
title('Harmonic + HD ZO Coordinate transfer Tar->Src')

subplot(1,2,2)
plotCoordinateTransfer(Src.SHAPE, Tar.SHAPE, pF_Tar_to_Src_HDZO)
hold on
plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
title('Harmonic + HD ZO Coordinate transfer Src->Tar')






%% Plot harmonic map lines -- ORIGINAL SHAPES

figure
subplot(1,2,1)
visualize_map_lines(Src.SHAPE.surface, Tar.SHAPE.surface, pF_Src_to_Tar_HDZO, Src_landmarks)
% hold on
% plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
title('Harmonic + HDZO Function transfer Tar->Src')

subplot(1,2,2)
visualize_map_lines(Tar.SHAPE.surface, Src.SHAPE.surface, pF_Tar_to_Src_HDZO, Tar_landmarks)
% hold on
% plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
title('Harmonic + HDZO Function transfer Src->Tar')



%% BIZZARE STUFF STARTS HERE

% 
% bizzare = false;
% 
% if bizzare
%     
%     
%     
% %% Plot the divided extended eigenfunctions
% 
% 
% 
%     nef = 2;
%     n_bound = 16;
% 
%     figure
%     subplot(1,2,1)
%     trimesh(Src_refined.SHAPE.surface.TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, Src_eF(:,nef,n_bound)./(Src_eF(:,1,n_bound)+1e-9)  , ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
%     title(sprintf('Lmk %d, DS divided eigenfunction %d',n_bound,nef))
% 
%     subplot(1,2,2)
%     trimesh(Tar_refined.SHAPE.surface.TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, Tar_eF(:,nef,n_bound)./(Tar_eF(:,1,n_bound)+1e-9) , ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
%     title(sprintf('Lmk %d, DS divided eigenfunction %d',n_bound,nef))
% 
%     
%     figure
%     subplot(1,2,1)
%     trimesh(Src_refined.SHAPE.surface.TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, Src_refined.SHAPE.W*(Src_eF(:,nef,n_bound)./(Src_eF(:,1,n_bound)+1e-9))  , ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
%     title(sprintf('Lmk %d, Lap of DS divided eigenfunction %d',n_bound,nef))
% 
%     subplot(1,2,2)
%     trimesh(Tar_refined.SHAPE.surface.TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, Tar_refined.SHAPE.W*(Tar_eF(:,nef,n_bound)./(Tar_eF(:,1,n_bound)+1e-9)) , ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
%     title(sprintf('Lmk %d, Lap of DS divided  eigenfunction %d',n_bound,nef))
% 
% 
%     
% end
% 


if false
    
    %% Plot the original shapes.
    
    figure
    subplot(1,2,1)
    trimesh(Src.SHAPE.surface.TRIV, Src.SHAPE.surface.X, Src.SHAPE.surface.Y, Src.SHAPE.surface.Z  , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'k'); axis equal; axis off; %colorbar;
    hold on
    plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
    title('Src')

    subplot(1,2,2)
    trimesh(Tar.SHAPE.surface.TRIV, Tar.SHAPE.surface.X, Tar.SHAPE.surface.Y, Tar.SHAPE.surface.Z , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'k'); axis equal; axis off; %colorbar;
    hold on
    plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')

    title('Tar')
    
    
    
    
    
    
end







if false
    
    
    
    %% Plot sign-log-abs eigenfunctions



    nef = 1;
    n_bound = 1;

    figure
    subplot(1,2,1)
    trimesh(Src_refined.SHAPE.surface.TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, sign(Src_eF(:,nef,n_bound)).*log(abs(Src_eF(:,nef,n_bound)))  , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
    title(sprintf('Lmk %d, sign-log-abs DS eigenfunction %d',n_bound,nef))

    subplot(1,2,2)
    trimesh(Tar_refined.SHAPE.surface.TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, sign(Tar_eF(:,nef,n_bound)).*log(abs(Tar_eF(:,nef,n_bound))) , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
    title(sprintf('Lmk %d, sign-log-abs DS eigenfunction %d',n_bound,nef))
  
    
    
    
end




%% Injectivity radius test

%     trimesh(Tar.SHAPE.surface.TRIV, Tar.SHAPE.surface.X, Tar.SHAPE.surface.Y, Tar.SHAPE.surface.Z, Tar.SHAPE.Gamma(Tar_landmarks(1),:) , ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
%     title(sprintf('Lmk %d, sign-log-abs DS eigenfunction %d',n_bound,nef))
% 









