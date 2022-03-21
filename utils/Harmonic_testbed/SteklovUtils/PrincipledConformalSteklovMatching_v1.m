%This implements the early October 2020 version of the approach, nameley a principled
%approach to getting a as conformal as possible map. It uses LB, Steklov
%and Dirichlet-Steklov bases.


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

Src.SHAPE = shape.compute('./Data/Shrec07_fourleg/381.off', ...
                        shape_settings, cache_settings);
load('./Data/Shrec07_fourleg/Landmarks/381') 

    
% Src.SHAPE = shape.compute('./Data/Shrec07_fourleg/394.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Shrec07_fourleg/Landmarks/394') 

% Src.SHAPE = shape.compute('./Data/Shrec07_fourleg/384.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Shrec07_fourleg/Landmarks/384') 

% Src_landmarks = landmarks([1, 6]); 
% Src_landmarks = landmarks([1, 6, 8, 10, 12, 14]); 
% Src_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13, 21, 15]); 
% Src_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13, 15]); 

Src_landmarks = landmarks; 
% 
%% Load target
Tar = struct;

Tar.SHAPE = shape.compute('./Data/Shrec07_fourleg/388.off', ...
                        shape_settings, cache_settings);
load('./Data/Shrec07_fourleg/Landmarks/388')   

% Tar.SHAPE = shape.compute('./Data/Shrec07_fourleg/387.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Shrec07_fourleg/Landmarks/387') 

% Tar.SHAPE = shape.compute('./Data/Shrec07_fourleg/385.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Shrec07_fourleg/Landmarks/385')   

% Tar_landmarks = landmarks([1, 6]); 
% Tar_landmarks = landmarks([1, 6, 8, 10, 12, 14]); 
% Tar_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13, 21, 15]);
% Tar_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13, 15]); 
Tar_landmarks = landmarks;

%% Steklov Settings -- Settings for the entire Steklov approach, including the Dirichlet LB portion of it

OrhogonalityTreatment = 'standard'; %TODO: use this. It will control the way the orthogonality conditions are implemented.

Steklov_settings = struct;

Steklov_settings.DS_num_eigs = 10; % Number of Steklov eigenfuncation/values PER LANDMARK. Do not set too high (20 is plenty).
Steklov_settings.Steklov_num_eigs = Steklov_settings.DS_num_eigs * length(Tar_landmarks); % Number of Steklov eigenfunctions.
                                                                                               % The constant eigenfunction is not included.


%------- Settings for refinement at the landmarks ----------------------------------------------------
Steklov_settings.numsubdivisions = 20*ones(size(Src_landmarks)); % 50*ones(size(Src_landmarks)); %Number of subdivisions per triangle at landmark.
Steklov_settings.radii_factor = 0.5; %This is a parameter that whose influence should be studied.
Steklov_settings.landmarks_radii = computeLandmarkRadii(Src, Src_landmarks, Tar, Tar_landmarks, Steklov_settings); %Radii of the disks at the landmarks.


%------- Settings for matching at the landmark circles ----------------------------------------------------

Steklov_settings.num_shifts = 200; %Number of different rotations attempted at a landmark circle.
Steklov_settings.num_it_ICP = 0;%100  %Number of ICP iterations at a landmark circle.

%------- Settings for general Dirichlet-Steklov matching ----------------------------------------------------

Steklov_settings.weight_Steklov = 1; % Weight for the Steklov spectrum term of the energy. Should be studied.
Steklov_settings.weight_Proper = 1; % Weight for the "Proper funcitonal map" term of the energy. Should be studied.
Steklov_settings.weight_Orthonormality = 1; % Weight for the "Orthonormal" term of the energy. Should be studied.
Steklov_settings.weight_Bijectivity = 1; % Weight for the "Bijectivity" term of the energy. Should be studied.

%------- Settings for Dirichlet Laplace Beltrami ZoomOut ----------------------------------------------------

Steklov_settings.num_LB_eigs = 100;
Steklov_settings.LB_ZO_start = 5;
Steklov_settings.LB_ZO_end = 100;
Steklov_settings.LB_ZO_step = 5;

Steklov_settings.weight_LB_Proper = 1; % Weight for the LB "Proper funcitonal map" term of the energy. Should be studied.
Steklov_settings.weight_LB_Orthonormality = 1; % Weight for the LB "Orthonormal" term of the energy. Should be studied.
Steklov_settings.weight_LB_Bijectivity = 1; % Weight for the LB "Bijectivity" term of the energy. Should be studied.


%% Refine Shapes and compute all relevant bases


Src_refined = ComputePrincipledSteklovAll_clean(Src, Src_landmarks, Steklov_settings);
Tar_refined = ComputePrincipledSteklovAll_clean(Tar, Tar_landmarks, Steklov_settings);



%CURRENTLY THE CONSTANT EIGENFUNCTION IS REMOVED FROM THE SET !!!

% Src_refined = ComputeCentralDirichletSteklovAll_clean(Src, Src_landmarks, Steklov_settings);
% Tar_refined = ComputeCentralDirichletSteklovAll_clean(Tar, Tar_landmarks, Steklov_settings);


%% Computing the boundary correspondence - Exhaustive initial search wrt directions given by harmonic functions

[FmapSrcTar, p2pTarSrc, FmapTarSrc, p2pSrcTar] = PrincipledSteklovBoundaryMatch_clean(Src_refined,  Tar_refined, Steklov_settings);


%% Find the full p2p map -- No LB eigenfunctions

% [fullp2pTarSrc, fullp2pSrcTar] = CentralSteklovFindFullP2P_clean(Src_refined, Src_refined.STEKLOV.evecs(1:Src.SHAPE.nv,:,:), Src_refined.STEKLOV.evals,...
%                                              Tar_refined, Tar_refined.STEKLOV.evecs(1:Tar.SHAPE.nv,:,:), Tar_refined.STEKLOV.evals,...
%                                              FmapSrcTar, FmapTarSrc, Steklov_settings);
%                                          
    
% sprintf('Just before finding full p2p')

[fullp2pTarSrc, fullp2pSrcTar] = PrincipledSteklovFindFullP2P_clean(Src_refined, Src_landmarks, Src_refined.STEKLOV.DS_evecs, Src_refined.STEKLOV.DS_evals,...
                                             Tar_refined, Tar_landmarks, Tar_refined.STEKLOV.DS_evecs, Tar_refined.STEKLOV.DS_evals,...
                                             FmapSrcTar, FmapTarSrc, Steklov_settings);

% sprintf('Just after finding full p2p')



%% Plot images of eigenfunctions Src -> Tar



%     nef = 3;
% 
%     figure
%     subplot(1,2,1)
%     trimesh(Src_refined.SHAPE.surface.TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, Src_refined.STEKLOV.evecs(:,nef)  , ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'k'); axis equal; axis off; colorbar;
%     title(sprintf('Steklov eigenfunction %d on Src',nef))
% 
%     subplot(1,2,2)
%     trimesh(Tar_refined.SHAPE.surface.TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, Tar_refined.STEKLOV.evecs*FmapSrcTar(:,nef) , ...
%     'FaceColor','interp', ...
%     'FaceAlpha', 1, 'EdgeColor', 'k'); axis equal; axis off; colorbar;
%     title('Image on Target via FmapSrcTar')

                                         
                                         
                                                         
%% Find the full p2p map -- Uses ZoomOut on LB eigenfunctions

% sprintf('Just before ZO')

% [fullp2pTarSrc_ZO, fullp2pSrcTar_ZO] = CentralSteklovFindFullP2P_ZO_clean(Src_refined, Src_landmarks, Tar_refined, Tar_landmarks, FmapSrcTar, FmapTarSrc, Steklov_settings);
[fullp2pTarSrc_ZO, fullp2pSrcTar_ZO] = Principled_findFullP2P_ZO_clean(Src_refined, Src_landmarks, Tar_refined, Tar_landmarks, FmapSrcTar, FmapTarSrc, Steklov_settings);
                      
% sprintf('Just after ZO')


%% Plot coord transfer -- REFINED SHAPES

figure
subplot(1,2,1)
plotCoordinateTransfer(Tar_refined.SHAPE, Src_refined.SHAPE, fullp2pSrcTar)
hold on
plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
title('Steklov Coordinate transfer Tar->Src')

subplot(1,2,2)
plotCoordinateTransfer(Src_refined.SHAPE, Tar_refined.SHAPE, fullp2pTarSrc)
hold on
plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
title('Steklov Coordinate transfer Src->Tar')



%% Plot coord transfer -- ZoomOut, original shape

figure
subplot(1,2,1)
plotCoordinateTransfer(Tar.SHAPE, Src.SHAPE, fullp2pSrcTar_ZO)
hold on
plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
title('ZO, original, Steklov Coordinate transfer Tar->Src')

subplot(1,2,2)
plotCoordinateTransfer(Src.SHAPE, Tar.SHAPE, fullp2pTarSrc_ZO)
hold on
plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
title('ZO, original, Steklov Coordinate transfer Src->Tar')

%% Plot map lines -- ORIGINAL SHAPES

figure
subplot(1,2,1)
visualize_map_lines(Src.SHAPE.surface, Tar.SHAPE.surface, fullp2pSrcTar_ZO, Src_landmarks)
% hold on
% plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
title('Steklov ZO Function transfer Tar->Src')

subplot(1,2,2)
visualize_map_lines(Tar.SHAPE.surface, Src.SHAPE.surface, fullp2pTarSrc_ZO, Tar_landmarks)
% hold on
% plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
title('Steklov ZO Function transfer Src->Tar')


%% Plot Steklov eigenfunctions



    nef = 1;

    figure
    subplot(1,2,1)
    trimesh(Src_refined.SHAPE.surface.TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, Src_refined.STEKLOV.Steklov_evecs(:,nef).^2  , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
    title(sprintf('Steklov eigenfunction %d',nef))

    subplot(1,2,2)
    trimesh(Tar_refined.SHAPE.surface.TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, Tar_refined.STEKLOV.Steklov_evecs(:,nef).^2 , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
    title(sprintf('Steklov eigenfunction %d',nef))



%% Plot Dirichlet-Steklov eigenfunctions



    nef = 20;
    n_bound = 1;

    figure
    subplot(1,2,1)
    trimesh(Src_refined.SHAPE.surface.TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, Src_refined.STEKLOV.DS_evecs(:,nef,n_bound).^2  , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
    title(sprintf('Lmk %d, DS eigenfunction %d',n_bound,nef))

    subplot(1,2,2)
    trimesh(Tar_refined.SHAPE.surface.TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, Tar_refined.STEKLOV.DS_evecs(:,nef,n_bound).^2 , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
    title(sprintf('Lmk %d, DS eigenfunction %d',n_bound,nef))

%% Plot Dirichlet Laplace eigenfunctions



    nef = 100;
   
    figure
    subplot(1,2,1)
    trimesh(Src_refined.SHAPE.surface.TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, Src_refined.STEKLOV.LB_evecs(:,nef)  , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
    title(sprintf('D-LB eigenfunction %d',nef))

    subplot(1,2,2)
    trimesh(Tar_refined.SHAPE.surface.TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, Tar_refined.STEKLOV.LB_evecs(:,nef) , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
    title(sprintf('D-LB eigenfunction %d',nef))




