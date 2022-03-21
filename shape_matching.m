restoredefaultpath;
addpath(genpath('./utils'));


mesh_dir = './data/'; 
lmk_dir = './data/';

shape_name_Src = 'source'
shape_name_Tar = 'target'

shape_settings = shape.Settings('computeMHB', true,...
                        'computeGeoDist', false,...
                        'MHB_settings', basis.MHB_Settings(20),...
                        'computeNormals', false);
%%
cache_settings = cache.Settings();
%% Load source
Src = struct;
Src.SHAPE = shape.compute([mesh_dir shape_name_Src '.off'], ...
               shape_settings, cache_settings);


%% Load target
Tar = struct;
Tar.SHAPE = shape.compute([mesh_dir shape_name_Tar '.off'], ...
               shape_settings, cache_settings);
%% Load landmarks
landmarks_obj = load([lmk_dir,shape_name_Src,'_landmarks.mat']);  
Src_landmarks = landmarks_obj.landmarks;

landmarks_obj = load([lmk_dir,shape_name_Tar,'_landmarks.mat']);  
Tar_landmarks = landmarks_obj.landmarks;
%% Visualize landmarks
figure;visualize.plot_vtx_on_mesh(Src.SHAPE,Src_landmarks);
figure;visualize.plot_vtx_on_mesh(Tar.SHAPE,Tar_landmarks);

%% Define the parameters of our methods
NN_type = 'fast';
InitialGuess = 'normal_derivatives';

DS_num_eigs = 10;

radii_factor = 0.5;

weight_Orthonormality = 1;
weight_Proper = 1;
weight_Bijectivity = 1;

num_LB_eigs = 120;
ZO_start = 20;
ZO_step = 5;
%% Compute LB eigenbasis and link landmarks to the source and target
Src.SHAPE = MESH.compute_LaplacianBasis(Src.SHAPE, num_LB_eigs);
Tar.SHAPE = MESH.compute_LaplacianBasis(Tar.SHAPE, num_LB_eigs);

Src.landmarks = Src_landmarks;
Tar.landmarks = Tar_landmarks;
%% Create the Steklov settings struct
num_landmarks = length(Src_landmarks);
Steklov_settings = compute_steklov_settings(...
num_landmarks, NN_type,InitialGuess,DS_num_eigs,radii_factor,...
weight_Orthonormality,weight_Proper,weight_Bijectivity,...
num_LB_eigs,ZO_start,ZO_step);
%% Compute the landmark radii
Steklov_settings.landmarks_radii = computeLandmarkRadii(Src, Src_landmarks, Tar, Tar_landmarks, Steklov_settings); %Radii of the disks at the landmarks.
%% Compute the map
% profile on
[Src_refined,Tar_refined,fullp2pTarSrc, ...
fullp2pSrcTar,fullp2pTarSrc_ZO, fullp2pSrcTar_ZO] = ...
compute_steklov(Src, Src_landmarks, ...
Tar, Tar_landmarks, Steklov_settings);
% profile viewer
%% Visualize the transfer
visualize.visualize_map_colors(Src.SHAPE,Tar.SHAPE,fullp2pSrcTar_ZO,'IfShowCoverage',false,'OverlayAxis','x');
title('Shape Matching using Dirichlet-Steklov Eigenfunctions');
view(0,90);