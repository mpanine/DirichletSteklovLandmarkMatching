% Tests the Steklov matching with landmark harmonic segmentation.

clear all; close all;

shape_settings = shape.Settings('computeMHB', true,...
                        'computeGeoDist', false,...
                        'MHB_settings', basis.MHB_Settings(20),...
                        'computeNormals', false);
cache_settings = cache.Settings('a','b', false);

%% Load source
Src = struct;
Src.SHAPE = shape.compute('./Data/TOSCA_Isometric/vtx_5k/cat0.off', ...
                        shape_settings, cache_settings);

load('./Data/TOSCA_Isometric/Landmarks/cat0')                    
Src_landmarks = landmarks(1:8);                    
                    
Src_basis = ComputeHarmonicBasis(Src.SHAPE, Src_landmarks);

[~, Src_predom] = sort(Src_basis,2,'descend');

Src_dom = Src_predom(:,1);

%% Load target
Tar = struct;
Tar.SHAPE = shape.compute('./Data/TOSCA_Isometric/vtx_5k/dog0.off', ...
                        shape_settings, cache_settings);

load('./Data/TOSCA_Isometric/Landmarks/dog0')                    
Tar_landmarks = landmarks(1:8);                    
                    
Tar_basis = ComputeHarmonicBasis(Tar.SHAPE, Tar_landmarks);

[~, Tar_predom] = sort(Tar_basis,2,'descend');

Tar_dom = Tar_predom(:,1);

%% Show the meshes

subplot(1,2,1)
trimesh(Src.SHAPE.surface.TRIV, Src.SHAPE.surface.X, Src.SHAPE.surface.Y, Src.SHAPE.surface.Z, Src_dom, ...
'FaceColor','flat', ...
'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
hold on
plot3(Src.SHAPE.surface.X(Src_landmarks), Src.SHAPE.surface.Y(Src_landmarks), Src.SHAPE.surface.Z(Src_landmarks),'.','Color','k','MarkerSize',20)

subplot(1,2,2)
trimesh(Tar.SHAPE.surface.TRIV, Tar.SHAPE.surface.X, Tar.SHAPE.surface.Y, Tar.SHAPE.surface.Z, Tar_dom, ...
'FaceColor','flat', ...
'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
hold on
plot3(Tar.SHAPE.surface.X(Tar_landmarks), Tar.SHAPE.surface.Y(Tar_landmarks), Tar.SHAPE.surface.Z(Tar_landmarks),'.','Color','k','MarkerSize',20)


%% Segment split test

landmark_num = 7;

segment = find(Src_dom == landmark_num);

[W, S, segment_TRIV, boundary_old, boundary_new, boundary_edges_old] = splitSegmentSteklov(Src.SHAPE, segment);

figure
trimesh(segment_TRIV, Src.SHAPE.surface.X, Src.SHAPE.surface.Y, Src.SHAPE.surface.Z, ...
'FaceColor','interp', ...
'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
 
hold on
plot3(Src.SHAPE.surface.X(Src_landmarks(landmark_num)), Src.SHAPE.surface.Y(Src_landmarks(landmark_num)), Src.SHAPE.surface.Z(Src_landmarks(landmark_num)),'.','Color','k','MarkerSize',20)
                    
                    