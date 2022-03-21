% Tests the Landmark Steklov Approach


%% Technical stuff

clear all;
close all;

shape_settings = shape.Settings('computeMHB', true,...
                        'computeGeoDist', false,...
                        'MHB_settings', basis.MHB_Settings(20),...
                        'computeNormals', false);
cache_settings = cache.Settings('a','b', false);




%% Load source
Src = struct;

% Src.SHAPE = shape.compute('./Data/bunny_conformal_short/sphere_10242.off', ...
%                         shape_settings, cache_settings);
% load('./Data/bunny_conformal_short/Landmarks/sphere_10242')

Src.SHAPE = shape.compute('./Data/Shrec07_fourleg/381.off', ...
                        shape_settings, cache_settings);
load('./Data/Shrec07_fourleg/Landmarks/381') 


% Src_landmarks = landmarks(1:6);
% Src_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13]);      
Src_landmarks = landmarks([1, 6, 8, 10, 12, 14]); 
                    
% Src_basis = ComputeHarmonicBasis(Src.SHAPE, Src_landmarks);

% [~, Src_predom] = sort(Src_basis,2,'descend');

% Src_dom = Src_predom(:,1);

%% Load target
Tar = struct;

% Tar.SHAPE = shape.compute('./Data/bunny_conformal_short/flow00000217.off', ...
%                         shape_settings, cache_settings);
% load('./Data/bunny_conformal_short/Landmarks/flow00000217') 

% Tar.SHAPE = shape.compute('./Data/bunny_conformal_short/sphere_10242.off', ...
%                         shape_settings, cache_settings);
% load('./Data/bunny_conformal_short/Landmarks/sphere_10242')


Tar.SHAPE = shape.compute('./Data/Shrec07_fourleg/388.off', ...
                        shape_settings, cache_settings);
load('./Data/Shrec07_fourleg/Landmarks/388')   

% Tar_landmarks = landmarks(1:6);
% Tar_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13]);  
Tar_landmarks = landmarks([1, 6, 8, 10, 12, 14]); 
                    
% Tar_basis = ComputeHarmonicBasis(Tar.SHAPE, Tar_landmarks);

% [~, Tar_predom] = sort(Tar_basis,2,'descend');

% Tar_dom = Tar_predom(:,1);

%% Plot shapes and landmarks

figure
subplot(1,2,1)
trimesh(Src.SHAPE.surface.TRIV, Src.SHAPE.surface.X, Src.SHAPE.surface.Y, Src.SHAPE.surface.Z, ...
'FaceColor','interp', ...
'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
hold on
plot3(Src.SHAPE.surface.X(Src_landmarks,:), Src.SHAPE.surface.Y(Src_landmarks,:), Src.SHAPE.surface.Z(Src_landmarks,:),'.','Color','r','MarkerSize',10)

title('Source')

subplot(1,2,2)
trimesh(Tar.SHAPE.surface.TRIV, Tar.SHAPE.surface.X, Tar.SHAPE.surface.Y, Tar.SHAPE.surface.Z, ...
'FaceColor','interp', ...
'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
hold on
plot3(Tar.SHAPE.surface.X(Tar_landmarks,:), Tar.SHAPE.surface.Y(Tar_landmarks,:), Tar.SHAPE.surface.Z(Tar_landmarks,:),'.','Color','r','MarkerSize',10)
title('Target')



%% Internal Settings

numsubdivisions = 5*ones(size(Src_landmarks)); %Number of subdivisions per triangle at landmark.
radii_factor = 0.5;

num_eigs = 100;


%% Landmark radii 

Src_shortest_edges = findShortestEdgesAtLandmarks(Src, Src_landmarks);
Tar_shortest_edges = findShortestEdgesAtLandmarks(Tar, Tar_landmarks);
    
        landmarks_radii = zeros(size(Src_shortest_edges));
        for land_ind = 1:length(landmarks_radii)
            
            landmarks_radii(land_ind) = min(Src_shortest_edges(land_ind), Tar_shortest_edges(land_ind) );
            
        end
    
        landmarks_radii = radii_factor * landmarks_radii;


%% Test without projection on original vertices

[Src_refined, Src_segment_TRIV, Src_evecs, Src_evals, Src_boundary_edges_refined, Src_boundary_old] = ComputeCentralSteklovBasis(Src, Src_landmarks, landmarks_radii, numsubdivisions, num_eigs);

[Tar_refined, Tar_segment_TRIV, Tar_evecs, Tar_evals, Tar_boundary_edges_refined, Tar_boundary_old] = ComputeCentralSteklovBasis(Tar, Tar_landmarks, landmarks_radii, numsubdivisions, num_eigs);


%% Plot the refined shape and its boundary



% plotShapeAndLandmarks(Src_refined.SHAPE, Src_segment, Src_boundary_old)
% title('Shape, segment and boundary')
% hold on
% 

%% Extend to full refined mesh by padding with NaNs

Src_nv = size(Src_refined.SHAPE.surface.VERT,1);
Src_segment = 1:Src_nv;
Src_segment(Src_landmarks) = [];

Src_eF = extendSegmentFunctions(Src_evecs, Src_segment, Src_nv);


Tar_nv = size(Tar_refined.SHAPE.surface.VERT,1);
Tar_segment = 1:Tar_nv;
Tar_segment(Tar_landmarks) = [];

Tar_eF = extendSegmentFunctions(Tar_evecs, Tar_segment, Tar_nv);


%% Plot eigenfunctions

nef = 8;
figure
subplot(1,2,1)
trimesh(Src_segment_TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, abs( Src_eF(:,nef) ), ...
'FaceColor','interp', ...
'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
title(sprintf('Abs Steklov Eigenfunction %d',nef))

subplot(1,2,2)
trimesh(Tar_segment_TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, abs( Tar_eF(:,nef) ), ...
'FaceColor','interp', ...
'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
title(sprintf('Abs Steklov Eigenfunction %d',nef))



%% Fmap Optimization

if false

    rf = @(F) reshape(F, [num_eigs, num_eigs]);

    max_eig = max(max(abs(diag(Src_evals))), max(abs(diag(Tar_evals))));

    Smask = (repmat(diag(Src_evals)/max_eig, [1,num_eigs]) - repmat(diag(Tar_evals)'/max_eig, [num_eigs,1])).^2; %Mask for the Steklov problem.

    Smask = Smask/norm(Smask);

    figure
    imagesc(Smask)
    title('Mask')


    E_SteklovCommutativity = @(F) sum(sum(( rf(F) .^2 .* Smask)/2));
    dE_SteklovCommutativity  = @(F) reshape(    rf(F).*Smask   ,[],1);

    E_orthonormal = @(F) sum(sum( (  ...
                                    rf(F)'*rf(F)...
                                    - eye(num_eigs)  ...
                              ).^2  ))/2;

    dE_orthonormal = @(F) 2 * reshape( ...
                                          rf(F) * ( rf(F)'*rf(F) - eye(num_eigs) ) ...
                                , [],1);

    % 
    % E = @(F) deal(  E_orthonormal(F) + E_SteklovCommutativity(F), ...
    %                 dE_orthonormal(F) + dE_SteklovCommutativity(F) );

    E = @(F) deal( E_SteklovCommutativity(F), ...
                   dE_SteklovCommutativity(F) );

    proj = @(F) F;

    Ini = 0.01*rand(num_eigs,num_eigs);
    % Ini(1,1) = 1;

    Ini = Ini(:);

    options.verbose = 1;
    options.maxIter = 1000;
    options.optTol = 1e-9; % 1e-5 default.
    options.progTol = 1e-10; % 1e-9 default.


    Fmap = minConf_PQN(E,...
                    Ini,...
                    proj,options);


    Fmap = reshape(Fmap,[num_eigs,num_eigs]) ;


    figure
    subplot(1,2,1)
    imagesc(Fmap)
    title('Fmap')
    subplot(1,2,2)
    imagesc(Fmap'*Fmap)
    title('F^T*F')

end

%% Get p2p map

% p2p = annquery(abs(Src_eF(:,2:9))', abs(Tar_eF(:,2:9))', 1);
% plotCoordinateTransfer(Tar_refined.SHAPE, Src_refined.SHAPE, p2p)
% 
% 
% figure
% visualize_map_lines(Tar_refined.SHAPE.surface, Src_refined.SHAPE.surface, p2p, 1)
% 


%% Prepare boundary for plotting

Src_boundary_reorder = boundaryReorder(Src_boundary_edges_refined);
Src_reorder_lengths = computeLengthAlongBoundary(Src_boundary_reorder, Src_refined.SHAPE);

Tar_boundary_reorder = boundaryReorder(Tar_boundary_edges_refined);
Tar_reorder_lengths = computeLengthAlongBoundary(Tar_boundary_reorder, Tar_refined.SHAPE);

%% Plot at boundary

funs_to_plot = 2:5;

plotFunsAtBoundaries(abs(Src_eF), Src_boundary_reorder, Src_reorder_lengths, funs_to_plot)

plotFunsAtBoundaries(abs(Tar_eF), Tar_boundary_reorder, Tar_reorder_lengths, funs_to_plot)


%% Plot Results

if false


    toplot = [2 3 4]+0; %Functions to plot

    Mapped_Src = Src_eF*Fmap;

    figure

    %     subplot(1,2,1)
        trimesh(Src.SHAPE.surface.TRIV, Mapped_Src(:,toplot(1)), Mapped_Src(:,toplot(2)), Mapped_Src(:,toplot(3)), ...
            Mapped_Src(:,1),...
            'FaceColor','flat', ...
            'FaceAlpha', 1, 'EdgeColor', 'b'); axis equal;    
    %     title('Embedding Src')

        hold on;

    %     subplot(1,2,2)
        trimesh(Tar.SHAPE.surface.TRIV, Tar_eF(:,toplot(1)), Tar_eF(:,toplot(2)), Tar_eF(:,toplot(3)), ...
            Tar_eF(:,1),...
            'FaceColor','flat', ...
            'FaceAlpha', 1, 'EdgeColor', 'r'); axis equal;    
        title(['Embedding funs ' num2str(toplot)])




end



%% Embedding plot

toplot = [2 3 4]+6; %Functions to plot

figure
    
%     subplot(1,2,1)
    trimesh(Src_refined.SHAPE.surface.TRIV,  Src_eF(:,toplot(1)) ,   Src_eF(:,toplot(2)) ,  Src_eF(:,toplot(3)) , ...
        Src_eF(:,1),...
        'FaceColor','flat', ...
        'FaceAlpha', 1, 'EdgeColor', 'b'); axis equal;    
%     title('Embedding Src')

    hold on;
        
%     subplot(1,2,2)
    trimesh(abs(  Tar_refined.SHAPE.surface.TRIV ),   Tar_eF(:,toplot(1)) ,  Tar_eF(:,toplot(2)) , Tar_eF(:,toplot(3)) , ...
        Tar_eF(:,1),...
        'FaceColor','flat', ...
        'FaceAlpha', 1, 'EdgeColor', 'r'); axis equal;    
    title(['Embedding funs ' num2str(toplot)])









