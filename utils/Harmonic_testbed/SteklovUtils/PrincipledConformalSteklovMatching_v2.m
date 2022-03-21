%THIS IS THE MOST RECENT VERSION, LIKELY FINAL.

%This implements the early November 2020 version of the approach, nameley a principled
%approach to getting a as conformal as possible map. It uses LB and Dirichlet-Steklov bases.
%The bases are unified into a large matrix, which simplifies the energy.


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

% Src.SHAPE = shape.compute('./Data/Shrec07_fourleg/383.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Shrec07_fourleg/Landmarks/383') 

% Src.SHAPE = shape.compute('./Data/Faust_remeshed/vtx_5k/tr_reg_000.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Faust_remeshed/Landmarks/tr_reg_000') 

% Src_landmarks = landmarks([1, 6]); 
Src_landmarks = landmarks([1, 6, 8, 10, 12]);
% Src_landmarks = landmarks([1, 6, 8, 10, 12, 14]); 
% Src_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13, 21, 15]); 
% Src_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13, 15]); 
% Src_landmarks = landmarks; 

% Src_landmarks = landmarks([5, 4, 1, 3, 2]);  %For FAUST.

%% Load target
Tar = struct;
% 
% Tar.SHAPE = shape.compute('./Data/Shrec07_fourleg/388.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Shrec07_fourleg/Landmarks/388')   

% Tar.SHAPE = shape.compute('./Data/Shrec07_fourleg/387.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Shrec07_fourleg/Landmarks/387') 

% Tar.SHAPE = shape.compute('./Data/Shrec07_fourleg/390.off', ...
%                         shape_settings, cache_settings);
% load('./Data/Shrec07_fourleg/Landmarks/390')   

Tar.SHAPE = shape.compute('./Data/Faust_remeshed/vtx_5k/tr_reg_022.off', ...
                        shape_settings, cache_settings);
load('./Data/Faust_remeshed/Landmarks/tr_reg_022') 


% Tar_landmarks = landmarks([1, 6]); 
% Tar_landmarks = landmarks([1, 6, 8, 10, 12]);
% Tar_landmarks = landmarks([1, 6, 8, 10, 12, 14]); 
% Tar_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13, 21, 15]);
% Tar_landmarks = landmarks([1, 6, 8, 10, 12, 14, 13, 15]); 
% Tar_landmarks = landmarks;

Tar_landmarks = landmarks([5, 4, 1, 3, 2]);  %For FAUST.

%% Steklov Settings -- Settings for our approach. Includes the 

Steklov_settings = struct;

%%% Specify the way to construct the matrices used in the nearest neighbor search.
    Steklov_settings.NN_type = 'principled'; %The correct way to do it.
%     Steklov_settings.NN_type = 'fast'; %Approximation. Fast, but unprincipled.


%%%% Declare the method for the initial guess. (Boundary correspondence) 
%     Steklov_settings.InitialGuess = 'trivial';  % Uneducated guess.
%     Steklov_settings.InitialGuess = 'landmark_harmonics'; % Match harmonic functions with Kronecker delta BC.
    Steklov_settings.InitialGuess = 'normal_derivatives'; % Uses normal derivatives of the above. Good setting.
%     Steklov_settings.InitialGuess = 'conformal_energy'; % Minimizes conformal energy. Principled, but possibly not as good (?).


Steklov_settings.DS_num_eigs = 10; % Number of Dirichlet-Steklov (DS) eigenfuncation/values PER LANDMARK. Do not set too high (10-20 is plenty).


%------- Settings for refinement at the landmarks ----------------------------------------------------
Steklov_settings.numsubdivisions = 10*ones(size(Src_landmarks)); % 50*ones(size(Src_landmarks)); %Number of subdivisions per triangle at landmark.
Steklov_settings.radii_factor = 0.5; %This is a parameter that whose influence should be studied.

%------- Settings for matching at the landmark circles ----------------------------------------------------

Steklov_settings.num_shifts = 200; %Number of different rotations (360 deg/num_shifts) attempted at a landmark circle. No need to test this, 200 is fine.

%------- Settings for general Dirichlet-Steklov matching ----------------------------------------------------

Steklov_settings.weight_Orthonormality = 1; % Weight for the "Orthonormal" term of the energy. Should be studied. HYPOTHESIS: this has to decrease as the number of landmarks increases.
Steklov_settings.weight_Proper = 1; % Weight for the "Proper funcitonal map" term of the energy. Should be studied. HYPOTHESIS: this and weight_Bijectivity can be the same.
Steklov_settings.weight_Bijectivity = 1; % Weight for the "Bijectivity" term of the energy. Should be studied. HYPOTHESIS: this and weight_Proper can be the same.

%------- Settings for Dirichlet Laplace Beltrami ZoomOut ----------------------------------------------------

Steklov_settings.num_LB_eigs = 400; % Number of LB eigenfunctions to compute.
Steklov_settings.LB_ZO_start = 5; % Start of the LB ZoomOut.
Steklov_settings.LB_ZO_end = 400; % End of the LB ZoomOut, normally the same as num_LB_eigs. - 400 works well.
Steklov_settings.LB_ZO_step = 30; % Step size used in ZoomOut. Can be set to any positive integer, as ZO_end is always used at the end, no matter the step size.


%% Automatic settings -- Not to be tested -- This simplifies some of the code.

Steklov_settings.num_landmarks = length(Src_landmarks);
Steklov_settings.num_DS_functions = length(Src_landmarks) * Steklov_settings.DS_num_eigs;

Steklov_settings.Steklov_num_eigs = Steklov_settings.DS_num_eigs * length(Tar_landmarks); % Total number of DS eigenfunctions used.
Steklov_settings.landmarks_radii = computeLandmarkRadii(Src, Src_landmarks, Tar, Tar_landmarks, Steklov_settings); %Radii of the disks at the landmarks.

%% ============================================== COMPUTATIONS START HERE ==============================================

%% Refine Shapes (landmark circles) and compute all relevant bases


Src_refined = ComputePrincipledSteklovAll_clean_v2(Src, Src_landmarks, Steklov_settings);
Tar_refined = ComputePrincipledSteklovAll_clean_v2(Tar, Tar_landmarks, Steklov_settings);


%% Compute the boundary correspondence -- Mutltiple methods available

switch Steklov_settings.InitialGuess
    
    case 'trivial' %Trivial circle matching.
        [fullp2pTarSrc, fullp2pSrcTar] = PrincipledSteklovBoundaryMatch_trivial(Src_refined, Tar_refined, Steklov_settings); 
        
    case 'landmark_harmonics' % Guess using landmark harmonics.
        [fullp2pTarSrc, fullp2pSrcTar] = PrincipledLandmarkHarmonicInitialGuess(Src_refined,  Tar_refined);
        
    case 'normal_derivatives' %Normal derivatives of landmark harmonics.
        [fullp2pTarSrc, fullp2pSrcTar] = PrincipledSteklovBoundaryMatch_clean_v3(Src_refined,  Tar_refined, Steklov_settings);
        
    case 'conformal_energy'
        [fullp2pTarSrc, fullp2pSrcTar] = PrincipledSteklovBoundaryMatch_clean_v2(Src_refined,  Tar_refined, Steklov_settings);
    
end




%% ZoomOut Procedure -- Produces the desired output

[fullp2pTarSrc_ZO, fullp2pSrcTar_ZO] = Principled_findFullP2P_ZO_clean_v2(Src_refined, Src_landmarks, Tar_refined, Tar_landmarks, fullp2pTarSrc, fullp2pSrcTar, Steklov_settings);




%% ============================ EVERYTHING BELOW IS PLOTS ========================================



%% Plot initial guess coord transfer -- Refined Shapes


ext_init_p2pTarSrc = reinsertBoundaryCorrespondence_Principled(fullp2pTarSrc, Tar_landmarks, Src_landmarks, Src_refined.STEKLOV.segment);
ext_init_p2pSrcTar = reinsertBoundaryCorrespondence_Principled(fullp2pSrcTar, Src_landmarks, Tar_landmarks, Tar_refined.STEKLOV.segment);

figure
subplot(1,2,1)
plotCoordinateTransfer(Tar_refined.SHAPE, Src_refined.SHAPE, ext_init_p2pSrcTar)
hold on
plot3(Tar_refined.SHAPE.surface.X(Tar_landmarks),Tar_refined.SHAPE.surface.Y(Tar_landmarks),Tar_refined.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
title('Initial Guess Tar->Src')

subplot(1,2,2)
plotCoordinateTransfer(Src_refined.SHAPE, Tar_refined.SHAPE, ext_init_p2pTarSrc)
hold on
plot3(Src_refined.SHAPE.surface.X(Src_landmarks),Src_refined.SHAPE.surface.Y(Src_landmarks),Src_refined.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
title('Initial Guess Src->Tar')



%% Plot coord transfer -- ZoomOut, original shape -- This is the solution
% 
figure
subplot(1,2,1)
plotCoordinateTransfer(Tar.SHAPE, Src.SHAPE, fullp2pSrcTar_ZO)
hold on
plot3(Tar.SHAPE.surface.X(Tar_landmarks),Tar.SHAPE.surface.Y(Tar_landmarks),Tar.SHAPE.surface.Z(Tar_landmarks),'.','MarkerSize',20,'color','r')
title('ZoomOut Tar->Src')

subplot(1,2,2)
plotCoordinateTransfer(Src.SHAPE, Tar.SHAPE, fullp2pTarSrc_ZO)
hold on
plot3(Src.SHAPE.surface.X(Src_landmarks),Src.SHAPE.surface.Y(Src_landmarks),Src.SHAPE.surface.Z(Src_landmarks),'.','MarkerSize',20,'color','r')
title('ZoomOut Src->Tar')

%% Plot map lines -- ORIGINAL SHAPES -- This is the solution

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



%% Plot Dirichlet-Steklov eigenfunctions



    nef = 1;
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



    nef = 1;
    
    figure
    subplot(1,2,1)
    trimesh(Src_refined.SHAPE.surface.TRIV, Src_refined.SHAPE.surface.X, Src_refined.SHAPE.surface.Y, Src_refined.SHAPE.surface.Z, Src_refined.STEKLOV.FullBasis(:,Steklov_settings.num_DS_functions + nef)  , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
    title(sprintf('D-LB eigenfunction %d',nef))

    subplot(1,2,2)
    trimesh(Tar_refined.SHAPE.surface.TRIV, Tar_refined.SHAPE.surface.X, Tar_refined.SHAPE.surface.Y, Tar_refined.SHAPE.surface.Z, Tar_refined.STEKLOV.FullBasis(:,Steklov_settings.num_DS_functions + nef) , ...
    'FaceColor','interp', ...
    'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
    title(sprintf('D-LB eigenfunction %d',nef))




