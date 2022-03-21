function [Steklov_settings] = compute_steklov_settings(...
    num_landmarks, NN_type,InitialGuess,DS_num_eigs,radii_factor,...
    weight_Orthonormality,weight_Proper,weight_Bijectivity,...
    num_LB_eigs,ZO_start,ZO_step)
%COMPUTE_STEKLOV_SETTINGS Summary of this function goes here
%   Detailed explanation goes here
    %% Steklov Settings -- Settings for our approach.
    
    Steklov_settings = struct;
    
    Steklov_settings.num_LB_margin = 10; % number of additional LB eigenpairs computed (if some of the eigenvalues happen to be below eig_threshold).
    Steklov_settings.eig_threshold = 1e-4; % threshold below which the eigenvalues are considered as null.

    %%% Specify the way to construct the matrices used in the nearest neighbor search.
    Steklov_settings.NN_type = NN_type;
%         Steklov_settings.NN_type = 'principled'; %The correct way to do it.
%         Steklov_settings.NN_type = 'fast'; %Approximation. Fast, but unprincipled.


    %%%% Declare the method for the initial guess. (Boundary correspondence) 
    Steklov_settings.InitialGuess = InitialGuess;
%         Steklov_settings.InitialGuess = 'trivial';  % Uneducated guess.
%         Steklov_settings.InitialGuess = 'landmark_harmonics'; % Match harmonic functions with Kronecker delta BC.
%         Steklov_settings.InitialGuess = 'normal_derivatives'; % Uses normal derivatives of the above. Good setting.
%         Steklov_settings.InitialGuess = 'conformal_energy'; % Minimizes conformal energy. Principled, but possibly not as good (?).

    
%     Steklov_settings.DS_num_eigs = 10; % Number of Dirichlet-Steklov (DS) eigenfuncation/values PER LANDMARK. Do not set too high (10-20 is plenty).
    Steklov_settings.DS_num_eigs = DS_num_eigs;

    %------- Settings for refinement at the landmarks ----------------------------------------------------
    Steklov_settings.numsubdivisions = 10*ones(num_landmarks,1); % 50*ones(size(Src_landmarks)); %Number of subdivisions per triangle at landmark.
%     Steklov_settings.radii_factor = 0.5; %This is a parameter that whose influence should be studied.
    Steklov_settings.radii_factor = radii_factor;

    %------- Settings for matching at the landmark circles ----------------------------------------------------

    Steklov_settings.num_shifts = 200; %Number of different rotations (360 deg/num_shifts) attempted at a landmark circle. No need to test this, 200 is fine.

    %------- Settings for general Dirichlet-Steklov matching ----------------------------------------------------

%     Steklov_settings.weight_Orthonormality = 1; % Weight for the "Orthonormal" term of the energy. Should be studied. HYPOTHESIS: this has to decrease as the number of landmarks increases.
%     Steklov_settings.weight_Proper = 1; % Weight for the "Proper funcitonal map" term of the energy. Should be studied. HYPOTHESIS: this and weight_Bijectivity can be the same.
%     Steklov_settings.weight_Bijectivity = 1; % Weight for the "Bijectivity" term of the energy. Should be studied. HYPOTHESIS: this and weight_Proper can be the same.

    Steklov_settings.weight_Orthonormality = weight_Orthonormality;
    Steklov_settings.weight_Proper = weight_Proper;
    Steklov_settings.weight_Bijectivity = weight_Bijectivity;
    %------- Settings for Dirichlet Laplace Beltrami ZoomOut ----------------------------------------------------

%     Steklov_settings.num_LB_eigs = 400; % Number of LB eigenfunctions to compute.
%     Steklov_settings.LB_ZO_start = 5; % Start of the LB ZoomOut.
%     Steklov_settings.LB_ZO_end = 400; % End of the LB ZoomOut, normally the same as num_LB_eigs. - 400 works well.
%     Steklov_settings.LB_ZO_step = 30; % Step size used in ZoomOut. Can be set to any positive integer, as ZO_end is always used at the end, no matter the step size.

    Steklov_settings.num_LB_eigs = num_LB_eigs;
    Steklov_settings.LB_ZO_start = ZO_start;
    Steklov_settings.LB_ZO_end = num_LB_eigs;
    Steklov_settings.LB_ZO_step = ZO_step;

    %% Automatic settings -- Not to be tested -- This simplifies some of the code.

    Steklov_settings.num_landmarks = num_landmarks;
    Steklov_settings.num_DS_functions = num_landmarks * Steklov_settings.DS_num_eigs;

    Steklov_settings.Steklov_num_eigs = Steklov_settings.DS_num_eigs * num_landmarks; % Total number of DS eigenfunctions used.
end

