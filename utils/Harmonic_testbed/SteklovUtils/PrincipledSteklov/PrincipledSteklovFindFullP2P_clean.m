function [fullp2pTarSrc, fullp2pSrcTar] = PrincipledSteklovFindFullP2P_clean(Src_refined, Src_landmarks, Src_eF, Src_evals, Tar_refined, Tar_landmarks, Tar_eF, Tar_evals, FmapSrcTar, FmapTarSrc, Steklov_settings)
%PrincipledSteklovFindFullP2P_clean: Converts the functional maps obtained at the boundary into point to point
%maps in the bulk.

%%TODO: finish implementing the principled energy

num_bounds = size(Src_eF,3);
num_eigs = Steklov_settings.DS_num_eigs; %The number of considerd eigenvalues.

SteklovWeight = sqrt(Steklov_settings.weight_Steklov);
ProperWeight = sqrt(Steklov_settings.weight_Proper);
OrthoWeight = sqrt(Steklov_settings.weight_Orthonormality);
BijectWeight = sqrt(Steklov_settings.weight_Bijectivity);


%% Normalize the DS eigenvectors

Src_evecs_combined = zeros(size(Src_eF,1), size(Src_eF,2)*size(Src_eF,3));
Tar_evecs_combined = zeros(size(Tar_eF,1), size(Tar_eF,2)*size(Tar_eF,3));

Src_evals_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);
Tar_evals_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);

FmapSrcTar_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);
FmapTarSrc_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);


for i = 1:num_bounds
        
       start_ind = (i-1)*num_eigs + 1;
       end_ind = i*num_eigs;
       
%        Src_evecs_combined(:, start_ind:end_ind ) = Src_eF(:,:,i);
%        Tar_evecs_combined(:, start_ind:end_ind ) = Tar_eF(:,:,i);
       
       Src_evecs_combined_normalized(:, start_ind:end_ind ) = Src_eF(:,:,i) * diag( diag( Src_refined.STEKLOV.DS_evals(:,:,i) ).^-0.5 );
       Tar_evecs_combined_normalized(:, start_ind:end_ind ) = Tar_eF(:,:,i) * diag( diag( Tar_refined.STEKLOV.DS_evals(:,:,i) ).^-0.5 );

       
       Src_evals_combined(start_ind:end_ind, start_ind:end_ind) = Src_evals(:,:,i);
       Tar_evals_combined(start_ind:end_ind, start_ind:end_ind) = Tar_evals(:,:,i);
    
       FmapSrcTar_combined(start_ind:end_ind, start_ind:end_ind) = FmapSrcTar{i};
       FmapTarSrc_combined(start_ind:end_ind, start_ind:end_ind) = FmapTarSrc{i};
end



% %     Src_sqrt_evals = diag( diag( Src_refined.STEKLOV.evals ).^0.5 );
%     Src_inv_sqrt_evals = diag( diag( Src_refined.STEKLOV.DS_evals ).^-0.5 );
% %     Src_inv_sqrt_evals(1,:) = []; %Remove the constant 
% %     Src_inv_sqrt_evals(:,1) = [];
%     
% %     Tar_sqrt_evals = diag( diag( Tar_refined.STEKLOV.evals ).^0.5 );
%     Tar_inv_sqrt_evals = diag( diag( Tar_refined.STEKLOV.DS_evals ).^-0.5 );
% 
%     Src_evecs_normalized = Src_refined.STEKLOV.DS_evecs * Src_inv_sqrt_evals; %DO NOT FORGET TO USE THESE IN THE CODE BELOW!!!
%     Tar_evecs_normalized = Tar_refined.STEKLOV.DS_evecs * Tar_inv_sqrt_evals;

  

%% Computing the fullp2pTarSrc map.

%     SrcFTarSrc = SteklovWeight * Src_eF(:,1:num_eigs)*Src_evals(1:num_eigs, 1:num_eigs) +... %Steklov term
%            ProperWeight * Src_eF(:,1:num_eigs)+... %"Proper" functional map term (pullback of p2p map)
%            OrthoWeight * Src_eF(:,1:num_eigs)*FmapSrcTar(1:num_eigs, 1:num_eigs)'+... %Orthonormality term.
%            BijectWeight * Src_eF(:,1:num_eigs)*FmapTarSrc(1:num_eigs, 1:num_eigs); %Bijectivity term
%        
%     TarFTarSrc = SteklovWeight * Tar_eF(:,1:num_eigs)*Tar_evals(1:num_eigs, 1:num_eigs)*FmapSrcTar(1:num_eigs, 1:num_eigs) +... %Steklov term
%            ProperWeight * Tar_eF(:,1:num_eigs)*FmapSrcTar(1:num_eigs, 1:num_eigs)+... %"Proper" functional map term (pullback of p2p map)
%            OrthoWeight * Tar_eF(:,1:num_eigs)+...%Orthonormality term.
%            Tar_eF(:,1:num_eigs); %Bijectivity term
%        
    SrcFTarSrc = ProperWeight * Src_evecs_combined_normalized+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Src_evecs_combined_normalized * FmapSrcTar_combined'+... %Orthonormality term.
           BijectWeight * Src_evecs_combined_normalized * FmapTarSrc_combined; %Bijectivity term
       
    TarFTarSrc = ProperWeight * Tar_evecs_combined_normalized * FmapSrcTar_combined+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Tar_evecs_combined_normalized +...%Orthonormality term.
           Tar_evecs_combined_normalized; %Bijectivity term
    
    
    
%     fullp2pTarSrc = annquery(SrcFTarSrc', TarFTarSrc', 1);
%     fullp2pTarSrc = SeparateInteriorExteriorNearestNeighbor(Src_refined, SrcFTarSrc, Tar_refined, TarFTarSrc);
    fullp2pTarSrc = SeparateInteriorBoundaryLandmarkNearestNeighbor(Src_refined, SrcFTarSrc, Src_landmarks, Tar_refined, TarFTarSrc, Tar_landmarks);
    
    
    %% Computing the fullp2pSrcTar map.

%     SrcFSrcTar = SteklovWeight * Src_eF(:,1:num_eigs)*Src_evals(1:num_eigs, 1:num_eigs)*FmapTarSrc(1:num_eigs, 1:num_eigs) +... %Steklov term
%            ProperWeight * Src_eF(:,1:num_eigs)*FmapTarSrc(1:num_eigs, 1:num_eigs)+... %"Proper" functional map term (pullback of p2p map)
%            OrthoWeight * Src_eF(:,1:num_eigs)+... %Orthonormality term.
%            BijectWeight * Src_eF(:,1:num_eigs); %Bijectivity term
%        
%     TarFSrcTar = SteklovWeight * Tar_eF(:,1:num_eigs)*Tar_evals(1:num_eigs, 1:num_eigs) +... %Steklov term
%            ProperWeight * Tar_eF(:,1:num_eigs)+... %"Proper" functional map term (pullback of p2p map)
%            OrthoWeight * Tar_eF(:,1:num_eigs)*FmapTarSrc(1:num_eigs, 1:num_eigs)'+...%Orthonormality term.
%            BijectWeight * Tar_eF(:,1:num_eigs)*FmapSrcTar(1:num_eigs, 1:num_eigs); %Bijectivity term
       
    SrcFSrcTar = ProperWeight * Src_evecs_combined_normalized * FmapTarSrc_combined+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Src_evecs_combined_normalized+... %Orthonormality term.
           BijectWeight * Src_evecs_combined_normalized; %Bijectivity term
       
    TarFSrcTar = ProperWeight * Tar_evecs_combined_normalized +... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Tar_evecs_combined_normalized  * FmapTarSrc_combined'+...%Orthonormality term.
           BijectWeight * Tar_evecs_combined_normalized * FmapSrcTar_combined; %Bijectivity term
    
    
%     fullp2pSrcTar = annquery(TarFSrcTar', SrcFSrcTar', 1);
%     fullp2pSrcTar = SeparateInteriorExteriorNearestNeighbor(Tar_refined, TarFSrcTar, Src_refined, SrcFSrcTar);
    fullp2pSrcTar = SeparateInteriorBoundaryLandmarkNearestNeighbor(Tar_refined, TarFSrcTar, Tar_landmarks, Src_refined, SrcFSrcTar, Src_landmarks);



end