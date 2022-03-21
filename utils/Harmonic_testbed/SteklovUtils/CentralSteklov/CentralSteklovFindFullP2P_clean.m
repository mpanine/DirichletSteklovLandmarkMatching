function [fullp2pTarSrc, fullp2pSrcTar] = CentralSteklovFindFullP2P_clean(Src_refined, Src_landmarks, Src_eF, Src_evals, Tar_refined, Tar_landmarks, Tar_eF, Tar_evals, FmapSrcTar, FmapTarSrc, Steklov_settings)
%CentralSteklovFindFullP2P_clean: Converts the functional maps obtained at the boundary into point to point
%maps in the bulk.

% num_bounds = size(Src_eF,3);
num_eigs = Steklov_settings.num_eigs; %The number of considerd eigenvalues.

SteklovWeight = sqrt(Steklov_settings.weight_Steklov);
ProperWeight = sqrt(Steklov_settings.weight_Proper);
OrthoWeight = sqrt(Steklov_settings.weight_Orthonormality);
BijectWeight = sqrt(Steklov_settings.weight_Bijectivity);


  

%% Computing the fullp2pTarSrc map.

    SrcFTarSrc = SteklovWeight * Src_eF(:,1:num_eigs)*Src_evals(1:num_eigs, 1:num_eigs) +... %Steklov term
           ProperWeight * Src_eF(:,1:num_eigs)+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Src_eF(:,1:num_eigs)*FmapSrcTar(1:num_eigs, 1:num_eigs)'+... %Orthonormality term.
           BijectWeight * Src_eF(:,1:num_eigs)*FmapTarSrc(1:num_eigs, 1:num_eigs); %Bijectivity term
       
    TarFTarSrc = SteklovWeight * Tar_eF(:,1:num_eigs)*Tar_evals(1:num_eigs, 1:num_eigs)*FmapSrcTar(1:num_eigs, 1:num_eigs) +... %Steklov term
           ProperWeight * Tar_eF(:,1:num_eigs)*FmapSrcTar(1:num_eigs, 1:num_eigs)+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Tar_eF(:,1:num_eigs)+...%Orthonormality term.
           Tar_eF(:,1:num_eigs); %Bijectivity term
    
    
%     fullp2pTarSrc = annquery(SrcFTarSrc', TarFTarSrc', 1);
%     fullp2pTarSrc = SeparateInteriorExteriorNearestNeighbor(Src_refined, SrcFTarSrc, Tar_refined, TarFTarSrc);
    fullp2pTarSrc = SeparateInteriorBoundaryLandmarkNearestNeighbor(Src_refined, SrcFTarSrc, Src_landmarks, Tar_refined, TarFTarSrc, Tar_landmarks);
    
    
    %% Computing the fullp2pSrcTar map.

    SrcFSrcTar = SteklovWeight * Src_eF(:,1:num_eigs)*Src_evals(1:num_eigs, 1:num_eigs)*FmapTarSrc(1:num_eigs, 1:num_eigs) +... %Steklov term
           ProperWeight * Src_eF(:,1:num_eigs)*FmapTarSrc(1:num_eigs, 1:num_eigs)+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Src_eF(:,1:num_eigs)+... %Orthonormality term.
           BijectWeight * Src_eF(:,1:num_eigs); %Bijectivity term
       
    TarFSrcTar = SteklovWeight * Tar_eF(:,1:num_eigs)*Tar_evals(1:num_eigs, 1:num_eigs) +... %Steklov term
           ProperWeight * Tar_eF(:,1:num_eigs)+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Tar_eF(:,1:num_eigs)*FmapTarSrc(1:num_eigs, 1:num_eigs)'+...%Orthonormality term.
           BijectWeight * Tar_eF(:,1:num_eigs)*FmapSrcTar(1:num_eigs, 1:num_eigs); %Bijectivity term
    
    
%     fullp2pSrcTar = annquery(TarFSrcTar', SrcFSrcTar', 1);
%     fullp2pSrcTar = SeparateInteriorExteriorNearestNeighbor(Tar_refined, TarFSrcTar, Src_refined, SrcFSrcTar);
    fullp2pSrcTar = SeparateInteriorBoundaryLandmarkNearestNeighbor(Tar_refined, TarFSrcTar, Tar_landmarks, Src_refined, SrcFSrcTar, Src_landmarks);



end