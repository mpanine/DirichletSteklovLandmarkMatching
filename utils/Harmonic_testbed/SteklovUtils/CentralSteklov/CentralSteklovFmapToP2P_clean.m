function [p2pTarSrc, p2pSrcTar] = CentralSteklovFmapToP2P_clean(Src_evecs_bound, Src_evals, Tar_evecs_bound, Tar_evals, FmapSrcTar, FmapTarSrc, Steklov_settings, ZOsize)
%Implements bijectivity condition as well as the others

%Fmap: functional map Src -> Tar 
%p2p: point-to-point map Tar -> Src

% ZOsize = Steklov_settings.num_eigs;
SteklovWeight = sqrt(Steklov_settings.weight_Steklov);
ProperWeight = sqrt(Steklov_settings.weight_Proper);
OrthoWeight = sqrt(Steklov_settings.weight_Orthonormality);
BijectWeight = sqrt(Steklov_settings.weight_Bijectivity);

for i = 1:length(Src_evecs_bound)
    
%     size(Fmap{i})
    
%% Computing the p2pTarSrc map.

    SrcFTarSrc = SteklovWeight * Src_evecs_bound{i}(:, 1:ZOsize)*Src_evals(1:ZOsize,1:ZOsize) +... %Steklov term
           ProperWeight * Src_evecs_bound{i}(:,1:ZOsize)+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Src_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar'+... %Orthonormality term.
           BijectWeight * Src_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc; %Bijectivity term
       
    TarFSrcTar = SteklovWeight * Tar_evecs_bound{i}(:, 1:ZOsize)*Tar_evals(1:ZOsize,1:ZOsize)*FmapSrcTar +... %Steklov term
           ProperWeight * Tar_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Tar_evecs_bound{i}(:,1:ZOsize)+...%Orthonormality term.
           BijectWeight * Tar_evecs_bound{i}(:,1:ZOsize); %Bijectivity term
    

    p2pTarSrc{i} = annquery(SrcFTarSrc', TarFSrcTar', 1);
    
    
    
    %% Computing the p2pSrcTar map.

    SrcFSrcTar = SteklovWeight * Src_evecs_bound{i}(:, 1:ZOsize)*Src_evals(1:ZOsize,1:ZOsize)*FmapTarSrc +... %Steklov term
           ProperWeight * Src_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Src_evecs_bound{i}(:,1:ZOsize)+... %Orthonormality term.
           BijectWeight * Src_evecs_bound{i}(:,1:ZOsize); %Bijectivity term
       
    TarFTarSrc = SteklovWeight * Tar_evecs_bound{i}(:, 1:ZOsize)*Tar_evals(1:ZOsize,1:ZOsize) +... %Steklov term
           ProperWeight * Tar_evecs_bound{i}(:,1:ZOsize)+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Tar_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc'+...%Orthonormality term.
           BijectWeight * Tar_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar; %Bijectivity term
       
 
    p2pSrcTar{i} = annquery(TarFTarSrc', SrcFSrcTar', 1);
 
    
end




end