function [p2pTarSrc, p2pSrcTar] = DirichletSteklovFmapToP2P_clean(Src_evecs_bound, Src_evals, Tar_evecs_bound, Tar_evals, FmapSrcTar, FmapTarSrc, Steklov_settings)
%Implements bijectivity condition as well as the others

%Fmap: functional map Src -> Tar 
%p2p: point-to-point map Tar -> Src

ZOsize = Steklov_settings.num_eigs;
SteklovWeight = sqrt(Steklov_settings.weight_Steklov);
ProperWeight = sqrt(Steklov_settings.weight_Proper);
OrthoWeight = sqrt(Steklov_settings.weight_Orthonormality);
BijectWeight = sqrt(Steklov_settings.weight_Bijectivity);

for i = 1:length(Src_evecs_bound)
    
%     size(Fmap{i})
    
%% Computing the p2pTarSrc map.

    SrcFTarSrc = SteklovWeight * Src_evecs_bound{i}(:, 1:ZOsize)*Src_evals(1:ZOsize,1:ZOsize,i) +... %Steklov term
           ProperWeight * Src_evecs_bound{i}(:,1:ZOsize)+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Src_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar{i}'+... %Orthonormality term.
           BijectWeight * Src_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc{i}; %Bijectivity term
       
    TarFSrcTar = SteklovWeight * Tar_evecs_bound{i}(:, 1:ZOsize)*Tar_evals(1:ZOsize,1:ZOsize,i)*FmapSrcTar{i} +... %Steklov term
           ProperWeight * Tar_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar{i}+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Tar_evecs_bound{i}(:,1:ZOsize)+...%Orthonormality term.
           BijectWeight * Tar_evecs_bound{i}(:,1:ZOsize); %Bijectivity term
    

    p2pTarSrc{i} = annquery(SrcFTarSrc', TarFSrcTar', 1);
    
    
    
    %% Computing the p2pSrcTar map.

    SrcFSrcTar = SteklovWeight * Src_evecs_bound{i}(:, 1:ZOsize)*Src_evals(1:ZOsize,1:ZOsize,i)*FmapTarSrc{i} +... %Steklov term
           ProperWeight * Src_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc{i}+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Src_evecs_bound{i}(:,1:ZOsize)+... %Orthonormality term.
           BijectWeight * Src_evecs_bound{i}(:,1:ZOsize); %Bijectivity term
       
    TarFTarSrc = SteklovWeight * Tar_evecs_bound{i}(:, 1:ZOsize)*Tar_evals(1:ZOsize,1:ZOsize,i) +... %Steklov term
           ProperWeight * Tar_evecs_bound{i}(:,1:ZOsize)+... %"Proper" functional map term (pullback of p2p map)
           OrthoWeight * Tar_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc{i}'+...%Orthonormality term.
           BijectWeight * Tar_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar{i}; %Bijectivity term
       
 
    p2pSrcTar{i} = annquery(TarFTarSrc', SrcFSrcTar', 1);
 
    
end




end