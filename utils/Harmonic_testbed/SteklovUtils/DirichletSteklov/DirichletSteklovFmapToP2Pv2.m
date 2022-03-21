function [p2pTarSrc, p2pSrcTar] = DirichletSteklovFmapToP2Pv2(Src_evecs_bound, Src_evals, Tar_evecs_bound, Tar_evals, FmapSrcTar, FmapTarSrc, ZOsize)
%Implements bijectivity condition as well as the others

%Fmap: functional map Src -> Tar 
%p2p: point-to-point map Tar -> Src

for i = 1:length(Src_evecs_bound)
    
%     size(Fmap{i})
    
%% Computing the p2pTarSrc map.

    SrcFTarSrc = Src_evecs_bound{i}(:, 1:ZOsize)*Src_evals(1:ZOsize,1:ZOsize,i) +... %Steklov term
           Src_evecs_bound{i}(:,1:ZOsize)+... %"Proper" functional map term (pullback of p2p map)
           Src_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar{i}'+... %Orthonormality term.
           Src_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc{i}; %Bijectivity term
       
    TarFSrcTar = Tar_evecs_bound{i}(:, 1:ZOsize)*Tar_evals(1:ZOsize,1:ZOsize,i)*FmapSrcTar{i} +... %Steklov term
           Tar_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar{i}+... %"Proper" functional map term (pullback of p2p map)
           Tar_evecs_bound{i}(:,1:ZOsize)+...%Orthonormality term.
           Tar_evecs_bound{i}(:,1:ZOsize); %Bijectivity term
    
%            SrcFTarSrc = Src_evecs_bound{i}(:,1:ZOsize)+... %"Proper" functional map term (pullback of p2p map)
%            Src_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar{i}'+... %Orthonormality term.
%            Src_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc{i}; %Bijectivity term
%        
%     TarFSrcTar = Tar_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar{i}+... %"Proper" functional map term (pullback of p2p map)
%            Tar_evecs_bound{i}(:,1:ZOsize)+...%Orthonormality term.
%            Tar_evecs_bound{i}(:,1:ZOsize); %Bijectivity term
       
    
    p2pTarSrc{i} = annquery(SrcFTarSrc', TarFSrcTar', 1);
    
    
    
    %% Computing the p2pSrcTar map.

    SrcFSrcTar = Src_evecs_bound{i}(:, 1:ZOsize)*Src_evals(1:ZOsize,1:ZOsize,i)*FmapTarSrc{i} +... %Steklov term
           Src_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc{i}+... %"Proper" functional map term (pullback of p2p map)
           Src_evecs_bound{i}(:,1:ZOsize)+... %Orthonormality term.
           Src_evecs_bound{i}(:,1:ZOsize); %Bijectivity term
       
    TarFTarSrc = Tar_evecs_bound{i}(:, 1:ZOsize)*Tar_evals(1:ZOsize,1:ZOsize,i) +... %Steklov term
           Tar_evecs_bound{i}(:,1:ZOsize)+... %"Proper" functional map term (pullback of p2p map)
           Tar_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc{i}'+...%Orthonormality term.
           Tar_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar{i}; %Bijectivity term
       
       
%            SrcFSrcTar = Src_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc{i}+... %"Proper" functional map term (pullback of p2p map)
%            Src_evecs_bound{i}(:,1:ZOsize); %+... %Orthonormality term.
%            Src_evecs_bound{i}(:,1:ZOsize); %Bijectivity term
       
%     TarFTarSrc = Tar_evecs_bound{i}(:,1:ZOsize)+... %"Proper" functional map term (pullback of p2p map)
%            Tar_evecs_bound{i}(:,1:ZOsize)*FmapTarSrc{i}';%+...%Orthonormality term.
           %Tar_evecs_bound{i}(:,1:ZOsize)*FmapSrcTar{i}; %Bijectivity term
    
    
    p2pSrcTar{i} = annquery(TarFTarSrc', SrcFSrcTar', 1);
 
    
end




end