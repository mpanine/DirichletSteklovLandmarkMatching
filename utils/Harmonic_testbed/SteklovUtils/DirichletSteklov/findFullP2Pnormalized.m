function [fullp2pTarSrc, fullp2pSrcTar] = findFullP2Pnormalized(Src_eF, Src_evals, Tar_eF, Tar_evals, FmapSrcTar, FmapTarSrc)
%IMPLEMENTS AN EXPERIMENTAL AND POORLY JUSTIFIED NORMALIZATION BEFORE
%EXTRACTING THE P2P MAP.

%Converts the functional maps obtained at the boundary into point to point
%maps in the bulk.

ep = 1e-9;

num_bounds = size(Src_eF,3);
num_eigs = size(Src_eF,2);


Src_evecs_combined = zeros(size(Src_eF,1), size(Src_eF,2)*size(Src_eF,3));
Tar_evecs_combined = zeros(size(Tar_eF,1), size(Tar_eF,2)*size(Tar_eF,3));

Src_evals_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);
Tar_evals_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);

FmapSrcTar_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);
FmapTarSrc_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);

% Src_normalization = spalloc(num_bounds*size(Src_eF,1), num_bounds*size(Src_eF,1), num_bounds*size(Src_eF,1)^2);
% Tar_normalization = spalloc(num_bounds*size(Tar_eF,1), num_bounds*size(Tar_eF,1), num_bounds*size(Tar_eF,1)^2);



% TEST THIS IF ABOVE IS SLOW
% Src_evals_combined = spalloc(num_bounds*num_eigs, num_bounds*num_eigs, num_bounds*num_eigs);
% Tar_evals_combined = spalloc(num_bounds*num_eigs, num_bounds*num_eigs, num_bounds*num_eigs);
% 
% FmapSrcTar_combined = spalloc(num_bounds*num_eigs, num_bounds*num_eigs, num_bounds*num_eigs*num_bounds*num_eigs);
% FmapTarSrc_combined = spalloc(num_bounds*num_eigs, num_bounds*num_eigs, num_bounds*num_eigs*num_bounds*num_eigs);




for i = 1:num_bounds
        
       start_ind = (i-1)*num_eigs + 1;
       end_ind = i*num_eigs;
       

       
       Src_evecs_combined(:, start_ind:end_ind ) = Src_eF(:,:,i)./(Src_eF(:,1,i) + ep);
       Tar_evecs_combined(:, start_ind:end_ind ) = Tar_eF(:,:,i)./(Tar_eF(:,1,i) + ep);
       

       
       Src_evals_combined(start_ind:end_ind, start_ind:end_ind) = Src_evals(:,:,i);
       Tar_evals_combined(start_ind:end_ind, start_ind:end_ind) = Tar_evals(:,:,i);
    
       FmapSrcTar_combined(start_ind:end_ind, start_ind:end_ind) = FmapSrcTar{i};
       FmapTarSrc_combined(start_ind:end_ind, start_ind:end_ind) = FmapTarSrc{i};
       
       
       
       Src_start = (i-1)*size(Src_eF,1) + 1;
       Src_end = i*size(Src_eF,1);
       
       Tar_start = (i-1)*size(Tar_eF,1) + 1;
       Tar_end = i*size(Tar_eF,1);
       
%        Src_normalization(Src_start:Src_end, Src_start:Src_end) = diag( (Src_eF(:,1,i) + ep)^-1);
%        Tar_normalization(Tar_start:Tar_end, Tar_start:Tar_end) = diag( (Tar_eF(:,1,i) + ep)^-1);
       
end


    
%% Computing the fullp2pTarSrc map.

    SrcFTarSrc = Src_evecs_combined;
       
    TarFSrcTar = Tar_evecs_combined*FmapSrcTar_combined;
       
       
%     SrcFTarSrc = Src_evecs_combined*Src_evals_combined +... %Steklov term
%            Src_evecs_combined+... %"Proper" functional map term (pullback of p2p map)
%            Src_evecs_combined*FmapSrcTar_combined'+... %Orthonormality term.
%            Src_evecs_combined*FmapTarSrc_combined; %Bijectivity term
%        
%     TarFSrcTar = Tar_evecs_combined*Tar_evals_combined*FmapSrcTar_combined +... %Steklov term
%            Tar_evecs_combined*FmapSrcTar_combined+... %"Proper" functional map term (pullback of p2p map)
%            Tar_evecs_combined+...%Orthonormality term.
%            Tar_evecs_combined; %Bijectivity term
    
    fullp2pTarSrc = annquery(SrcFTarSrc', TarFSrcTar', 1);
    
    
    
    %% Computing the fullp2pSrcTar map.

    SrcFSrcTar = Src_evecs_combined*FmapTarSrc_combined;
       
    TarFTarSrc = Tar_evecs_combined;
       
       
%            SrcFSrcTar = Src_evecs_combined*Src_evals_combined*FmapTarSrc_combined +... %Steklov term
%            Src_evecs_combined*FmapTarSrc_combined+... %"Proper" functional map term (pullback of p2p map)
%            Src_evecs_combined+... %Orthonormality term.
%            Src_evecs_combined; %Bijectivity term
%        
%     TarFTarSrc = Tar_evecs_combined*Tar_evals_combined +... %Steklov term
%            Tar_evecs_combined+... %"Proper" functional map term (pullback of p2p map)
%            Tar_evecs_combined*FmapTarSrc_combined'+...%Orthonormality term.
%            Tar_evecs_combined*FmapSrcTar_combined; %Bijectivity term
    
    
    fullp2pSrcTar = annquery(TarFTarSrc', SrcFSrcTar', 1);






















end