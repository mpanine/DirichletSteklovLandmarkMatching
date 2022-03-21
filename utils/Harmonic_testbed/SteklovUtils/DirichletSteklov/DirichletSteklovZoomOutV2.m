function [FmapSrcTar, p2pTarSrc, FmapTarSrc, p2pSrcTar] = DirichletSteklovZoomOutV2(Src_evecs_bound, Src_evals, Src_Mass_Boundary, Tar_evecs_bound, Tar_evals, Tar_Mass_Boundary)
%SteklovZoomOutV1: Second version of the Dirichlet-Steklov ZoomOut:
%Implements Dirichlet-to-Neumann commutativity, orthogonality and
%bijectivity

toplot = true; %Set to true for plots.

ZO_start = 3;

%evecs_bound: eigenvectors on the boundary (cell array over boundary components)
%evals: eigenvalues as a 3d array (diagonal matrices for each boundary component)

% max_size = 20;

Src_evals_mod = Src_evals;
Tar_evals_mod = Tar_evals;


%     for i = 1:length(Src_evecs_bound) % Optional resolvent formulation. DOES NOT WORK RIGHT NOW.
%         Src_evals_mod(:,:,i) = diag( diag(Src_evals(:,:,i) + 10).^-1  );
%         Tar_evals_mod(:,:,i) = diag( diag(Tar_evals(:,:,i) + 10).^-1  );
%     end

%% Guess initial map. Version 1.
% 
% for i = 1:length(Src_evecs_bound) %Guess the initial Fmap
% 
%     Src_guess{i} = sum(Src_Mass_Boundary{i}*Src_evecs_bound{i},1); %Unit function on relevant boundary expanded in basis.
%     Tar_guess{i} = sum(Tar_Mass_Boundary{i}*Tar_evecs_bound{i},1);
%     Fmap{i} = Tar_guess{i}(1,1)' * Src_guess{i}(1,1);
% 
% end
% 
% p2p = DirichletSteklovFmapToP2P(Src_evecs_bound, Src_evals_mod, Tar_evecs_bound, Tar_evals_mod, Fmap, 1);


%% Guess initial map. Version 2: assume circular-ish correspondence

for i = 1:1:length(Src_evecs_bound)
    
    Src_nv = size(Src_evecs_bound{i},1);
    Tar_nv = size(Tar_evecs_bound{i},1);
    
    Src_line = mod( (1:Src_nv)/Src_nv + .25 , 1);
    Tar_line = (1:Tar_nv)/Tar_nv;
    p2pTarSrc{i} = annquery(Src_line,Tar_line,1);
    p2pSrcTar{i} = annquery(Tar_line, Src_line,1);

end

%%

if toplot
    figure('Units','normalized','Position',[0.25 0.25 0.5 0.5])
    bc_to_plot = 6;
    

    clf
    subplot(2,3,1)
    imagesc(zeros(size(Src_evals,1))); colorbar; title('Fmap Src->Tar, It: 0')
    subplot(2,3,2)
    imagesc(zeros(size(Src_evals,1))); colorbar; title('Fmap Tar->Src, It: 0')
    subplot(2,3,3)
    imagesc(zeros(size(Src_evals,1))); colorbar; title('Inverse Check, It: 0')
    pause(1e-9)
    
    subplot(2,3,4)
    plot(p2pTarSrc{bc_to_plot},'.'); title(sprintf('p2p Tar->Src, It: %d',0))
    subplot(2,3,5)
    plot(p2pSrcTar{bc_to_plot},'.'); title(sprintf('p2p Src->Tar, It: %d',0))
    subplot(2,3,6)
    plot(p2pTarSrc{bc_to_plot}(p2pSrcTar{bc_to_plot}),'.'); title(sprintf('p2p Tar->Tar inv check, It: %d',0))
    pause(2)
    
    

    
    
end












for it = ZO_start:size(Src_evals_mod,2)
    
    FmapSrcTar = SteklovDirichletFmapFromP2P(Src_evecs_bound, Tar_evecs_bound, Tar_Mass_Boundary, p2pTarSrc, it);
    FmapTarSrc = SteklovDirichletFmapFromP2P(Tar_evecs_bound, Src_evecs_bound, Src_Mass_Boundary, p2pSrcTar, it);
%     p2pTarSrc = DirichletSteklovFmapToP2P(Src_evecs_bound, Src_evals_mod, Tar_evecs_bound, Tar_evals_mod, FmapSrcTar, it);    
    
    [p2pTarSrc, p2pSrcTar] = DirichletSteklovFmapToP2Pv2(Src_evecs_bound, Src_evals, Tar_evecs_bound, Tar_evals, FmapSrcTar, FmapTarSrc, it);
    
    
    
    if toplot
        clf
        subplot(2,3,1)
        imagesc(FmapSrcTar{bc_to_plot}); colorbar; title(sprintf('Fmap Src->Tar, It: %d',it))
        subplot(2,3,2)
        imagesc(FmapTarSrc{bc_to_plot}); colorbar; title(sprintf('Fmap Tar->Src, It: %d',it))
        subplot(2,3,3)
        imagesc(FmapSrcTar{bc_to_plot}*FmapTarSrc{bc_to_plot}); colorbar; title(sprintf('Inverse Check, It: %d',it))
        pause(1e-9)
        
        subplot(2,3,4)
        plot(p2pTarSrc{bc_to_plot},'.'); title(sprintf('p2p Tar->Src, It: %d',it))
        subplot(2,3,5)
        plot(p2pSrcTar{bc_to_plot},'.'); title(sprintf('p2p Src->Tar, It: %d',it))
        subplot(2,3,6)
        plot(p2pTarSrc{bc_to_plot}(p2pSrcTar{bc_to_plot}),'.'); title(sprintf('p2p Tar->Tar inv check, It: %d',it))
        pause(.5)
        
        
    end
    
end






















end