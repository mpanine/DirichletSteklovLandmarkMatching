function [FmapSrcTar, p2pTarSrc, FmapTarSrc, p2pSrcTar] = DirichletSteklovZoomOutV4(Src_refined, Src_evecs_bound, Src_evals, Src_Mass_Boundary, Src_refined_boundaries, Src_normalTriangles, Src_harmonic_basis, Tar_refined, Tar_evecs_bound, Tar_evals, Tar_Mass_Boundary, Tar_refined_boundaries, Tar_normalTriangles, Tar_harmonic_basis)
%SteklovZoomOutV3: Third version of the Dirichlet-Steklov ZoomOut:
%Implements Dirichlet-to-Neumann commutativity, orthogonality and
%bijectivity

%THIS VERSION IMPLEMENTS AN EXHAUSTIVE SEARCH FOR THE OPTIMAL INITIAL
%GUESS BASED ON NORMAL DERIVATIVE MATCHING. THEN IT PERFORMS A FEW ICP-LIKE ITERATIONS TO SMOOTH EVERYTHING OUT.


toplot = false; %Set to true for plots.


num_shifts = 200; % Determines the number of different shifts tried. Includes the trivial shift
num_it_ICP = 10;

%% Preliminaries

num_bounds = length(Src_evecs_bound);
num_eigs = size(Src_evals,2);

%% Optimize the initial correspondence via a shifting of the orientation

best_shifts = zeros(num_bounds,1); %Stores the answer of the first optimization.

for b = 1:num_bounds % Iterate over each boundary
    
    
    
    Src_nv = size(Src_evecs_bound{b},1);
    Tar_nv = size(Tar_evecs_bound{b},1);
        
    Src_ND = zeros(Src_nv,num_bounds-1);
    Tar_ND = zeros(Tar_nv,num_bounds-1);
    
    inds = [1:(b-1) (b+1):num_bounds];
    
    for i = 1:(num_bounds-1)
        
%         Src_ND(:,i) = normalDerivativeAtLandmarkBoundary(Src_refined.SHAPE.surface, Src_harmonic_basis(:,inds(i)), Src_refined_boundaries{b}, Src_Mass_Boundary{b}, Src_landmarks(b), Src.SHAPE.nf);
%         Tar_ND(:,i) = normalDerivativeAtLandmarkBoundary(Tar_refined.SHAPE.surface, Tar_harmonic_basis(:,inds(i)), Tar_refined_boundaries{b}, Tar_Mass_Boundary{b}, Tar_landmarks(b), Tar.SHAPE.nf);
        
        Src_ND(:,i) = normalDerivativeAtLandmarkBoundaryBoundaryTriangles(Src_refined.SHAPE.surface, Src_normalTriangles{b}, Src_harmonic_basis(:,inds(i)), Src_refined_boundaries{b}, Src_Mass_Boundary{b});
        Tar_ND(:,i) = normalDerivativeAtLandmarkBoundaryBoundaryTriangles(Tar_refined.SHAPE.surface, Tar_normalTriangles{b}, Tar_harmonic_basis(:,inds(i)), Tar_refined_boundaries{b}, Tar_Mass_Boundary{b});
        
    end
    
    
    
    Tar_line = (1:Tar_nv)/Tar_nv;
    
    bestE = Inf;
    
    for s = 0:(num_shifts - 1)
        
        Src_line = mod( (1:Src_nv)/Src_nv + s/num_shifts , 1);
        
        p2pTarSrc_pre = annquery(Src_line, Tar_line,1);
        p2pSrcTar_pre = annquery(Tar_line, Src_line,1);
        

%         FmapSrcTar_pre = Tar_evecs_bound{b}' * Tar_Mass_Boundary{b} * Src_evecs_bound{b}(p2pTarSrc_pre,:);
%         FmapTarSrc_pre = Src_evecs_bound{b}' * Src_Mass_Boundary{b} * Tar_evecs_bound{b}(p2pSrcTar_pre,:);
%      
%         E = DS_Fmap_energy(Src_evals(:,:,b), Tar_evals(:,:,b), FmapSrcTar_pre, FmapTarSrc_pre) + ...
%             DS_Fmap_energy(Tar_evals(:,:,b), Src_evals(:,:,b), FmapTarSrc_pre, FmapSrcTar_pre);

        Src_on_Tar = Src_ND(p2pTarSrc_pre,:);
        Tar_on_Src = Tar_ND(p2pSrcTar_pre,:);
        
        E = (Tar_ND - Src_on_Tar)' * Tar_Mass_Boundary{b} * (Tar_ND - Src_on_Tar) +...
            (Src_ND - Tar_on_Src)' * Src_Mass_Boundary{b} * (Src_ND - Tar_on_Src);       
        
        
        if E < bestE
            
            bestE = E;            
            best_shifts(b) = s; 
%             current_best = s
            
            p2pTarSrc{b} = p2pTarSrc_pre;
            p2pSrcTar{b} = p2pSrcTar_pre;
            
        end


    end
    
   
    
    
end


 






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












for it = 1:num_it_ICP
    
    FmapSrcTar = SteklovDirichletFmapFromP2P(Src_evecs_bound, Tar_evecs_bound, Tar_Mass_Boundary, p2pTarSrc, num_eigs);
    FmapTarSrc = SteklovDirichletFmapFromP2P(Tar_evecs_bound, Src_evecs_bound, Src_Mass_Boundary, p2pSrcTar, num_eigs);
%     p2pTarSrc = DirichletSteklovFmapToP2P(Src_evecs_bound, Src_evals_mod, Tar_evecs_bound, Tar_evals_mod, FmapSrcTar, it);    
    
    [p2pTarSrc, p2pSrcTar] = DirichletSteklovFmapToP2Pv2(Src_evecs_bound, Src_evals, Tar_evecs_bound, Tar_evals, FmapSrcTar, FmapTarSrc, num_eigs);
    
    
    
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